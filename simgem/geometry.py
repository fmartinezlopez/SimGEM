"""
Unit-cell geometry and mesh (gmsh/OpenCASCADE).

The cell is a rectangle that is a translational period of every foil's hole
pattern (p x p for square, p x p*sqrt(3) for hex, or a common supercell), with
translational-periodic side walls. Each foil's holes are placed on its own
lattice, shifted by its offset; every lattice hole that touches the cell is
cut, so a hole straddling a wall reappears on the opposite wall and the two
walls are geometrically identical (needed for a conformal periodic mesh).

Copper layers are real volumes, cut out of the domain: their surfaces become
the electrode boundaries and electrons that hit them leave the field map.
With copper_thickness_um = 0 the electrodes are the bare kapton faces.

Lengths in mm (Elmer/Garfield are told the mesh is in mm).
"""
import math
import os

from .config import lattice

EPS = 1e-6           # geometric matching tolerance (mm)
CLEARANCE = 1e-3     # min distance between a hole edge and a wall/corner (mm)


def wrap(v, period):
    """Fold into [-period/2, period/2)."""
    return v - period * math.floor(v / period + 0.5)


def wrap_offset(foil, off):
    """Offsets are only defined modulo the foil's own pattern period."""
    _, _, cx, cy = lattice(foil)
    return (wrap(off[0], cx), wrap(off[1], cy))


def hole_radii(foil):
    """Radii (mm) of the kapton hole at bottom/middle/top and of the copper
    holes at bottom/top."""
    h = foil["hole"]
    rO = h["outer_diameter_um"] / 2000.0
    rI = h["inner_diameter_um"] / 2000.0
    prof = h["profile"]
    if prof == "double_conical":
        kb, km, kt = rO, rI, rO
    elif prof == "cylindrical":
        kb = km = kt = rO
    else:  # conical
        kb, kt = (rI, rO) if h["wide_side"] == "top" else (rO, rI)
        km = (kb + kt) / 2
    if h.get("copper_diameter_um"):
        cb = ct = h["copper_diameter_um"] / 2000.0
    else:
        cb, ct = kb, kt
    return {"kapton": (kb, km, kt), "copper": (cb, ct)}


def hole_centres(foil, offset, cell):
    """Centres (mm) of every lattice hole of `foil` that touches the cell."""
    a1, a2, _, _ = lattice(foil)
    Lx, Ly = cell
    rad = hole_radii(foil)
    R = max(max(rad["kapton"]), max(rad["copper"]))
    p = foil["pitch_um"] / 1000.0
    n = int(math.ceil(max(Lx, Ly) / p)) + 3
    out = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            x = offset[0] + i * a1[0] + j * a2[0]
            y = offset[1] + i * a1[1] + j * a2[1]
            if abs(x) < Lx / 2 + R and abs(y) < Ly / 2 + R:
                out.append((x, y))
    return out


def _too_close(foil, offset, cell):
    """True if a hole edge nearly touches a wall or corner (would create
    slivers / a degenerate periodic face)."""
    Lx, Ly = cell
    rad = hole_radii(foil)
    radii = set(rad["kapton"]) | set(rad["copper"])
    for x, y in hole_centres(foil, offset, cell):
        for r in radii:
            for w, c in ((Lx / 2, x), (Lx / 2, -x), (Ly / 2, y), (Ly / 2, -y)):
                if abs((w - c) - r) < CLEARANCE:
                    return True
            for sx in (-1, 1):
                for sy in (-1, 1):
                    d = math.hypot(x - sx * Lx / 2, y - sy * Ly / 2)
                    if d > EPS and abs(d - r) < CLEARANCE:
                        return True
    return False


def safe_offsets(stack, log=print):
    """Per-foil offsets wrapped into the pattern period and nudged (by <= a
    few um) away from hole/wall tangencies that gmsh cannot mesh cleanly."""
    out = {}
    for f in stack.foils:
        off = wrap_offset(f, stack.offsets[f["name"]])
        tries = 0
        while _too_close(f, off, stack.cell):
            tries += 1
            if tries > 20:
                raise RuntimeError(f"{f['name']}: could not find a clean hole position near offset {off}")
            off = wrap_offset(f, (off[0] + 0.7e-3, off[1] + 0.4e-3))
        if tries:
            log(f"  {f['name']}: offset nudged by {tries} step(s) to avoid a hole/wall tangency "
                f"-> ({off[0] * 1e3:+.2f}, {off[1] * 1e3:+.2f}) um")
        out[f["name"]] = off
    return out


def _kapton_hole(occ, foil, x, y, z0, z1):
    kb, km, kt = hole_radii(foil)["kapton"]
    t = z1 - z0
    prof = foil["hole"]["profile"]
    if prof == "cylindrical":
        return [(3, occ.addCylinder(x, y, z0 - EPS, 0, 0, t + 2 * EPS, kb))]
    if prof == "conical":
        return [(3, occ.addCone(x, y, z0, 0, 0, t, kb, kt))]
    c1 = occ.addCone(x, y, z0, 0, 0, t / 2, kb, km)
    c2 = occ.addCone(x, y, z0 + t / 2, 0, 0, t / 2, km, kt)
    fused, _ = occ.fuse([(3, c1)], [(3, c2)])
    return fused


def build_mesh(stack, offsets, outdir, show_gui=False, log=print):
    """Build the cell, mesh it and write <outdir>/stack.msh.
    Returns a dict describing the mesh (physical names -> tags, stats)."""
    import gmsh

    Lx, Ly = stack.cell
    H = stack.height
    m = stack.cfg["mesh"]
    rmin = min(hole_radii(f)["kapton"][1] for f in stack.foils)
    lc_min, lc_max = m["fine_um"] / 1000.0, m["coarse_um"] / 1000.0
    d_min = (m["fine_distance_um"] / 1000.0) if m["fine_distance_um"] else 2 * rmin
    d_max = (m["coarse_distance_um"] / 1000.0) if m["coarse_distance_um"] else 6 * rmin

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("stack")
    occ = gmsh.model.occ
    box = occ.addBox(-Lx / 2, -Ly / 2, 0, Lx, Ly, H)

    parts = []  # (kind, foil name, side, dimtags)
    for f in stack.foils:
        name = f["name"]
        zk0, zk1 = stack.kapton[name]
        zb, zt = stack.z[name + "_bot"], stack.z[name + "_top"]
        centres = hole_centres(f, offsets[name], stack.cell)
        slab = occ.addBox(-Lx / 2, -Ly / 2, zk0, Lx, Ly, zk1 - zk0)
        tools = []
        for x, y in centres:
            tools += _kapton_hole(occ, f, x, y, zk0, zk1)
        kap, _ = occ.cut([(3, slab)], tools)
        parts.append(("kapton", name, None, kap))
        if f["copper_thickness_um"] > 0:
            cb, ct = hole_radii(f)["copper"]
            for side, z0, z1, r in (("bot", zb, zk0, cb), ("top", zk1, zt, ct)):
                cu = occ.addBox(-Lx / 2, -Ly / 2, z0, Lx, Ly, z1 - z0)
                holes = [(3, occ.addCylinder(x, y, z0 - EPS, 0, 0, z1 - z0 + 2 * EPS, r))
                         for x, y in centres]
                cut, _ = occ.cut([(3, cu)], holes)
                parts.append(("copper", name, side, cut))

    tools = [dt for p in parts for dt in p[3]]
    _, out_map = occ.fragment([(3, box)], tools)
    occ.synchronize()

    # out_map[0] -> pieces of the box; out_map[1 + k] -> pieces of tools[k]
    k = 1
    kapton_vols, copper_vols = {}, []
    for kind, name, side, dts in parts:
        vols = []
        for _ in dts:
            vols += [t for d, t in out_map[k] if d == 3]
            k += 1
        if kind == "kapton":
            kapton_vols[name] = vols
        else:
            copper_vols.append((name, side, vols))
    solid = {t for v in kapton_vols.values() for t in v} | {t for _, _, v in copper_vols for t in v}
    gas_vols = [t for d, t in out_map[0] if d == 3 and t not in solid]

    def on_wall(s):
        b = gmsh.model.getBoundingBox(2, s)
        return ((abs(b[0] - b[3]) < EPS and abs(abs(b[0]) - Lx / 2) < EPS) or
                (abs(b[1] - b[4]) < EPS and abs(abs(b[1]) - Ly / 2) < EPS))

    def plane(zv):
        return [t for _, t in gmsh.model.getEntitiesInBoundingBox(
            -Lx, -Ly, zv - EPS, Lx, Ly, zv + EPS, 2)]

    electrodes = {"drift": plane(H), "readout": plane(0.0)}
    for f in stack.foils:
        name = f["name"]
        if f["copper_thickness_um"] > 0:
            for n, side, vols in copper_vols:
                if n != name:
                    continue
                surfs = {abs(t) for _, t in gmsh.model.getBoundary(
                    [(3, v) for v in vols], combined=False, oriented=False)}
                electrodes[f"{name}_{side}"] = sorted(s for s in surfs if not on_wall(s))
        else:
            electrodes[f"{name}_bot"] = plane(stack.z[name + "_bot"])
            electrodes[f"{name}_top"] = plane(stack.z[name + "_top"])

    # copper is a conductor: remove its volumes (keep their surfaces) and any
    # surface left bounding nothing (the copper faces on the side walls)
    if copper_vols:
        occ.remove([(3, v) for _, _, vols in copper_vols for v in vols], recursive=False)
        occ.synchronize()
        orphans = [(2, s) for _, s in gmsh.model.getEntities(2)
                   if len(gmsh.model.getAdjacencies(2, s)[0]) == 0]
        if orphans:
            occ.remove(orphans, recursive=False)
            occ.synchronize()

    surf_all = [t for _, t in gmsh.model.getEntities(2)]
    for n, s in electrodes.items():
        missing = [t for t in s if t not in surf_all]
        if not s or missing:
            raise RuntimeError(f"electrode '{n}': no surfaces found (geometry construction failed)")

    # ---- physical groups -------------------------------------------------
    phys = {}
    gmsh.model.addPhysicalGroup(3, gas_vols, tag=1, name="gas")
    phys["gas"] = 1
    for i, f in enumerate(stack.foils):
        nm = f"kapton_{f['name']}"
        gmsh.model.addPhysicalGroup(3, kapton_vols[f["name"]], tag=2 + i, name=nm)
        phys[nm] = 2 + i
    tag = 100
    for n in stack.electrodes():
        gmsh.model.addPhysicalGroup(2, electrodes[n], tag=tag, name=n)
        phys[n] = tag
        tag += 1

    # ---- periodic side walls --------------------------------------------
    def wall(axis, val):
        if axis == 0:
            box_ = (val - EPS, -Ly, -EPS, val + EPS, Ly, H + EPS)
        else:
            box_ = (-Lx, val - EPS, -EPS, Lx, val + EPS, H + EPS)
        return [t for _, t in gmsh.model.getEntitiesInBoundingBox(*box_, 2)]

    for axis, L, lo_name, hi_name in ((0, Lx, "x_lo", "x_hi"), (1, Ly, "y_lo", "y_hi")):
        lo, hi = wall(axis, -L / 2), wall(axis, L / 2)
        pairs = []
        for s in hi:
            bh = gmsh.model.getBoundingBox(2, s)
            match = None
            for t in lo:
                bl = list(gmsh.model.getBoundingBox(2, t))
                bl[axis] += L
                bl[axis + 3] += L
                if all(abs(a - b) < 10 * EPS for a, b in zip(bh, bl)):
                    match = t
                    break
            if match is None:
                raise RuntimeError(f"periodic wall {hi_name}: surface {s} has no partner "
                                   f"on {lo_name} (geometry not periodic)")
            pairs.append((s, match))
        if len(pairs) != len(lo):
            raise RuntimeError(f"periodic walls {lo_name}/{hi_name} have different surfaces")
        T = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        T[3 if axis == 0 else 7] = L
        gmsh.model.mesh.setPeriodic(2, [p[0] for p in pairs], [p[1] for p in pairs], T)
        gmsh.model.addPhysicalGroup(2, lo, tag=tag, name=lo_name)
        phys[lo_name] = tag
        gmsh.model.addPhysicalGroup(2, hi, tag=tag + 1, name=hi_name)
        phys[hi_name] = tag + 1
        tag += 2

    # ---- mesh size: fine near the foils, coarse in the gaps --------------
    refine = set()
    for vols in kapton_vols.values():
        refine |= {abs(t) for _, t in gmsh.model.getBoundary(
            [(3, v) for v in vols], combined=False, oriented=False)}
    for n, s in electrodes.items():
        if n not in ("drift", "readout"):
            refine |= set(s)
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "SurfacesList", sorted(refine))
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", lc_min)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", lc_max)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", d_min)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", d_max)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)
    for opt in ("MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"):
        gmsh.option.setNumber("Mesh." + opt, 0)

    gmsh.option.setNumber("Mesh.ElementOrder", 2)       # Garfield wants 10-node tets
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # what ElmerGrid reads
    gmsh.model.mesh.generate(3)

    os.makedirs(outdir, exist_ok=True)
    msh = os.path.join(outdir, "stack.msh")
    gmsh.write(msh)
    _, etags, _ = gmsh.model.mesh.getElements(3)
    n_el = sum(len(t) for t in etags)
    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    log(f"  mesh: {n_el} tetrahedra, {n_nodes} nodes -> {msh}")
    if show_gui:
        gmsh.fltk.run()
    gmsh.finalize()
    return {"physical": phys, "elements": n_el, "nodes": n_nodes}
