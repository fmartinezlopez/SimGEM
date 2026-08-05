"""
Geometry + field-map construction for the parametric triple-GEM stack.
Imported by gem.py; not run directly.

Modeling choices (see README):
  * holes aligned across the three foils, square cell, mirror-periodic walls
  * copper = zero-thickness electrode on the kapton top/bottom faces
  * drift gap truncated to `drift_gap_modeled_mm`, carrying the true drift field
Lengths inside gmsh are in mm (ComponentElmer is told "mm").
"""
import json, os, re, subprocess, sys, contextlib, datetime, random, math
from pathlib import Path


# ----------------------------------------------------------------- logging utils
class _Tee:
    """Write to several streams at once (console + a log file)."""
    def __init__(self, *streams):
        self._s = streams

    def write(self, x):
        for s in self._s:
            s.write(x); s.flush()

    def flush(self):
        for s in self._s:
            s.flush()

    def fileno(self):
        return self._s[0].fileno()


@contextlib.contextmanager
def stage_log(rundir, name):
    """Tee all stdout/stderr of a stage into rundir/<name>.log."""
    os.makedirs(rundir, exist_ok=True)
    path = os.path.join(rundir, f"{name}.log")
    f = open(path, "w")
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = _Tee(old_out, f), _Tee(old_err, f)
    print(f"===== {name}  {datetime.datetime.now().isoformat(timespec='seconds')} =====")
    try:
        yield path
    finally:
        sys.stdout, sys.stderr = old_out, old_err
        f.close()
        old_out.write(f"[log -> {path}]\n")


def run_streamed(cmd, cwd=None):
    """Run a subprocess, streaming its output through the current stdout
    (so it lands in both the console and the active stage log)."""
    cmd = [str(c) for c in cmd]
    print("+ " + " ".join(cmd))
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in p.stdout:
        sys.stdout.write(line)
    p.wait()
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd)


# ----------------------------------------------------------------- config helpers
def load_config(path):
    with open(path) as f:
        return json.load(f)


def run_dir(config_path):
    """Every config works in its own folder: runs/<config stem>/."""
    return os.path.join("runs", Path(config_path).stem)


def gas_path(cfg):
    """Gas-file path from the config, or auto-named from the gas parameters."""
    gp = cfg["gas"].get("gas_file")
    if gp:
        return gp
    g = cfg["gas"]
    name = f"ar{int(g['ar_fraction'])}_ch4{int(g['ch4_fraction'])}_{int(g['pressure_bar'])}bar.gas"
    return os.path.join("gas", name)


# ----------------------------------------------------------------- misalignment
# A `misalignment` block in the config switches the whole run from the legacy
# aligned MIRROR cell to a translational-PERIODIC cell whose three foils can be
# shifted laterally. With NO such block the geometry, BCs and Garfield
# periodicity are byte-for-byte the legacy behaviour (mirror walls, no periodic
# BC). Offsets are per-foil (dx, dy); only the RELATIVE offsets affect the gain,
# but per-foil is the natural way to express both fixed and random scenarios.
FOILS = ("gem1", "gem2", "gem3")


def has_misalignment(cfg):
    """True if this config uses the periodic (misalignment-capable) path."""
    return cfg.get("misalignment") is not None


def _pitch_mm(cfg):
    return cfg["geometry"]["pitch_um"] / 1000.0


def _unit_to_mm(value, unit, pitch_mm):
    u = (unit or "pitch_fraction").lower()
    if u in ("pitch_fraction", "pitch", "frac", "fraction"):
        return value * pitch_mm
    if u in ("um", "micron", "microns", "µm"):
        return value / 1000.0
    if u == "mm":
        return float(value)
    raise ValueError(f"unknown misalignment unit '{unit}' "
                     "(use 'pitch_fraction', 'um' or 'mm')")


def _wrap_mm(d, pitch_mm):
    """Fold an offset into [-pitch/2, pitch/2]; the infinite pattern is identical
    under a full-pitch shift, and the 9-image hole cut only covers one pitch."""
    return d - pitch_mm * round(d / pitch_mm)


def resolve_offsets(cfg, index=0):
    """Return (offsets_mm, meta) for one geometry.

      offsets_mm = {'gem1': (dx, dy), 'gem2': (...), 'gem3': (...)}  in mm
      meta       = small dict describing how they were drawn (for the manifest)

    `index` selects the draw in random mode and is deterministic given the
    config's seed, so a run is fully reproducible and re-runnable.
    """
    mis = cfg.get("misalignment") or {}
    mode = mis.get("mode", "aligned").lower()
    pitch = _pitch_mm(cfg)

    if mode == "aligned":
        return {f: (0.0, 0.0) for f in FOILS}, {"mode": "aligned", "index": index}

    if mode == "fixed":
        unit = mis.get("unit", "pitch_fraction")
        raw = mis.get("offsets", {})
        off = {}
        for f in FOILS:
            dx, dy = raw.get(f, [0.0, 0.0])
            off[f] = (_wrap_mm(_unit_to_mm(dx, unit, pitch), pitch),
                      _wrap_mm(_unit_to_mm(dy, unit, pitch), pitch))
        return off, {"mode": "fixed", "unit": unit, "index": index}

    if mode == "random":
        has_um = "sigma_um" in mis
        has_fr = "sigma_pitch_fraction" in mis
        if has_um == has_fr:
            raise ValueError("random misalignment needs EXACTLY ONE of "
                             "'sigma_um' or 'sigma_pitch_fraction'")
        sigma_mm = (mis["sigma_um"] / 1000.0) if has_um \
            else mis["sigma_pitch_fraction"] * pitch
        base = int(mis.get("seed", 12345))
        # deterministic, well-separated stream per config index
        rng = random.Random(base * 1_000_003 + index)
        off = {}
        for f in FOILS:                            # all three foils independent
            dx = _wrap_mm(rng.gauss(0.0, sigma_mm), pitch)
            dy = _wrap_mm(rng.gauss(0.0, sigma_mm), pitch)
            off[f] = (dx, dy)
        return off, {"mode": "random", "sigma_mm": sigma_mm,
                     "seed": base, "index": index}

    raise ValueError(f"unknown misalignment mode '{mode}' "
                     "(use 'aligned', 'fixed' or 'random')")


def write_periodicity(rundir, mode):
    """Marker read by gain_stack/vis_stack to pick mirror vs translational
    periodicity in Garfield. 'mirror' (legacy) or 'translate'."""
    with open(os.path.join(rundir, "periodicity"), "w") as f:
        f.write(mode + "\n")


def offsets_summary(offsets_mm, pitch_mm):
    """One-line human string of the per-foil offsets, in pitch fractions."""
    return "  ".join(
        f"{f}=({offsets_mm[f][0] / pitch_mm:+.3f},{offsets_mm[f][1] / pitch_mm:+.3f})p"
        for f in FOILS)


def compute_potentials(cfg):
    """Absolute electrode potentials (V). Readout = 0 (highest); potentials
    DECREASE upward so the field points up and electrons drift DOWN to the
    readout, as in a real detector (drift cathode most negative)."""
    st, fl = cfg["stack"], cfg["fields"]
    cm = lambda mm: mm / 10.0
    H_ind, H_t2 = cm(st["induction_gap_mm"]), cm(st["transfer_gap_2_mm"])
    H_t1, H_dr = cm(st["transfer_gap_1_mm"]), cm(st["drift_gap_modeled_mm"])
    E_ind = fl["induction_field_kVcm"] * 1e3
    E_t2 = fl["transfer_2_field_kVcm"] * 1e3
    E_t1 = fl["transfer_1_field_kVcm"] * 1e3
    E_dr = fl["drift_field_kVcm"] * 1e3
    V = {"readout": 0.0}
    V["gem3_bot"] = V["readout"] - E_ind * H_ind
    V["gem3_top"] = V["gem3_bot"] - fl["dV_gem_bot_V"]
    V["gem2_bot"] = V["gem3_top"] - E_t2 * H_t2
    V["gem2_top"] = V["gem2_bot"] - fl["dV_gem_mid_V"]
    V["gem1_bot"] = V["gem2_top"] - E_t1 * H_t1
    V["gem1_top"] = V["gem1_bot"] - fl["dV_gem_top_V"]
    V["drift"] = V["gem1_top"] - E_dr * H_dr
    return V


def z_landmarks(cfg):
    """z of every interface (mm), z=0 at readout, increasing upward."""
    tK = cfg["geometry"]["kapton_thickness_um"] / 1000.0
    st = cfg["stack"]
    z = {"readout": 0.0}
    z["gem3_bot"] = st["induction_gap_mm"]
    z["gem3_top"] = z["gem3_bot"] + tK
    z["gem2_bot"] = z["gem3_top"] + st["transfer_gap_2_mm"]
    z["gem2_top"] = z["gem2_bot"] + tK
    z["gem1_bot"] = z["gem2_top"] + st["transfer_gap_1_mm"]
    z["gem1_top"] = z["gem1_bot"] + tK
    z["drift"] = z["gem1_top"] + st["drift_gap_modeled_mm"]
    return z, tK


# ----------------------------------------------------------------- geometry / mesh
def build_mesh(cfg, z, tK, rundir, offsets_mm=None, periodic=False, show_gui=False):
    """Build the unit-cell mesh.

    periodic=False (default, legacy): aligned holes on the z-axis, mirror-
      symmetric side walls left untagged (Elmer's natural BC == mirror). Byte-
      for-byte the original behaviour; `offsets_mm` is ignored.
    periodic=True: each foil's hole is placed at its own (dx, dy) offset and the
      four side walls are tied together with a TRANSLATIONAL periodic mesh, so a
      rigid per-foil shift is represented exactly on a single pitch-sized cell.
    """
    import gmsh
    g = cfg["geometry"]
    p = g["pitch_um"] / 1000.0
    rO = g["hole_outer_diameter_um"] / 2000.0
    rI = g["hole_inner_diameter_um"] / 2000.0
    lcMin = cfg["mesh"]["lc_fine_um"] / 1000.0
    lcMax = cfg["mesh"]["lc_coarse_um"] / 1000.0
    Htot = z["drift"]

    gmsh.initialize()
    gmsh.model.add("stack")
    occ = gmsh.model.occ
    occ.addBox(-p / 2, -p / 2, 0, p, p, Htot, tag=1)

    foil_z = [("gem1", z["gem1_bot"]), ("gem2", z["gem2_bot"]), ("gem3", z["gem3_bot"])]
    if not periodic:
        # --- legacy: one aligned hole per foil on the axis --------------------
        for _, zb in foil_z:
            slab = occ.addBox(-p / 2, -p / 2, zb, p, p, tK)
            c1 = occ.addCone(0, 0, zb,          0, 0, tK / 2, rO, rI)
            c2 = occ.addCone(0, 0, zb + tK / 2, 0, 0, tK / 2, rI, rO)
            hole, _ = occ.fuse([(3, c1)], [(3, c2)])
            occ.cut([(3, slab)], hole)
    else:
        # --- periodic: offset hole + its 8 periodic images per foil ----------
        # Cutting the 3x3 lattice of images makes a hole that straddles a wall
        # reappear on the opposite wall, which is also what makes the two walls
        # geometrically identical -> a conformal periodic mesh (see setPeriodic).
        off = offsets_mm or {f: (0.0, 0.0) for f in FOILS}
        for name, zb in foil_z:
            dx, dy = off[name]
            slab = occ.addBox(-p / 2, -p / 2, zb, p, p, tK)
            holes = []
            for ix in (-1, 0, 1):
                for iy in (-1, 0, 1):
                    ox, oy = dx + ix * p, dy + iy * p
                    c1 = occ.addCone(ox, oy, zb,          0, 0, tK / 2, rO, rI)
                    c2 = occ.addCone(ox, oy, zb + tK / 2, 0, 0, tK / 2, rI, rO)
                    h, _ = occ.fuse([(3, c1)], [(3, c2)])
                    holes.extend(h)
            occ.cut([(3, slab)], holes)

    occ.fragment(gmsh.model.occ.getEntities(3), [])
    occ.synchronize()

    gas_vols, kap_vols = [], []
    for _, t in gmsh.model.getEntities(3):
        bb = gmsh.model.getBoundingBox(3, t)
        (kap_vols if (bb[5] - bb[2]) <= tK * 1.5 else gas_vols).append(t)
    gmsh.model.addPhysicalGroup(3, gas_vols, tag=1); gmsh.model.setPhysicalName(3, 1, "gas")
    gmsh.model.addPhysicalGroup(3, kap_vols, tag=2); gmsh.model.setPhysicalName(3, 2, "kapton")

    def faces_at(zval, eps=1e-3):
        return [t for _, t in gmsh.model.getEntitiesInBoundingBox(
            -p, -p, zval - eps, p, p, zval + eps, 2)]

    for name, zval, ptag in [
            ("drift", z["drift"], 11), ("readout", 0.0, 12),
            ("gem1_top", z["gem1_top"], 21), ("gem1_bot", z["gem1_bot"], 22),
            ("gem2_top", z["gem2_top"], 23), ("gem2_bot", z["gem2_bot"], 24),
            ("gem3_top", z["gem3_top"], 25), ("gem3_bot", z["gem3_bot"], 26)]:
        s = faces_at(zval)
        if not s:
            print(f"  WARNING: no surface for electrode '{name}' at z={zval:.4f} mm")
            continue
        gmsh.model.addPhysicalGroup(2, s, tag=ptag); gmsh.model.setPhysicalName(2, ptag, name)

    kap_surfs = set()
    for k in kap_vols:
        for _, t in gmsh.model.getBoundary([(3, k)], oriented=False):
            kap_surfs.add(abs(t))
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "SurfacesList", list(kap_surfs))
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", lcMin)
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", lcMax)
    gmsh.model.mesh.field.setNumber(ft, "DistMin", 2 * rI)
    gmsh.model.mesh.field.setNumber(ft, "DistMax", 6 * rI)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)
    for opt in ("MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"):
        gmsh.option.setNumber("Mesh." + opt, 0)

    # ---- translational-periodic side walls (periodic path only) -------------
    # Tag the four side walls as physical surfaces (so ElmerGrid carries them as
    # named boundaries for the periodic BC) and force the +x/+y walls to be the
    # mesh image of the -x/-y walls under a one-pitch translation, so the
    # discretisation is conformal and Elmer's periodic BC is exact.
    if periodic:
        e = 1e-4

        def wall(axis, val):
            if axis == "x":
                box = (val - e, -p, -e, val + e, p, Htot + e)
            else:
                box = (-p, val - e, -e, p, val + e, Htot + e)
            return [t for _, t in gmsh.model.getEntitiesInBoundingBox(*box, 2)]

        xlo, xhi = wall("x", -p / 2), wall("x", p / 2)
        ylo, yhi = wall("y", -p / 2), wall("y", p / 2)
        for ptag, (nm, s) in zip((31, 32, 33, 34),
                                 [("x_lo", xlo), ("x_hi", xhi),
                                  ("y_lo", ylo), ("y_hi", yhi)]):
            if s:
                gmsh.model.addPhysicalGroup(2, s, tag=ptag)
                gmsh.model.setPhysicalName(2, ptag, nm)
            else:
                print(f"  WARNING: no surface found for periodic wall '{nm}'")
        # affine 4x4 (row-major) translations of +pitch in x and in y
        tX = [1, 0, 0, p,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1]
        tY = [1, 0, 0, 0,  0, 1, 0, p,  0, 0, 1, 0,  0, 0, 0, 1]
        if xlo and xhi:
            gmsh.model.mesh.setPeriodic(2, xhi, xlo, tX)
        if ylo and yhi:
            gmsh.model.mesh.setPeriodic(2, yhi, ylo, tY)

    gmsh.option.setNumber("Mesh.ElementOrder", 2)      # Garfield needs 10-node tets
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # ElmerGrid-friendly
    gmsh.model.mesh.generate(3)
    try:
        gmsh.model.mesh.optimize("HighOrderFast")
    except Exception:
        pass

    os.makedirs(rundir, exist_ok=True)
    msh = os.path.join(rundir, "stack.msh")
    gmsh.write(msh)
    els = gmsh.model.mesh.getElements(3)[1]
    print(f"  mesh: {len(els[0]) if els else 0} volume elements -> {msh}")
    if show_gui:
        gmsh.fltk.run()
    gmsh.finalize()
    return msh


# ----------------------------------------------------------------- Elmer
def parse_mesh_names(path):
    names = {}
    if os.path.exists(path):
        for line in open(path):
            m = re.search(r"\$\s*(\w+)\s*=\s*(\d+)", line)
            if m:
                names[m.group(1)] = int(m.group(2))
    return names


def write_sif(cfg, V, names, rundir, periodic=False):
    eps_k = cfg["geometry"]["kapton_permittivity"]
    gid, kid = names.get("gas", 1), names.get("kapton", 2)
    bcs, n = [], 0
    for nm in ["drift", "readout", "gem1_top", "gem1_bot",
               "gem2_top", "gem2_bot", "gem3_top", "gem3_bot"]:
        if nm not in names:
            print(f"  WARNING: boundary '{nm}' not in mesh.names; skipping its BC")
            continue
        n += 1
        bcs.append(f"Boundary Condition {n}\n"
                   f"  Target Boundaries(1) = {names[nm]}\n"
                   f"  Potential = {V[nm]:.4f}\nEnd\n")

    # --- translational periodic walls (periodic path only) -------------------
    # x_hi is the periodic slave of x_lo (translate maps slave -> master by
    # -pitch in x); likewise y. 'Periodic BC Potential = Logical True' ties the
    # potential across the walls. If the zero-offset validation shows a wrong
    # field, flip the sign of the Translate vectors (master/slave convention).
    if periodic:
        pmm = _pitch_mm(cfg)
        pairs = [("x_lo", "x_hi", (-pmm, 0.0, 0.0)),
                 ("y_lo", "y_hi", (0.0, -pmm, 0.0))]
        for lo, hi, tr in pairs:
            if lo not in names or hi not in names:
                print(f"  WARNING: periodic walls '{lo}'/'{hi}' not in mesh.names;"
                      " periodic BC skipped (field will be wrong for offsets!)")
                continue
            n += 1
            master = n
            bcs.append(f"Boundary Condition {n}\n"
                       f"  Target Boundaries(1) = {names[lo]}\nEnd\n")
            n += 1
            bcs.append(f"Boundary Condition {n}\n"
                       f"  Target Boundaries(1) = {names[hi]}\n"
                       f"  Periodic BC = {master}\n"
                       f"  Periodic BC Potential = Logical True\n"
                       f"  Periodic BC Translate(3) = {tr[0]:.6f} {tr[1]:.6f} {tr[2]:.6f}\n"
                       f"End\n")
    sif = f"""Header
  Mesh DB "." "stack"
End
Simulation
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output File = "case.result"
  Post File = "case.vtu"
End
Constants
  Permittivity Of Vacuum = 8.8542e-12
End
Body 1
  Target Bodies(1) = {gid}
  Equation = 1
  Material = 1
End
Body 2
  Target Bodies(1) = {kid}
  Equation = 1
  Material = 2
End
Material 1
  Relative Permittivity = 1.0
End
Material 2
  Relative Permittivity = {eps_k}
End
Equation 1
  Active Solvers(1) = 1
End
Solver 1
  Equation = "Stat Elec Solver"
  Procedure = "StatElecSolve" "StatElecSolver"
  Variable = Potential
  Variable DOFs = 1
  Calculate Electric Field = True
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU1
  Linear System Max Iterations = 2000
  Linear System Convergence Tolerance = 1.0e-9
  Steady State Convergence Tolerance = 5.0e-8
End

{''.join(bcs)}"""
    with open(os.path.join(rundir, "stack.sif"), "w") as f:
        f.write(sif)


def write_dielectrics(cfg, names, rundir):
    gid, kid = names.get("gas", 1), names.get("kapton", 2)
    eps_k = cfg["geometry"]["kapton_permittivity"]
    with open(os.path.join(rundir, "dielectrics.dat"), "w") as f:
        f.write(f"2\n{gid} 1.0\n{kid} {eps_k}\n")


def solve_field(cfg, V, rundir, skip_grid=False, periodic=False):
    if not skip_grid:
        run_streamed(["ElmerGrid", "14", "2", "stack.msh", "-autoclean", "-names"],
                     cwd=rundir)
    names = parse_mesh_names(os.path.join(rundir, "stack", "mesh.names"))
    if not names:
        print("  WARNING: mesh.names empty; using assumed IDs.")
        names = {"gas": 1, "kapton": 2, "drift": 11, "readout": 12,
                 "gem1_top": 21, "gem1_bot": 22, "gem2_top": 23,
                 "gem2_bot": 24, "gem3_top": 25, "gem3_bot": 26}
        if periodic:
            names.update({"x_lo": 31, "x_hi": 32, "y_lo": 33, "y_hi": 34})
    write_sif(cfg, V, names, rundir, periodic=periodic)
    write_dielectrics(cfg, names, rundir)
    run_streamed(["ElmerSolver", "stack.sif"], cwd=rundir)


def build_field(cfg, rundir, mesh_only=False, no_mesh=False, gui=False,
                offsets_mm=None):
    """Full geometry+field build for one config into its run folder.

    A `misalignment` block in the config switches on the translational-periodic
    path; without it, the legacy aligned mirror cell is built unchanged. For the
    periodic path `offsets_mm` may be passed in (the mc driver does this per
    config); if omitted it is resolved from the config (fixed/aligned modes).
    """
    periodic = has_misalignment(cfg)
    if periodic and offsets_mm is None:
        offsets_mm, _ = resolve_offsets(cfg, 0)

    z, tK = z_landmarks(cfg)
    V = compute_potentials(cfg)
    print("Electrode potentials (V):")
    for k in ["drift", "gem1_top", "gem1_bot", "gem2_top", "gem2_bot",
              "gem3_top", "gem3_bot", "readout"]:
        print(f"  {k:10s} {V[k]:10.2f}")
    if periodic:
        print("Periodicity: TRANSLATIONAL (misalignment-capable)")
        print("Foil offsets: " + offsets_summary(offsets_mm, _pitch_mm(cfg)))
    else:
        print("Periodicity: mirror (aligned cell, legacy)")
    print(f"Stack height: {z['drift']:.3f} mm   ->   {rundir}/")
    os.makedirs(rundir, exist_ok=True)
    if no_mesh:
        print("Reusing existing mesh; re-solving field only.")
    else:
        build_mesh(cfg, z, tK, rundir, offsets_mm=offsets_mm,
                   periodic=periodic, show_gui=gui)
    if mesh_only:
        write_periodicity(rundir, "translate" if periodic else "mirror")
        print("Mesh-only: stopping.")
        return
    solve_field(cfg, V, rundir, skip_grid=no_mesh, periodic=periodic)
    write_periodicity(rundir, "translate" if periodic else "mirror")
    print(f"Field map ready in {rundir}/stack/.")
