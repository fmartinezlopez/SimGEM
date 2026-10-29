"""
Field maps: Elmer input, solve, validation, and the content-addressed cache.

A field map depends only on the geometry, the electrode potentials and the
mesh settings, so it is stored under fields/<hash>/ (hash of exactly those
inputs) and shared by every config, scan point and misalignment geometry
that needs it. A map is only ever used after it has passed `check_field`.
"""
import datetime
import json
import os
import re
import shutil
import time

from . import geometry
from .util import run_streamed

# Relative to the working directory (the install may be read-only, e.g. a SIF).
FIELDS_DIR = os.path.abspath(os.environ.get("SIMGEM_FIELDS", "fields"))
PERIODIC_TOL_V = 1e-2      # max allowed potential mismatch across periodic walls


def field_dir(h):
    return os.path.join(FIELDS_DIR, h)


def load_meta(h):
    p = os.path.join(field_dir(h), "field.json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def is_ready(h):
    meta = load_meta(h)
    return bool(meta and meta.get("check", {}).get("ok"))


# ----------------------------------------------------------------- Elmer input
def parse_mesh_names(path):
    names = {}
    with open(path) as f:
        for line in f:
            m = re.search(r"\$\s*(\w+)\s*=\s*(\d+)", line)
            if m:
                names[m.group(1)] = int(m.group(2))
    return names


def write_elmer_input(stack, names, d, vtu=False):
    bodies = [("gas", 1.0)] + [(f"kapton_{f['name']}", f["kapton_permittivity"])
                               for f in stack.foils]
    with open(os.path.join(d, "dielectrics.dat"), "w") as fh:
        fh.write(f"{len(bodies)}\n")
        for n, eps in bodies:
            fh.write(f"{names[n]} {eps}\n")

    blocks = []
    for i, (n, eps) in enumerate(bodies, 1):
        blocks.append(f"Body {i}\n  Target Bodies(1) = {names[n]}\n"
                      f"  Equation = 1\n  Material = {i}\nEnd\n"
                      f"Material {i}\n  Relative Permittivity = {eps}\nEnd\n")
    bc = 0
    for n in stack.electrodes():
        bc += 1
        blocks.append(f"Boundary Condition {bc}\n  Target Boundaries(1) = {names[n]}\n"
                      f"  Potential = {stack.V[n]:.6f}\nEnd\n")
    # Translational periodicity. No 'Periodic BC Translate': Elmer finds the
    # mapping itself. (An explicit vector with the wrong sign makes Elmer
    # match 0 nodes and silently produce a garbage field.)
    for lo, hi in (("x_lo", "x_hi"), ("y_lo", "y_hi")):
        bc += 1
        master = bc
        blocks.append(f"Boundary Condition {bc}\n  Target Boundaries(1) = {names[lo]}\nEnd\n")
        bc += 1
        blocks.append(f"Boundary Condition {bc}\n  Target Boundaries(1) = {names[hi]}\n"
                      f"  Periodic BC = {master}\n  Periodic BC Potential = Logical True\nEnd\n")

    post = '  Post File = "case.vtu"\n' if vtu else ""
    sif = f"""Header
  Mesh DB "." "stack"
End
Simulation
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output File = "case.result"
{post}End
Constants
  Permittivity Of Vacuum = 8.8542e-12
End
Equation 1
  Active Solvers(1) = 1
End
Solver 1
  Equation = "Stat Elec Solver"
  Procedure = "StatElecSolve" "StatElecSolver"
  Variable = Potential
  Variable DOFs = 1
  Calculate Electric Field = {"True" if vtu else "False"}
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU1
  Linear System Max Iterations = 5000
  Linear System Convergence Tolerance = 1.0e-10
  Linear System Abort Not Converged = True
  Steady State Convergence Tolerance = 5.0e-8
End

{"".join(blocks)}"""
    with open(os.path.join(d, "stack.sif"), "w") as fh:
        fh.write(sif)


# ----------------------------------------------------------------- validation
def read_nodes(d):
    nodes = {}
    with open(os.path.join(d, "stack", "mesh.nodes")) as f:
        for line in f:
            p = line.split()
            nodes[int(p[0])] = (float(p[2]), float(p[3]), float(p[4]))
    return nodes


def read_potential(d):
    """Node id -> potential from Elmer's case.result (first 'Perm:' block,
    which is what Garfield reads too). Values follow the node order of the
    permutation table."""
    with open(os.path.join(d, "stack", "case.result")) as f:
        for line in f:
            if line.startswith("Perm:"):
                n = int(line.split()[2])
                break
        else:
            raise RuntimeError("case.result has no 'Perm:' block")
        ids = [int(next(f).split()[0]) for _ in range(n)]
        vals = [float(next(f)) for _ in range(n)]
    return dict(zip(ids, vals))


def check_field(stack, d, solver_log):
    """Validate a solved field map. Returns a dict with 'ok' and diagnostics."""
    res = {"ok": True, "problems": []}
    with open(solver_log) as f:
        txt = f.read()
    m = re.search(r"Points not found:\s*(\d+)", txt)
    if m and int(m.group(1)) > 0:
        res["ok"] = False
        res["problems"].append(f"Elmer periodic projector could not match {m.group(1)} wall nodes")
    if "ALL DONE" not in txt:
        res["ok"] = False
        res["problems"].append("ElmerSolver did not finish")
        return res

    nodes = read_nodes(d)
    V = read_potential(d)
    Lx, Ly = stack.cell
    tol = 1e-7
    for axis, L, nm in ((0, Lx, "x"), (1, Ly, "y")):
        o = [i for i in range(3) if i != axis]
        lo, hi = {}, []
        for nid, c in nodes.items():
            if abs(c[axis] + L / 2) < tol:
                lo[(round(c[o[0]], 6), round(c[o[1]], 6))] = V[nid]
            elif abs(c[axis] - L / 2) < tol:
                hi.append(((round(c[o[0]], 6), round(c[o[1]], 6)), V[nid]))
        diffs = [abs(v - lo[k]) for k, v in hi if k in lo]
        frac = len(diffs) / max(len(hi), 1)
        dmax = max(diffs) if diffs else float("inf")
        res[f"wall_{nm}"] = {"nodes": len(hi), "matched_fraction": round(frac, 5),
                             "max_dV": dmax}
        if frac < 0.99:
            res["ok"] = False
            res["problems"].append(f"{nm}-walls: only {100 * frac:.1f}% of nodes pair up (mesh not periodic)")
        if dmax > PERIODIC_TOL_V:
            res["ok"] = False
            res["problems"].append(f"{nm}-walls: potential differs by up to {dmax:.3g} V "
                                   "across the periodic boundary")
    vmin, vmax = min(V.values()), max(V.values())
    res["V_range"] = [vmin, vmax]
    lo_ok = min(stack.V.values()) - 1e-3 * abs(stack.V["drift"])
    if vmin < lo_ok or vmax > 1e-3 * abs(stack.V["drift"]):
        res["ok"] = False
        res["problems"].append(f"potential range [{vmin:.1f}, {vmax:.1f}] V exceeds the electrode potentials")
    return res


# ----------------------------------------------------------------- build
def _filter_elmergrid(line):
    return not line.startswith("Untreated command")


def build(stack, h=None, force=False, vtu=False, keep_msh=False, log=print):
    """Ensure the field map for `stack` exists in the cache; return its hash."""
    h = h or stack.field_hash()
    final = field_dir(h)
    if is_ready(h) and not force:
        log(f"Field map {h} is cached ({final}).")
        return h
    if os.path.exists(final):
        shutil.rmtree(final)
    os.makedirs(FIELDS_DIR, exist_ok=True)
    tmp = os.path.join(FIELDS_DIR, f".build-{h}-{os.getpid()}")
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)

    t0 = time.time()
    log(f"Building field map {h} in {tmp}")
    offsets = geometry.safe_offsets(stack, log=log)
    info = geometry.build_mesh(stack, offsets, tmp, log=log)
    t_mesh = time.time() - t0
    run_streamed(["ElmerGrid", "14", "2", "stack.msh", "-autoclean", "-names"],
                 cwd=tmp, log_path=os.path.join(tmp, "elmergrid.log"), keep=_filter_elmergrid,
                 echo=False)
    names = parse_mesh_names(os.path.join(tmp, "stack", "mesh.names"))
    missing = [n for n in info["physical"] if n not in names]
    if missing:
        raise RuntimeError(f"ElmerGrid lost physical groups {missing}")
    write_elmer_input(stack, names, tmp, vtu=vtu)
    solver_log = os.path.join(tmp, "elmersolver.log")
    run_streamed(["ElmerSolver", "stack.sif"], cwd=tmp, log_path=solver_log, echo=False)
    t_solve = time.time() - t0 - t_mesh

    log("Checking the solution (periodic walls, potential range) ...")
    chk = check_field(stack, tmp, solver_log)
    for k in ("wall_x", "wall_y"):
        if k in chk:
            w = chk[k]
            log(f"  {k}: {w['nodes']} nodes, {100 * w['matched_fraction']:.2f}% paired, "
                f"max |dV| = {w['max_dV']:.2e} V")
    if not keep_msh:
        os.remove(os.path.join(tmp, "stack.msh"))

    meta = {
        "hash": h,
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "key": stack.field_key(),
        "cell_mm": list(stack.cell),
        "z_mm": stack.z,
        "kapton_mm": stack.kapton,
        "potentials_V": stack.V,
        "offsets_mm": {k: list(v) for k, v in offsets.items()},
        "holes_mm": {f["name"]: geometry.hole_centres(f, offsets[f["name"]], stack.cell)
                     for f in stack.foils},
        "mesh": {"elements": info["elements"], "nodes": info["nodes"]},
        "timing_s": {"mesh": round(t_mesh, 1), "solve": round(t_solve, 1)},
        "check": chk,
    }
    with open(os.path.join(tmp, "field.json"), "w") as f:
        json.dump(meta, f, indent=2)
    if not chk["ok"]:
        bad = os.path.join(FIELDS_DIR, f".failed-{h}")
        if os.path.exists(bad):
            shutil.rmtree(bad)
        os.rename(tmp, bad)
        raise RuntimeError("field map FAILED validation:\n  " + "\n  ".join(chk["problems"])
                           + f"\n(kept for inspection in {bad})")
    try:
        os.rename(tmp, final)
    except OSError:      # someone else built it meanwhile
        shutil.rmtree(tmp)
    log(f"Field map {h} ready ({t_mesh:.0f} s mesh, {t_solve:.0f} s solve).")
    return h


def preview_mesh(stack, outdir, gui=False, log=print):
    """Mesh only (no solve, not cached) for a look in the gmsh GUI."""
    offsets = geometry.safe_offsets(stack, log=log)
    geometry.build_mesh(stack, offsets, outdir, show_gui=gui, log=log)


def list_cache():
    out = []
    if not os.path.isdir(FIELDS_DIR):
        return out
    for h in sorted(os.listdir(FIELDS_DIR)):
        meta = load_meta(h)
        if not meta:
            continue
        size = sum(os.path.getsize(os.path.join(dp, fn))
                   for dp, _, fns in os.walk(field_dir(h)) for fn in fns)
        out.append((h, meta, size))
    return out
