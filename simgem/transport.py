"""
Electron transport in the gaps and the misalignment-sensitivity estimate.

Electrons leave a foil at its hole positions and diffuse laterally while
crossing a transfer gap of length L: sigma_T = D_T * sqrt(L). How much of the
upper foil's hole pattern is still imprinted on the arrival positions at the
lower foil is set by the shortest reciprocal-lattice vector G of the patterns:

    memory = exp(-G^2 sigma_T^2 / 2)       (G = 2pi/p square, 4pi/(sqrt3 p) hex)

memory ~ 1: electrons arrive where the holes above put them -> the relative
alignment of the two foils matters. memory << 1: arrivals are uniform over the
lower foil's cell whatever the offset -> alignment cannot change the gain.
The same estimate is given for ions (thermal diffusion), relevant for ion
back-flow.
"""
import csv
import math
import os
import subprocess

from .config import min_reciprocal
from .runner import derive_seed, gas_args
from .util import binary, ncpu

K_B_EV = 8.617333e-5


def classify(mem):
    if mem < 1e-3:
        return "negligible"
    if mem < 0.1:
        return "small"
    return "SIGNIFICANT"


def run(cfg, stack, outdir, electrons=300, distance_um=None):
    """Microscopic drift in every gap's field (one process per distinct
    field); returns (rows by field, report dict)."""
    os.makedirs(outdir, exist_ok=True)
    fields = sorted({round(g["field_kV_cm"] * 1e3, 6) for g in cfg["gaps"]})
    dist_mm = (distance_um / 1000.0) if distance_um else min(
        0.5, min(g["thickness_mm"] for g in cfg["gaps"]))
    exe = binary("simgem_transport")
    by_field = {}
    pending = list(enumerate(fields))
    while pending:                      # at most ncpu() processes at a time
        batch, pending = pending[:ncpu()], pending[ncpu():]
        procs = []
        for i, E in batch:
            out = os.path.join(outdir, f".transport_{i}.csv")
            cmd = [exe] + gas_args(cfg, penning=False) + [
                "--fields", repr(E), "--distance", f"{dist_mm / 10:.6f}",
                "--electrons", str(electrons), "--seed", str(derive_seed("transport", E)),
                "--out", out]
            procs.append((E, out, subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                                   stderr=subprocess.STDOUT, text=True)))
        for E, out, p in procs:
            txt, _ = p.communicate()
            if p.returncode != 0:
                raise RuntimeError(f"simgem_transport failed:\n{txt[-2000:]}")
            with open(out) as f:
                by_field[E] = next(csv.DictReader(f))     # one field per process
            os.remove(out)
    rows = [by_field[E] for E in fields]
    path = os.path.join(outdir, "transport.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for E in fields:
            w.writerow(by_field[E])
    return by_field, sensitivity(cfg, stack, by_field)


def sensitivity(cfg, stack, by_field):
    T = cfg["gas"]["temperature_K"]
    gaps = []
    for i, g in enumerate(cfg["gaps"]):
        E = g["field_kV_cm"] * 1e3
        r = by_field[round(E, 6)]
        L_cm = g["thickness_mm"] / 10.0
        DT = float(r["D_T_sqrt_cm"])
        sig_e = DT * math.sqrt(L_cm) * 10.0                         # mm
        sig_i = math.sqrt(2 * K_B_EV * T * L_cm / E) * 10.0          # mm
        entry = {"gap": g["name"], "thickness_mm": g["thickness_mm"], "field_V_cm": E,
                 "v_drift_um_ns": float(r["v_drift_cm_ns"]) * 1e4,
                 "D_T_um_sqrt_cm": DT * 1e4, "sigma_T_e_um": sig_e * 1e3,
                 "sigma_T_ion_um": sig_i * 1e3,
                 "arrived_fraction": int(r["arrived"]) / max(int(r["electrons"]), 1)}
        if 0 < i < len(cfg["gaps"]) - 1:            # transfer gap: foil above & below
            above, below = stack.foils[i - 1], stack.foils[i]
            G = max(min_reciprocal(above), min_reciprocal(below))  # 1/mm
            entry["between"] = [above["name"], below["name"]]
            entry["sigma_over_pitch"] = sig_e / (below["pitch_um"] / 1000.0)
            entry["memory_electrons"] = math.exp(-0.5 * (G * sig_e) ** 2)
            entry["memory_ions"] = math.exp(-0.5 * (G * sig_i) ** 2)
        gaps.append(entry)
    return {"gaps": gaps}


def format_report(rep):
    lines = ["Transport in the gaps (microscopic, uniform field):",
             f"  {'gap':<11}{'E [V/cm]':>10}{'v_d [um/ns]':>13}{'D_T [um/rcm]':>14}"
             f"{'sigma_T [um]':>14}{'arrived':>9}"]
    for g in rep["gaps"]:
        lines.append(f"  {g['gap']:<11}{g['field_V_cm']:>10.0f}{g['v_drift_um_ns']:>13.3f}"
                     f"{g['D_T_um_sqrt_cm']:>14.0f}{g['sigma_T_e_um']:>14.1f}"
                     f"{100 * g['arrived_fraction']:>8.1f}%")
    tg = [g for g in rep["gaps"] if "memory_electrons" in g]
    if tg:
        lines.append("Hole-pattern memory across transfer gaps (1 = fully imprinted, 0 = washed out):")
        for g in tg:
            me, mi = g["memory_electrons"], g["memory_ions"]
            lines.append(f"  {g['gap']:<11}{g['between'][0]}->{g['between'][1]}: "
                         f"sigma_T/pitch = {g['sigma_over_pitch']:.2f}, "
                         f"electrons {me:.1e} ({classify(me)}), "
                         f"ions {mi:.1e} ({classify(mi)})")
        worst = max(g["memory_electrons"] for g in tg)
        if worst < 1e-3:
            lines.append("  => foil misalignment cannot measurably change the electron gain here "
                         "(it can still matter for ion back-flow if ion memory is large).")
        elif worst < 0.1:
            lines.append("  => misalignment effects on the gain are expected to be small.")
        else:
            lines.append("  => misalignment CAN change the gain in this configuration.")
    return "\n".join(lines)
