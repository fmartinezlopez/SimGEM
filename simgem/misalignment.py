"""The foil-misalignment Monte Carlo (parameter scans live in sweep.py)."""
import csv
import json
import math
import os
import random

from . import config as C
from . import field, runner, transport
from .geometry import wrap_offset
from .config import lattice
from .analysis import chi2_sf


def _cached_summary(outdir, key, events, label="gain"):
    """A previous result in `outdir` with the same run key and length."""
    p = os.path.join(outdir, f"{label}_summary.json")
    if os.path.exists(p):
        with open(p) as f:
            s = json.load(f)
        if s.get("run_key") == key and s.get("events") == events:
            return s
    return None


# ----------------------------------------------------------------- misalignment MC
def draw_offsets(cfg, mc, index):
    """Offsets (mm) of geometry `index`: nominal + a random perturbation for
    the selected foils. Deterministic in (mc.seed, index)."""
    rng = random.Random(runner.derive_seed("mc", mc["seed"], index))
    chosen = set(mc["foils"] or [f["name"] for f in cfg["foils"]])
    out = {}
    for f in cfg["foils"]:
        nx, ny = (v / 1000.0 for v in f["offset_um"])
        if f["name"] in chosen:
            if mc["distribution"] == "gaussian":
                s = mc["sigma_um"] / 1000.0
                dx, dy = rng.gauss(0, s), rng.gauss(0, s)
            else:
                _, _, cx, cy = lattice(f)
                dx, dy = rng.uniform(-cx / 2, cx / 2), rng.uniform(-cy / 2, cy / 2)
            nx, ny = nx + dx, ny + dy
        out[f["name"]] = wrap_offset(f, (nx, ny))
    return out


def mc_statistics(geoms, ref):
    """Separate geometry-induced spread from avalanche statistics.
    geoms: list of per-geometry summaries; ref: reference summary or None."""
    m = [g["effective_gain"] for g in geoms]
    e = [g["effective_gain_err"] for g in geoms]
    K = len(m)
    mean = sum(m) / K
    s2 = sum((x - mean) ** 2 for x in m) / (K - 1)          # observed variance of means
    w2 = sum(x * x for x in e) / K                          # expected from statistics alone
    tau2 = max(s2 - w2, 0.0)
    # rough 1-sigma uncertainty on the variance estimate (not below the
    # statistical level, or a lucky small spread gives a too-tight limit)
    tau2_err = math.sqrt(2.0 / (K - 1)) * max(s2, w2)
    chi2 = sum((x - mean) ** 2 / (err * err) for x, err in zip(m, e) if err > 0)
    # error on the mean: the observed spread, but never less than the
    # statistical one (with few geometries s2 can come out accidentally small)
    mean_err = math.sqrt(max(s2, w2) / K)
    out = {
        "geometries": K,
        "mean_of_geometries": mean, "mean_of_geometries_err": mean_err,
        "spread_observed_rms": math.sqrt(s2),
        "spread_statistical_rms": math.sqrt(w2),
        "spread_geometry_rms": math.sqrt(tau2),
        "spread_geometry_rms_upper95": math.sqrt(tau2 + 1.645 * tau2_err),
        "chi2": chi2, "ndof": K - 1, "p_value_no_geometry_effect": chi2_sf(chi2, K - 1),
    }
    if ref and ref.get("events"):
        d = mean - ref["effective_gain"]
        de = math.sqrt(mean_err ** 2 + ref["effective_gain_err"] ** 2)
        out.update({"reference": ref["effective_gain"], "reference_err": ref["effective_gain_err"],
                    "shift": d, "shift_err": de, "shift_sigma": d / de if de > 0 else 0.0,
                    "shift_relative": d / ref["effective_gain"] if ref["effective_gain"] else 0.0})
    return out


def format_mc(st, mc):
    dist = (f"gaussian sigma = {mc['sigma_um']} um" if mc["distribution"] == "gaussian"
            else "uniform over the cell")
    L = ["=" * 68, f"MISALIGNMENT MC  ({st['geometries']} geometries, {dist}, "
         f"foils: {', '.join(mc['foils']) if mc['foils'] else 'all'})", "-" * 68,
         f"mean over geometries : {st['mean_of_geometries']:.4g} +/- {st['mean_of_geometries_err']:.2g}"]
    if "reference" in st:
        L.append(f"nominal geometry     : {st['reference']:.4g} +/- {st['reference_err']:.2g}")
        L.append(f"shift                : {st['shift']:+.3g} +/- {st['shift_err']:.2g} "
                 f"({st['shift_sigma']:+.1f} sigma, {100 * st['shift_relative']:+.2f}%)")
    L += ["-" * 68,
          f"spread of geometry means: observed RMS {st['spread_observed_rms']:.3g}, "
          f"expected from avalanche statistics {st['spread_statistical_rms']:.3g}",
          f"  -> geometry-induced RMS {st['spread_geometry_rms']:.3g} "
          f"(95% upper limit {st['spread_geometry_rms_upper95']:.3g})",
          f"  chi2/ndof = {st['chi2']:.1f}/{st['ndof']}, "
          f"p(no geometry effect) = {st['p_value_no_geometry_effect']:.3g}",
          "=" * 68]
    return "\n".join(L)


def run_mc(cfg, rundir, geometries=None, events=None, reference_events=None,
                    jobs=None, skip_transport=False, force=False):
    mc = cfg.get("misalignment_mc")
    if not mc:
        raise C.ConfigError("this config has no 'misalignment_mc' block")
    K = geometries or mc["geometries"]
    M = events or mc["events_per_geometry"]
    Nref = reference_events or mc["reference_events"] or K * M
    mdir = os.path.join(rundir, "mc")
    os.makedirs(mdir, exist_ok=True)
    nominal = C.Stack(cfg)

    if not skip_transport:
        print("Checking whether misalignment can matter at all (transport in the gaps)...")
        _, rep = transport.run(cfg, nominal, mdir)
        print(transport.format_report(rep))
        with open(os.path.join(mdir, "transport.json"), "w") as f:
            json.dump(rep, f, indent=2)

    geoms = []
    names = [f["name"] for f in cfg["foils"]]
    for i in range(K):
        off = draw_offsets(cfg, mc, i)
        st = C.Stack(cfg, offsets_mm=off)
        desc = "  ".join(f"{n}=({off[n][0] * 1e3:+.1f},{off[n][1] * 1e3:+.1f})um" for n in names)
        print(f"\n--- geometry {i + 1}/{K}: {desc}")
        h = field.build(st)
        gdir = os.path.join(mdir, f"g{i:04d}")
        seed = runner.derive_seed(cfg["simulation"]["seed"], "mc", i)
        s = None if force else _cached_summary(gdir, runner.run_key(cfg, h, seed), M)
        if s:
            print("  (cached result)")
        else:
            s = runner.run_gain(cfg, st, h, gdir, M, jobs=jobs, seed=seed)
        # offsets as meshed (a hole/wall tangency may have nudged them by ~um)
        used = {k: tuple(v) for k, v in field.load_meta(h)["offsets_mm"].items()}
        s = dict(s, index=i, offsets_mm=used)
        geoms.append(s)
        print(f"  effective gain {s['effective_gain']:.4g} +/- {s['effective_gain_err']:.2g}")

    print(f"\n--- nominal geometry: {Nref} events ---")
    h0 = field.build(nominal)
    rdir = os.path.join(mdir, "nominal")
    key0 = runner.run_key(cfg, h0, cfg["simulation"]["seed"])
    ref = None if force else _cached_summary(rdir, key0, Nref)
    ref = ref or runner.run_gain(cfg, nominal, h0, rdir, Nref, jobs=jobs)

    with open(os.path.join(mdir, "mc_geometries.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["geometry"] + [f"{n}_{a}_um" for n in names for a in ("dx", "dy")] +
                   ["field_hash", "events", "effective_gain", "effective_gain_err",
                    "absolute_gain", "absolute_gain_err"])
        for g in geoms:
            w.writerow([g["index"]] + [f"{v * 1e3:.3f}" for n in names for v in g["offsets_mm"][n]] +
                       [g["field_hash"], g["events"], g["effective_gain"], g["effective_gain_err"],
                        g["absolute_gain"], g["absolute_gain_err"]])
    stats = mc_statistics(geoms, ref)
    text = format_mc(stats, mc)
    print("\n" + text)
    with open(os.path.join(mdir, "mc_summary.json"), "w") as f:
        json.dump({"settings": mc, "statistics": stats}, f, indent=2)
    with open(os.path.join(mdir, "mc_summary.txt"), "w") as f:
        f.write(text + "\n")
    print(f"Wrote {mdir}/mc_geometries.csv, mc_summary.json, mc_summary.txt")
    return stats
