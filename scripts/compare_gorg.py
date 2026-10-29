#!/usr/bin/env python3
"""
Compare a SimGEM gain-curve scan with the GORG measurements.

    python scripts/compare_gorg.py runs/scans/gorg/gain_curve_14p5psi \\
        [--gorg ../GORG] [--out runs/scans/gorg/gain_curve_14p5psi/compare]

Reads the measured gains from GORG/db/gain_results.csv (gas, gas ratio and
nominal pressure taken from the simulated configs), keeps runs with quality
'ok' and no fit flags, and writes:

  <out>.png  gain vs V_GEM (data and simulation, log scale) + simulation/data ratio
  <out>.csv  per simulated point: V_GEM, simulated gain, measured gain(s), ratio

and prints exponential fits G = exp(a + b V) of both curves: b is the
logarithmic slope (the gain factor per 100 V is exp(100 b)).

The simulated gain is the zero-inflated-Polya gain (mean over all events),
the quantity GORG's photopeak-based gain measures. Needs numpy, pandas,
matplotlib (e.g. the base conda env).
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from simgem import analysis, corrections  # noqa: E402  (pure python)
from simgem.config import BAR_PER_PSI as PSI, GAS_DISPLAY  # noqa: E402



def load_sim(scan_dir, include_partial=False, corrected=False):
    """Finished points; with include_partial also points that were interrupted
    (events read from the workers' part files; marked partial=True). With
    corrected, gains include the config's empirical corrections
    (simgem/corrections.py; factor 1 for configs without them)."""
    rows = []
    for pdir in sorted(d for d in glob.glob(os.path.join(scan_dir, "p[0-9][0-9][0-9]_*")) if os.path.isdir(d)):
        try:
            cfg = json.load(open(os.path.join(pdir, "config.json")))["resolved"]
        except FileNotFoundError:
            continue
        partial = not os.path.exists(os.path.join(pdir, "gain_summary.json"))
        if partial:
            if not include_partial:
                continue
            vals = []
            for part in sorted(glob.glob(os.path.join(pdir, ".gain_w*_*.csv"))):
                with open(part) as f:
                    lines = f.read().splitlines()[1:]
                vals += [int(l.split(",")[4]) for l in lines if l.count(",") == 6]
            if len(vals) < 20:
                continue
            fit = analysis.fit_gain(vals)
            if not fit.get("ok"):
                fit = {"ok": True, "gain": fit["mean_all"], "gain_err": fit["mean_all_err"]}
            summ = {"events": len(vals)}
        else:
            summ = json.load(open(os.path.join(pdir, "gain_summary.json")))
            fit = {}
            if os.path.exists(os.path.join(pdir, "gain_fit.json")):
                fit = json.load(open(os.path.join(pdir, "gain_fit.json")))
        ok = fit.get("ok") and "gain" in fit
        k = corrections.factor(cfg) if corrected else 1.0
        comps = cfg["gas"]["components"]
        tot = sum(comps.values())
        rows.append({
            "V_gem": cfg["hv"]["V_gem"],
            "gain": k * (fit["gain"] if ok else summ["effective_gain"]),
            "gain_err": k * (fit["gain_err"] if ok else summ["effective_gain_err"]),
            "gain_correction": k,
            "events": summ["events"],
            "partial": partial,
            "pressure_psi": round(cfg["gas"]["pressure_bar"] / PSI, 1),
            "gas": "/".join(GAS_DISPLAY.get(k, k) for k in comps),
            "gas_ratio": "/".join(f"{100 * v / tot:g}" for v in comps.values()),
            "penning_r": cfg["gas"]["penning"]["r"],
        })
    if not rows:
        sys.exit(f"no finished points in {scan_dir}")
    return pd.DataFrame(rows).sort_values("V_gem")


def load_data(gorg, gas, ratio, psi, adc_offset=None):
    """GORG gains (quality ok, no fit flags). With adc_offset (channels), the zero
    offset of the pulse-height scale is removed from both peaks before the
    source/pulser ratio, exactly as gorgana.analysis.systematics.apply_adc_offset
    (August 2025 reference; GORG measures +18.5 +- 0.8 ch from the escape peak)."""
    g = pd.read_csv(os.path.join(gorg, "db", "gain_results.csv"))
    # nominal pressure (folder); the measured one only if no run has that nominal value
    g = g[(g.gas == gas) & (g.gas_ratio == ratio)]
    sel = g[np.isclose(g.p_nominal_psi, psi)]
    if sel.empty:
        sel = g[np.isclose(g.pressure_psi, psi)]
    if "quality" in sel:
        sel = sel[sel.quality.fillna("ok") == "ok"]
    if "fit_flags" in sel:
        sel = sel[sel.fit_flags.isna()]
    sel = sel.copy()
    sel["gain_uncorrected"] = sel["gain"]
    if adc_offset:
        offsets = {"aug2025": float(adc_offset)}
        try:
            from gorgana.analysis.systematics import apply_adc_offset
            sel = apply_adc_offset(sel, offsets)
        except ImportError:                      # same formula, without gorgana
            off = sel["reference"].map(offsets).fillna(0.0)
            scale = ((sel.source_mu - off) / (sel.pulser_mu - off)) / (sel.source_mu / sel.pulser_mu)
            sel["gain"] = sel["gain"] * scale.where(np.isfinite(scale), 1.0)
    return sel[["session_date", "V_gem", "gain", "gain_uncorrected"]].sort_values("V_gem")


def expfit(V, G, err=None):
    """Weighted least squares of ln G = a + b V; returns a, b, b_err."""
    V, lnG = np.asarray(V, float), np.log(np.asarray(G, float))
    w = np.ones_like(V) if err is None else (np.asarray(G) / np.asarray(err)) ** 2
    A = np.vstack([np.ones_like(V), V]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    a, b = cov @ (A.T @ (w * lnG))
    if err is None:                         # scale by the scatter
        res = lnG - (a + b * V)
        cov = cov * (res @ res) / max(len(V) - 2, 1)
    return a, b, math.sqrt(cov[1, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan_dirs", nargs="+", help="one or more gain-curve scan folders (e.g. one per pressure)")
    ap.add_argument("--gorg", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                   "..", "..", "GORG"))
    ap.add_argument("--adc-offset", type=float, metavar="CH",
                    help="correct the GORG gains for a zero offset of the pulse-height scale "
                         "(channels; GORG's escape-peak value for August 2025 is 18.5)")
    ap.add_argument("--corrected", action="store_true",
                    help="show simulated gains with the configs' empirical corrections "
                         "(charging-up, per-foil scale; see simgem/corrections.py)")
    ap.add_argument("--data-at-sim", action="store_true",
                    help="use only the measured runs at the simulated V_GEM values")
    ap.add_argument("--include-partial", action="store_true",
                    help="also show unfinished points (open markers; not used in the fits)")
    ap.add_argument("--out", help="output stem (default <scan_dir>/compare_gorg, or "
                                  "<parent>/compare_gorg_<pressures> for several scans)")
    args = ap.parse_args()

    curves = []
    for d in args.scan_dirs:
        sim_all = load_sim(d, args.include_partial, args.corrected)
        sim = sim_all[~sim_all.partial]
        gas, ratio, psi = sim.gas.iloc[0], sim.gas_ratio.iloc[0], sim.pressure_psi.iloc[0]
        data = load_data(args.gorg, gas, ratio, psi, args.adc_offset)
        if args.data_at_sim:
            data = data[data.V_gem.isin(sim.V_gem)]
        if data.empty:
            sys.exit(f"no GORG gains for {gas} {ratio} at {psi} psi")
        meas = data.groupby("V_gem")["gain"].agg(["mean", "std", "count"])
        meas["correction"] = (data.groupby("V_gem")["gain"].mean()
                              / data.groupby("V_gem")["gain_uncorrected"].mean())
        tab = sim_all.merge(meas, left_on="V_gem", right_index=True, how="left").rename(
            columns={"mean": "gain_data", "std": "gain_data_spread", "count": "n_data_runs"})
        tab["ratio_sim_data"] = tab["gain"] / tab["gain_data"]
        tab["ratio_err"] = tab["gain_err"] / tab["gain_data"]
        fs, fd = expfit(sim.V_gem, sim.gain, sim.gain_err), expfit(data.V_gem, data.gain)
        curves.append(dict(dir=d, sim=sim, sim_all=sim_all, data=data, tab=tab, gas=gas, ratio=ratio, psi=psi, fs=fs, fd=fd))

    if args.out:
        out = args.out
    elif len(curves) == 1:
        out = os.path.join(args.scan_dirs[0], "compare_gorg")
    else:
        tag = "_".join(f"{c['psi']:g}".replace(".", "p") for c in curves)
        out = os.path.join(os.path.dirname(os.path.abspath(args.scan_dirs[0])), f"compare_gorg_{tag}psi")
    if args.adc_offset and not args.out:
        out += f"_adcoffset{args.adc_offset:g}".replace(".", "p")
    if args.corrected and not args.out:
        out += "_corrected"
    pd.concat([c["tab"] for c in curves]).to_csv(out + ".csv", index=False)

    for c in curves:
        sim, tab = c["sim"], c["tab"]
        (a_s, b_s, eb_s), (a_d, b_d, eb_d) = c["fs"], c["fd"]
        print(f"\n{c['gas']} {c['ratio']} at {c['psi']:g} psi, Penning r = {sim.penning_r.iloc[0]}  "
              + (f"(simulation x{sim.gain_correction.iloc[0]:.3g}, empirical corrections)  "
                 if args.corrected else "")
              + f"({len(sim)} simulated points, {len(c['data'])} measured runs)")
        corr_col = args.adc_offset is not None
        print(f"{'V_gem':>7} {'sim gain':>10} {'+/-':>8} {'data gain':>10} {'sim/data':>9} {'+/-':>7}"
              + (f" {'data corr.':>10}" if corr_col else ""))
        for _, r in tab.iterrows():
            print(f"{r.V_gem:>7.0f} {r.gain:>10.4g} {r.gain_err:>8.2g} {r.gain_data:>10.4g} "
                  f"{r.ratio_sim_data:>9.3f} {r.ratio_err:>7.3f}"
                  + (f" {100 * (r.correction - 1):>+9.1f}%" if corr_col and pd.notna(r.correction) else "")
                  + (f"   PARTIAL ({r.events} events; not in fits)" if r.partial else ""))
        m = tab.ratio_sim_data.notna() & ~tab.partial
        w = 1 / tab.ratio_err[m] ** 2
        mean_ratio = (tab.ratio_sim_data[m] * w).sum() / w.sum()
        print(f"mean sim/data {mean_ratio:.3f} +/- {w.sum() ** -0.5:.3f}")
        print(f"slope (ln G per kV): simulation {1e3 * b_s:.2f} +/- {1e3 * eb_s:.2f}, "
              f"data {1e3 * b_d:.2f} +/- {1e3 * eb_d:.2f}   "
              f"(gain x{math.exp(100 * b_s):.2f} vs x{math.exp(100 * b_d):.2f} per 100 V)")
        mid = 0.5 * (sim.V_gem.min() + sim.V_gem.max())
        dV = ((a_d + b_d * mid) - (a_s + b_s * mid)) / b_s
        print(f"at {mid:.0f} V the simulation needs {dV:+.0f} V more V_GEM to reach the measured gain")

    plot(curves, out, args.adc_offset, args.corrected)


def plot(curves, out, adc_offset=None, corrected=False):
    """GORGAna style (gorgana.style: fonts, ticks, pressure colour/marker registry,
    experiment label, PDF + PNG) when GORG's package is importable, else plain matplotlib.
    Marker fill = source (open: GORG, filled: SimGEM, half: unfinished SimGEM point);
    colour + marker shape = pressure."""
    import matplotlib
    matplotlib.use("Agg")
    try:
        from gorgana import style as gs
    except ImportError:
        gs = None
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    def pstyle(k, psi):
        if gs:
            return gs.pressure_style(psi)
        return f"C{k}", "osD^v"[k % 5]

    tex = gs.tex if gs else (lambda t: t)
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.25], "hspace": 0.05})
    ys = []
    for k, c in enumerate(curves):
        col, mk = pstyle(k, c["psi"])
        sim, data, tab, part = c["sim"], c["data"], c["tab"], c["sim_all"][c["sim_all"].partial]
        (a_s, b_s, _), (a_d, b_d, _) = c["fs"], c["fd"]
        ax.plot(data.V_gem, data.gain, mk, ls="none", color=col, mfc="white", mew=1.4, zorder=3)
        ax.errorbar(sim.V_gem, sim.gain, yerr=sim.gain_err, fmt=mk, color=col, mfc=col, mew=1.4, zorder=4)
        if len(part):
            ax.errorbar(part.V_gem, part.gain, yerr=part.gain_err, fmt=mk, color=col, mew=1.4,
                        fillstyle="left", mfc=col, mfcalt="white", zorder=4)
        vd = np.linspace(data.V_gem.min(), data.V_gem.max(), 50)
        vs = np.linspace(sim.V_gem.min(), max(sim.V_gem.max(), part.V_gem.max() if len(part) else 0), 50)
        ax.plot(vd, np.exp(a_d + b_d * vd), "-", color=col, lw=1.2, zorder=2)
        ax.plot(vs, np.exp(a_s + b_s * vs), "--", color=col, lw=1.2, zorder=2)
        ys += list(data.gain) + list(sim.gain)
        m = tab.ratio_sim_data.notna() & ~tab.partial
        axr.errorbar(tab.V_gem[m], tab.ratio_sim_data[m], yerr=tab.ratio_err[m], fmt=mk, color=col,
                     mfc=col, mew=1.4)
        pt = tab[tab.partial & tab.ratio_sim_data.notna()]
        if len(pt):
            axr.errorbar(pt.V_gem, pt.ratio_sim_data, yerr=pt.ratio_err, fmt=mk, color=col, mew=1.4,
                         fillstyle="left", mfc=col, mfcalt="white")

    ax.set_yscale("log")
    lo, hi = min(ys), max(ys)
    ax.set_ylim(lo / 2, hi * 60)              # room for the label and legends above the data
    ax.set_ylabel("Effective gain")
    axr.set_xlabel(r"$V_\mathrm{GEM}$ (V)")
    axr.set_ylabel(r"Sim./data")
    top = max(c["tab"].ratio_sim_data.max() + c["tab"].ratio_err.max() for c in curves)
    axr.set_ylim(0, max(0.25, 1.3 * top))
    if corrected:
        axr.axhline(1.0, color="0.5", lw=0.8, ls=":", zorder=1)

    if gs:
        gs.add_experiment_label(ax, "Work In Progress")
        chem = {"CH4": "CH$_4$", "CO2": "CO$_2$", "CF4": "CF$_4$"}
        gas = ":".join(chem.get(g, g) for g in curves[0]["gas"].split("/"))   # Ar:CH4 92:8
        note = f", ADC offset {adc_offset:g} ch corrected" if adc_offset else ""
        ax.text(0.03, 0.865, tex(f"{gas} {curves[0]['ratio'].replace('/', ':')}{note}"), transform=ax.transAxes,
                ha="left", va="top", fontsize="small")
    grey = "0.3"
    what = [Line2D([], [], ls="-", lw=1.2, marker="o", mfc="white", mew=1.4, color=grey, label="GORG"),
            Line2D([], [], ls="--", lw=1.2, marker="o", mfc=grey, mew=1.4, color=grey,
                   label=tex("SimGEM, corrected") if corrected else "SimGEM")]
    if any(c["sim_all"].partial.any() for c in curves):
        what.append(Line2D([], [], ls="none", marker="o", mew=1.4, color=grey, fillstyle="left",
                           mfc=grey, mfcalt="white", label=tex("SimGEM, unfinished")))
    where = []
    for k, c in enumerate(curves):
        col, mk = pstyle(k, c["psi"])
        where.append(Line2D([], [], ls="none", marker=mk, ms=8, color=col, mfc=col, label=f"{c['psi']:g} psi"))
    leg = ax.legend(handles=what, loc="upper left", bbox_to_anchor=(0.0, 0.80), fontsize="small")
    ax.add_artist(leg)
    ax.legend(handles=where, loc="upper right", fontsize="small")
    if gs:
        gs.savefig(fig, out)
    else:
        fig.savefig(out + ".png", dpi=200, bbox_inches="tight")
    print(f"\n-> {out}.png" + (f", {out}.pdf" if gs else "") + f", {out}.csv")


if __name__ == "__main__":
    main()
