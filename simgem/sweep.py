"""
Parameter scans driven by scan files (scans/*.json).

A scan file names a base config and the parameters to vary; every point is
an ordinary config (base + overrides) with its own result folder. Example:

{
  "description": "GEM voltage x pressure, transfer fields scaled with pressure",
  "base": "../configs/baseline.json",          // relative to the scan file
  "set":  {"simulation.avalanche_size_limit": 2000000},   // fixed, every point
  "mode": "grid",                              // "grid" (all combinations) or "zip" (paired)
  "axes": [
    {"name": "dV", "path": "foils[*].dV_V", "values": [520, 540, 560]},
    {"name": "p",  "path": "gas.pressure_bar", "range": {"start": 6, "stop": 10, "step": 2}},
    {"name": "mix", "path": "gas.components",
       "values": [{"ar": 98, "ch4": 2}, {"ar": 95, "ch4": 5}], "labels": ["CH4_2", "CH4_5"]},
    {"name": "split",                          // several paths moving together = one axis
       "set": {"foils[0].dV_V": [560, 540], "foils[2].dV_V": [520, 540]}}
  ],
  "derived": {                                 // computed per point from the axes
    "gaps[1].field_kV_cm": "base('gaps[1].field_kV_cm') * p / 8"
  },
  "points": [{"label": "nominal", "set": {}}], // optional extra points, taken literally
                                               // ("set" + the fixed "set"; no "derived")
  "events": 500,                               // per point (default: base simulation.events)
  "analysis": {"polya": true, "column": "effective_gain"}   // "threshold": n for a plain
}                                                             // Polya above n (default: none)

Ranges: {"start", "stop", "step"} (inclusive) or {"start", "stop", "num"},
optionally with "log": true (geometric spacing). In "derived" expressions the
axis names stand for their values (their index if a value is not a number),
`base('path')` reads the base config, and math functions are available.

Running is resumable: field maps come from the cache, and points that already
have results for the same run key are not rerun. A failing point is recorded
and the scan goes on (unless --stop-on-error).
"""
import ast
import copy
import csv
import itertools
import json
import math
import os
import re
import time
import traceback

from . import analysis, corrections, field, runner
from . import config as C
from .util import ncpu

SPEC_KEYS = {"description", "base", "set", "mode", "axes", "derived", "points", "events",
             "seed", "jobs", "analysis", "output"}
ANALYSIS_DEFAULTS = {"polya": True, "threshold": None, "column": "effective_gain",
                     "max_value": None, "bin_width": None}


class ScanError(ValueError):
    pass


# ----------------------------------------------------------------- safe expressions
_FUNCS = {k: getattr(math, k) for k in ("sqrt", "exp", "log", "log10", "sin", "cos", "tan",
                                        "pi", "e", "floor", "ceil")}
_FUNCS.update({"abs": abs, "min": min, "max": max, "round": round})


def safe_eval(expr, names):
    """Evaluate an arithmetic expression with the given names; no attribute
    access, no imports, only whitelisted functions."""
    tree = ast.parse(expr, mode="eval")
    ok = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load,
          ast.Call, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod, ast.FloorDiv,
          ast.USub, ast.UAdd, ast.IfExp, ast.Compare, ast.Lt, ast.Gt, ast.LtE, ast.GtE,
          ast.Eq, ast.NotEq, ast.BoolOp, ast.And, ast.Or)
    for node in ast.walk(tree):
        if not isinstance(node, ok):
            raise ScanError(f"expression '{expr}': '{type(node).__name__}' is not allowed")
        if isinstance(node, ast.Call) and not (isinstance(node.func, ast.Name)
                                               and (node.func.id in _FUNCS or node.func.id == "base")):
            raise ScanError(f"expression '{expr}': only math functions and base() can be called")
    env = dict(_FUNCS)
    env.update(names)
    return eval(compile(tree, "<expr>", "eval"), {"__builtins__": {}}, env)


# ----------------------------------------------------------------- spec -> points
def _expand_range(r):
    if "num" in r:
        a, b, n = r["start"], r["stop"], int(r["num"])
        if n < 1:
            raise ScanError("range num must be >= 1")
        if r.get("log"):
            if a <= 0 or b <= 0:
                raise ScanError("log range needs positive start/stop")
            return [a * (b / a) ** (i / (n - 1)) if n > 1 else a for i in range(n)]
        return [a + (b - a) * i / (n - 1) if n > 1 else a for i in range(n)]
    a, b, s = r["start"], r["stop"], r["step"]
    if s == 0 or (b - a) / s < 0:
        raise ScanError(f"range {r}: step has the wrong sign or is 0")
    n = int(math.floor((b - a) / s + 1e-9)) + 1
    vals = [a + i * s for i in range(n)]
    return [round(v, 12) for v in vals]


def _fmt(v):
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        return f"{v:g}"
    if isinstance(v, (int, str)) or v is None:
        return str(v)
    if isinstance(v, dict):
        return "_".join(f"{k}{_fmt(x)}" for k, x in v.items())
    if isinstance(v, list):
        return "-".join(_fmt(x) for x in v)
    return str(v)


def _slug(s):
    return re.sub(r"[^A-Za-z0-9.=+-]+", "_", s).strip("_")[:80]


def parse_axes(spec):
    """-> list of axes: {name, labels[i], values[i], overrides[i] -> [(path, value)]}"""
    axes = []
    for k, a in enumerate(spec.get("axes", [])):
        name = a.get("name") or f"axis{k}"
        if not re.fullmatch(r"[A-Za-z_]\w*", name):
            raise ScanError(f"axis name '{name}' must be an identifier (used in 'derived')")
        if "path" in a:
            if ("values" in a) == ("range" in a):
                raise ScanError(f"axis '{name}': give exactly one of 'values' or 'range'")
            vals = a["values"] if "values" in a else _expand_range(a["range"])
            ovs = [[(a["path"], v)] for v in vals]
        elif "set" in a:
            lens = {len(v) for v in a["set"].values()}
            if len(lens) != 1:
                raise ScanError(f"axis '{name}': all lists in 'set' must have the same length")
            n = lens.pop()
            vals = [{p: a["set"][p][i] for p in a["set"]} for i in range(n)]
            ovs = [[(p, a["set"][p][i]) for p in a["set"]] for i in range(n)]
        else:
            raise ScanError(f"axis '{name}': needs 'path' + 'values'/'range', or 'set'")
        labels = a.get("labels") or [_fmt(v) if "path" in a else
                                     ",".join(_fmt(x) for x in v.values()) for v in vals]
        if len(labels) != len(vals):
            raise ScanError(f"axis '{name}': {len(labels)} labels for {len(vals)} values")
        if not vals:
            raise ScanError(f"axis '{name}' has no values")
        axes.append({"name": name, "path": a.get("path"), "values": vals, "labels": labels,
                     "overrides": ovs})
    if len({a["name"] for a in axes}) != len(axes):
        raise ScanError("axis names must be unique")
    return axes


def expand(spec, base_raw):
    """All points of a scan: list of dicts with index, label, coords, overrides."""
    axes = parse_axes(spec)
    mode = spec.get("mode", "grid")
    if mode not in ("grid", "zip"):
        raise ScanError("mode must be 'grid' or 'zip'")
    if mode == "zip" and len({len(a["values"]) for a in axes}) > 1:
        raise ScanError("mode 'zip' needs all axes to have the same number of values")
    fixed = list((spec.get("set") or {}).items())
    derived = spec.get("derived") or {}

    combos = []
    if axes:
        idx = (itertools.product(*[range(len(a["values"])) for a in axes]) if mode == "grid"
               else [(i,) * len(axes) for i in range(len(axes[0]["values"]))])
        combos = list(idx)
    points = []
    for combo in combos:
        coords, ovs, names = {}, list(fixed), {}
        for a, i in zip(axes, combo):
            coords[a["name"]] = a["labels"][i]
            ovs += a["overrides"][i]
            v = a["values"][i]
            names[a["name"]] = v if isinstance(v, (int, float)) and not isinstance(v, bool) else i
        points.append({"coords": coords, "names": names, "overrides": ovs,
                       "label": "_".join(f"{n}={l}" for n, l in coords.items())})
    for p in spec.get("points") or []:
        if "label" not in p:
            raise ScanError("every explicit point needs a 'label'")
        ovs = list(fixed) + [(k, v) for k, v in (p.get("set") or {}).items()]
        points.append({"coords": {"point": p["label"]}, "names": {}, "overrides": ovs,
                       "label": p["label"], "explicit": True})
    if not points:
        raise ScanError("the scan has no points (give 'axes' and/or 'points')")

    def base_get(path):
        return C.get_path(base_raw, path)

    for k, p in enumerate(points):
        p["index"] = k
        names = dict(p["names"], base=base_get)
        for path, expr in (derived.items() if not p.get("explicit") else ()):
            try:
                val = safe_eval(expr, names) if isinstance(expr, str) else expr
            except ScanError:
                raise
            except Exception as e:
                raise ScanError(f"derived '{path}' = '{expr}' failed at point {p['label']}: {e}")
            p["overrides"].append((path, val))
            p.setdefault("derived", {})[path] = val
        p["dir"] = f"p{k:03d}_{_slug(p['label'])}"
    return axes, points


# ----------------------------------------------------------------- loading
def load_spec(path):
    with open(path) as f:
        spec = json.load(f)
    unknown = set(spec) - SPEC_KEYS
    if unknown:
        raise ScanError(f"unknown key(s) in scan file: {sorted(unknown)}")
    if "base" not in spec:
        raise ScanError("scan file needs 'base' (path of the base config)")
    base = spec["base"]
    if not os.path.isabs(base):
        cand = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), base))
        base = cand if os.path.exists(cand) else base
    spec["base"] = base
    spec["_file"] = os.path.abspath(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
    # scans/gorg/x.json -> output runs/scans/gorg/x
    spec["_name"] = stem if parent == "scans" else f"{parent}/{stem}"
    return spec


def is_scan_file(path):
    try:
        with open(path) as f:
            d = json.load(f)
    except (OSError, json.JSONDecodeError):
        return False
    return isinstance(d, dict) and "base" in d and ("axes" in d or "points" in d)


def parse_only(text, n):
    """'0,3-5,9' -> sorted set of indices."""
    out = set()
    for part in text.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out |= set(range(int(a), int(b) + 1))
        elif part:
            out.add(int(part))
    bad = [i for i in out if i < 0 or i >= n]
    if bad:
        raise ScanError(f"--only: no point(s) {bad} (scan has {n})")
    return sorted(out)


def resolve_points(spec, extra_overrides=()):
    """Expand and resolve every point to a full config; collects all errors."""
    with open(spec["base"]) as f:
        base_raw = json.load(f)
    for ov in extra_overrides:
        C.set_path(base_raw, *C.parse_override(ov))
    axes, points = expand(spec, base_raw)
    errors = []
    for p in points:
        raw = copy.deepcopy(base_raw)
        try:
            for path, v in p["overrides"]:
                C.set_path(raw, path, v)
            p["cfg"] = C.resolve(raw)
            p["stack"] = C.Stack(p["cfg"])
            p["field_hash"] = p["stack"].field_hash()
            p["events"] = spec.get("events") or p["cfg"]["simulation"]["events"]
        except C.ConfigError as e:
            errors.append(f"point {p['index']} ({p['label']}): {e}")
    return axes, points, errors


# ----------------------------------------------------------------- running
def _result_row(axes, p, s, fit, status, err=""):
    row = {"index": p["index"], "label": p["label"], "dir": p["dir"]}
    for a in axes:
        row[a["name"]] = p["coords"].get(a["name"], "")
    for path, v in (p.get("derived") or {}).items():
        row[path] = v
    row.update({"status": status, "events": (s or {}).get("events", ""),
                "effective_gain": (s or {}).get("effective_gain", ""),
                "effective_gain_err": (s or {}).get("effective_gain_err", ""),
                "absolute_gain": (s or {}).get("absolute_gain", ""),
                "absolute_gain_err": (s or {}).get("absolute_gain_err", ""),
                "zero_fraction": (s or {}).get("zero_fraction", ""),
                "truncated_events": (s or {}).get("truncated_events", "")})
    if fit is not None:
        row.update({"gain_fit": fit.get("gain", ""), "gain_fit_err": fit.get("gain_err", ""),
                    "survival": fit.get("survival", ""),
                    "polya_G0": fit.get("G0", ""), "polya_G0_err": fit.get("G0_err", ""),
                    "polya_theta": fit.get("theta", ""), "polya_theta_err": fit.get("theta_err", ""),
                    "polya_chi2_ndof": (f"{fit['chi2']:.1f}/{fit['ndof']}" if fit.get("ok") else ""),
                    "above_threshold": fit.get("efficiency", "")})
    k = corrections.factor(p["cfg"]) if p.get("cfg") else 1.0
    if k != 1.0 and s:
        g, ge = ((fit["gain"], fit["gain_err"]) if fit is not None and fit.get("ok") and "gain" in fit
                 else (s["effective_gain"], s["effective_gain_err"]))
        row.update({"gain_correction": k, "gain_corrected": k * g, "gain_corrected_err": k * ge})
    row.update({"field_hash": p.get("field_hash", ""), "error": err})
    return row


def _cached(pdir, key, events):
    f = os.path.join(pdir, "gain_summary.json")
    if os.path.exists(f):
        with open(f) as fh:
            s = json.load(fh)
        if s.get("run_key") == key and s.get("events") == events:
            return s
    return None


def run(spec, outdir=None, extra_overrides=(), events=None, jobs=None, only=None, dry_run=False,
        force=False, analyze_only=False, stop_on_error=False, plot=True):
    if events:
        spec = dict(spec, events=events)
    axes, points, errors = resolve_points(spec, extra_overrides)
    outdir = outdir or spec.get("output") or os.path.join(
        os.environ.get("SIMGEM_RUNS", "runs"), "scans", spec["_name"])
    outdir = os.path.abspath(outdir)
    ana = dict(ANALYSIS_DEFAULTS, **(spec.get("analysis") or {}))
    jobs = jobs or spec.get("jobs") or ncpu()
    sel = parse_only(only, len(points)) if only else list(range(len(points)))

    # ---- plan -------------------------------------------------------------
    print(f"Scan '{spec['_name']}': {len(points)} point(s), base {spec['base']}")
    if spec.get("description"):
        print(f"  {spec['description']}")
    if errors:
        raise ScanError("invalid point(s):\n  " + "\n  ".join(errors))
    need_fields = []
    w = max(len(p["label"]) for p in points)
    print(f"  {'#':>4}  {'point':<{w}}  {'field':<12} {'events':>7}  status")
    for i in sel:
        p = points[i]
        h = p["field_hash"]
        fstate = "cached" if field.is_ready(h) else "NEW"
        if fstate == "NEW" and h not in need_fields:
            need_fields.append(h)
        key = runner.run_key(p["cfg"], h, p["cfg"]["simulation"]["seed"])
        done = _cached(os.path.join(outdir, p["dir"]), key, p["events"]) is not None
        print(f"  {i:>4}  {p['label']:<{w}}  {h} {p['events']:>7}  "
              f"{'done' if done and not force else 'to run'}{'' if fstate == 'cached' else ' (new field)'}")
    n_run = sum(1 for i in sel if force or _cached(
        os.path.join(outdir, points[i]["dir"]),
        runner.run_key(points[i]["cfg"], points[i]["field_hash"], points[i]["cfg"]["simulation"]["seed"]),
        points[i]["events"]) is None)
    print(f"  -> {len(need_fields)} field map(s) to build, {n_run} point(s) to simulate, "
          f"{jobs} workers, output {outdir}")
    if dry_run:
        return None
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "scan.json"), "w") as f:
        json.dump({"spec": {k: v for k, v in spec.items() if not k.startswith("_")},
                   "spec_file": spec["_file"], "extra_overrides": list(extra_overrides),
                   "points": [{"index": p["index"], "label": p["label"], "dir": p["dir"],
                               "coords": p["coords"], "overrides": p["overrides"],
                               "field_hash": p["field_hash"], "events": p["events"]}
                              for p in points]}, f, indent=2, default=str)

    # ---- fields first (each built once; later points then only simulate) --
    failed_fields = {}
    if not analyze_only:
        for k, h in enumerate(need_fields):
            p = next(points[i] for i in sel if points[i]["field_hash"] == h)
            print(f"\n=== field {k + 1}/{len(need_fields)}: {h} (first used by '{p['label']}')")
            try:
                field.build(p["stack"], h=h)
            except Exception as e:      # noqa: BLE001 - record and go on
                failed_fields[h] = str(e)
                print(f"  FIELD FAILED: {e}")
                if stop_on_error:
                    raise

    # ---- points ------------------------------------------------------------
    rows = []
    t0 = time.time()
    selected = set(sel)
    n = -1
    for i, p in enumerate(points):
        pdir = os.path.join(outdir, p["dir"])
        h = p["field_hash"]
        cfg = p["cfg"]
        key = runner.run_key(cfg, h, cfg["simulation"]["seed"])
        if i not in selected:           # keep the table complete: existing results only
            s = _cached(pdir, key, p["events"])
            fit = None
            fp = os.path.join(pdir, "gain_fit.json")
            if s and os.path.exists(fp):
                with open(fp) as f:
                    fit = json.load(f)
            rows.append(_result_row(axes, p, s, fit, "cached" if s else "not run"))
            continue
        n += 1
        print(f"\n=== point {i} ({n + 1}/{len(sel)}): {p['label']}")
        s, fit, status, err = None, None, "", ""
        try:
            s = None if (force and not analyze_only) else _cached(pdir, key, p["events"])
            if s:
                status = "cached"
                print("  (result exists)")
            elif analyze_only:
                status = "missing"
            elif h in failed_fields:
                status, err = "field-failed", failed_fields[h]
            else:
                os.makedirs(pdir, exist_ok=True)
                with open(os.path.join(pdir, "config.json"), "w") as f:
                    json.dump({"scan": spec["_file"], "point": p["label"],
                               "overrides": p["overrides"], "resolved": cfg}, f, indent=2, default=str)
                s = runner.run_gain(cfg, p["stack"], h, pdir, p["events"], jobs=jobs)
                status = "done"
            if s:
                print(runner.format_summary(s, p["label"]))
                if ana["polya"]:
                    fit = analysis.analyze_run(os.path.join(pdir, "gain_events.csv"), pdir,
                                               column=ana["column"], threshold=ana["threshold"],
                                               max_value=ana["max_value"], bin_width=ana["bin_width"],
                                               plot=plot, quiet=True)
                    if fit.get("ok"):
                        print((f"  fitted gain = {fit['gain']:.4g} +/- {fit['gain_err']:.2g}, "
                               if "gain" in fit else "") +
                              f"Polya G0 = {fit['G0']:.4g} +/- {fit['G0_err']:.2g}, "
                              f"theta = {fit['theta']:.3g} +/- {fit['theta_err']:.2g}")
        except KeyboardInterrupt:
            raise
        except Exception as e:          # noqa: BLE001 - one bad point must not kill the scan
            status, err = "failed", str(e).splitlines()[0] if str(e) else type(e).__name__
            print(f"  POINT FAILED: {e}")
            with open(os.path.join(outdir, f"{p['dir']}.error.txt"), "w") as f:
                f.write(traceback.format_exc())
            if stop_on_error:
                raise
        rows.append(_result_row(axes, p, s, fit, status, err))
        _write_table(outdir, rows)
        if n + 1 < len(sel) and status == "done":
            el = time.time() - t0
            print(f"  [{n + 1}/{len(sel)} points, {el / 60:.1f} min elapsed]")

    _write_table(outdir, rows)
    print_table(axes, rows)
    if plot:
        png = plot_scan(axes, rows, os.path.join(outdir, "scan.png"), ana)
        if png:
            print(f"Plot  -> {png}")
    print(f"Table -> {os.path.join(outdir, 'scan.csv')}")
    bad = [r for r in rows if r["status"] in ("failed", "field-failed")]
    if bad:
        print(f"WARNING: {len(bad)} point(s) failed: " + ", ".join(str(r["index"]) for r in bad))
    return rows


def _write_table(outdir, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(outdir, "scan.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def print_table(axes, rows):
    names = [a["name"] for a in axes] or ["label"]
    print("\n" + "  ".join(f"{n:>12}" for n in names) +
          f"  {'mean gain':>11} {'+/-':>8}  {'fit gain':>10} {'+/-':>8}  status")
    for r in rows:
        cols = [str(r.get(n, r["label"]))[:12] for n in names]
        g = r.get("effective_gain")
        G = r.get("gain_fit")
        print("  ".join(f"{c:>12}" for c in cols) +
              (f"  {g:>11.4g} {r['effective_gain_err']:>8.2g}" if g != "" else f"  {'-':>11} {'':>8}") +
              (f"  {G:>10.4g} {r['gain_fit_err']:>8.2g}" if G not in ("", None) else f"  {'-':>10} {'':>8}") +
              f"  {r['status']}")


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def plot_scan(axes, rows, path, ana):
    """Gain vs the first axis (one line per value of the second axis; further
    axes are folded into the line label). Needs matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    ok = [r for r in rows if r.get("effective_gain") not in ("", None)]
    if not ok:
        return None
    fig, ax = plt.subplots(figsize=(7.5, 5))
    if axes:
        x0 = axes[0]["name"]
        numeric = all(_num(r[x0]) is not None for r in ok)
        groups = {}
        for r in ok:
            key = ", ".join(f"{a['name']}={r[a['name']]}" for a in axes[1:]) or None
            groups.setdefault(key, []).append(r)
        xticks = None
        if not numeric:
            cats = list(dict.fromkeys(r[x0] for r in ok))
            xticks = cats
        for k, (key, rs) in enumerate(groups.items()):
            xs = [(_num(r[x0]) if numeric else xticks.index(r[x0])) for r in rs]
            order = sorted(range(len(xs)), key=lambda i: xs[i])
            xs = [xs[i] for i in order]
            rs = [rs[i] for i in order]
            c = f"C{k % 10}"
            ax.errorbar(xs, [r["effective_gain"] for r in rs], yerr=[r["effective_gain_err"] for r in rs],
                        fmt="o-", color=c, ms=4, lw=1, capsize=2,
                        label=("mean (all events)" + (f", {key}" if key else "")))
            G = [(x, r) for x, r in zip(xs, rs) if _num(r.get("gain_fit")) is not None]
            if G:
                ax.errorbar([x for x, _ in G], [r["gain_fit"] for _, r in G],
                            yerr=[r["gain_fit_err"] for _, r in G], fmt="s--", color=c, ms=3,
                            lw=0.8, mfc="none", capsize=2,
                            label=("fitted gain (zero-inflated Polya)" + (f", {key}" if key else "")))
        ax.set_xlabel(x0 + (f"  ({axes[0]['path']})" if axes[0].get("path") else ""))
        if xticks:
            ax.set_xticks(range(len(xticks)))
            ax.set_xticklabels(xticks, rotation=30, ha="right")
    else:
        xs = list(range(len(ok)))
        ax.errorbar(xs, [r["effective_gain"] for r in ok], yerr=[r["effective_gain_err"] for r in ok],
                    fmt="o", label="mean (all events)")
        ax.set_xticks(xs)
        ax.set_xticklabels([r["label"] for r in ok], rotation=30, ha="right")
    vals = [r["effective_gain"] for r in ok if r["effective_gain"] > 0]
    if vals and max(vals) / min(vals) > 20:
        ax.set_yscale("log")
    ax.set_ylabel(ana["column"].replace("_", " "))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def adhoc_spec(config_path, param, values, outdir):
    """`gem.py scan CONFIG --param P --values ...` as a one-axis scan spec."""
    name = "scan_" + _slug(param.replace("[", "_").replace("]", "").replace("*", "all"))
    return {"base": os.path.abspath(config_path), "axes": [{"name": "x", "path": param,
            "values": [C.parse_override(f"x={v}")[1] for v in values]}],
            "_file": os.path.abspath(config_path), "_name": name,
            "output": os.path.join(outdir, name)}
