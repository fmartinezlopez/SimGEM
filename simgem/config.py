"""
Configuration: loading, defaults, validation, command-line overrides and the
quantities derived from a config (electrode heights and potentials, the
periodic cell, the field-map hash).

Conventions
-----------
* Foils and gaps are listed TOP (drift side) -> BOTTOM (readout side).
  There is one more gap than foils: drift, transfer..., induction.
* z = 0 at the readout plane, increasing upward. All lengths in a config carry
  their unit in the key name (_um, _mm); internally everything is in mm.
* Potentials decrease upward (readout = 0 V, drift cathode most negative), so
  electrons drift down to the readout.
* A gap thickness is the distance between the outer copper surfaces bounding
  it. The drift gap is the *modeled* slice above the top foil, carrying the
  true drift field (the rest of a real drift volume is uniform and adds nothing
  to the gain).
"""
import copy
import hashlib
import json
import math
import re
from pathlib import Path

from . import corrections

# Bump when the geometry/field construction changes in a way that makes old
# cached field maps invalid.
FIELD_MODEL_VERSION = 2

FOIL_DEFAULTS = {
    "pattern": "hex",               # "hex" (standard GEM) or "square"
    "pitch_um": 140.0,
    "kapton_thickness_um": 50.0,
    "copper_thickness_um": 5.0,     # 0 -> zero-thickness electrodes on the kapton
    "kapton_permittivity": 3.5,
    "offset_um": [0.0, 0.0],        # lateral shift of this foil's hole pattern
    "hole": {
        # "double_conical": outer diameter at both kapton faces, inner at mid-plane
        # "conical":        outer diameter on `wide_side`, inner on the other face
        # "cylindrical":    outer diameter throughout
        "profile": "double_conical",
        "outer_diameter_um": 70.0,
        "inner_diameter_um": 50.0,
        "wide_side": "top",
        "copper_diameter_um": None,  # None -> same as the kapton hole at that face
    },
}

DEFAULTS = {
    "description": "",
    "notes": [],                     # free text (assumptions, sources); not used
    "gas": {
        "components": None,          # required, e.g. {"ar": 98, "ch4": 2}
        "pressure_bar": 1.01325,
        "temperature_K": 293.15,
        # r: null (off), "auto" (Garfield's parameterisation) or a number.
        # match: null, or "table" = run the microscopic avalanche with the r that
        #   reproduces the Townsend coefficient of the table formula
        #   alpha (1 + r nu_exc / nu_ion), with which published r values are
        #   defined (Sahin et al.); see simgem/penning.py. Calibrated once per
        #   gas at match_fields_kV_cm (null -> 30/45/60/80 kV/cm x p/1 atm).
        "penning": {"r": None, "lambda": 0.0, "gas": None,
                    "match": None, "match_fields_kV_cm": None, "match_events": 400},
        # Optional Magboltz table (`gem.py gastable`); NOT used by the
        # microscopic avalanche, which works from the cross sections directly.
        "table": {"field_min_V_cm": 100.0, "field_max_V_cm": 1.0e6,
                  "field_points": 30, "max_electron_energy_eV": 400.0,
                  "collisions_1e7": 10},
    },
    "foil": FOIL_DEFAULTS,
    "foils": None,                   # required
    "gaps": None,                    # required
    "mesh": {"fine_um": 6.0, "coarse_um": 60.0,
             "fine_distance_um": None,     # None -> inner hole diameter
             "coarse_distance_um": None},  # None -> 3 x inner hole diameter
    "simulation": {
        "events": 1000,
        "avalanche_size_limit": 500000,
        "start_depth_um": None,      # below the drift plane; None -> middle of the drift gap
        "readout_threshold_um": 10.0,
        "primary_energy_eV": 0.1,
        "seed": 1,
        "stall_timeout_s": None,     # None -> 120 s without a collision = hung worker
    },
    "misalignment_mc": None,
    "hv": None,
    # Optional empirical gain corrections (charging-up, per-foil scale), applied
    # to the reported gains only: null, "literature" or a mapping; see
    # simgem/corrections.py.
    "corrections": None,
}

MC_DEFAULTS = {
    "distribution": "gaussian",      # "gaussian" (sigma_um) or "uniform" (whole cell)
    "sigma_um": 20.0,
    "foils": None,                   # names of the foils to perturb; None -> all
    "geometries": 20,
    "events_per_geometry": 500,
    "reference_events": None,        # None -> geometries x events_per_geometry
    "seed": 1,
}

# Optional "hv" block: the stack is powered like a real detector, by a GEM supply
# V_gem across a resistive divider (GEM1 top -> readout) plus a drift supply.
# Foil voltages and gap fields are then DERIVED and must not be given by hand.
#   "hv": {
#     "V_gem": 2600,                 # V across the divider = |potential of GEM1 top|
#     "V_drift_gap": 500,            # V_drift - V_gem, across the TRUE drift gap
#                                    #   (or "V_drift": absolute drift-electrode voltage)
#     "drift_length_mm": 18.0,       # true drift gap; drift field = V_drift_gap / length
#     "divider": {"gem1": 0.109, "transfer1": 0.2383, ...}   # share of V_gem per foil
#                                    #   and per gap below the top foil (names); sums to 1
#     or "divider_resistance_kohm": {"gem1": 550, "transfer1": 1000, ...}
#                                    #   resistance of each element; shares = R / sum(R)
#   }                                #   (unloaded divider: the GEM electrodes draw no current)
HV_KEYS = {"V_gem", "V_drift_gap", "V_drift", "drift_length_mm", "divider",
           "divider_resistance_kohm", "shares"}

PATTERNS = ("hex", "square")
PROFILES = ("double_conical", "conical", "cylindrical")


class ConfigError(ValueError):
    pass


# ----------------------------------------------------------------- dict helpers
def deep_merge(base, over):
    """Recursively merge `over` into a copy of `base` (dicts only; lists and
    scalars in `over` replace)."""
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _tokens(path):
    """'foils[1].dV_V' / 'foils.1.dV_V' / 'foils[*].dV_V' -> ['foils', 1, 'dV_V']."""
    toks = []
    for part in re.split(r"\.(?![^\[]*\])", path):
        m = re.fullmatch(r"([^\[\]]*)((?:\[[^\]]+\])*)", part)
        if not m:
            raise ConfigError(f"bad path '{path}'")
        if m.group(1):
            toks.append(int(m.group(1)) if m.group(1).isdigit() else m.group(1))
        for idx in re.findall(r"\[([^\]]+)\]", m.group(2)):
            toks.append("*" if idx == "*" else int(idx))
    return toks


def set_path(cfg, path, value):
    """Set a (possibly nested / indexed / wildcarded) key in a raw config."""
    def rec(node, toks):
        t, rest = toks[0], toks[1:]
        keys = range(len(node)) if t == "*" else [t]
        for k in keys:
            if not rest:
                node[k] = copy.deepcopy(value)
                continue
            if isinstance(node, dict) and k not in node:
                node[k] = {}
            try:
                rec(node[k], rest)
            except (KeyError, IndexError, TypeError):
                raise ConfigError(f"cannot set '{path}': no element '{k}'")
    toks = _tokens(path)
    if not toks:
        raise ConfigError("empty override path")
    rec(cfg, toks)


def get_path(cfg, path):
    """Read a nested / indexed key ('foils[1].dV_V'); '*' returns a list."""
    def rec(node, toks):
        if not toks:
            return node
        t, rest = toks[0], toks[1:]
        if t == "*":
            return [rec(x, rest) for x in node]
        try:
            return rec(node[t], rest)
        except (KeyError, IndexError, TypeError):
            raise ConfigError(f"no value at '{path}'")
    return rec(cfg, _tokens(path))


def parse_override(text):
    """'key.path=value' with value parsed as JSON when possible."""
    if "=" not in text:
        raise ConfigError(f"override '{text}' must look like path=value")
    path, raw = text.split("=", 1)
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        value = raw
    return path.strip(), value


# ----------------------------------------------------------------- loading
def load(path, overrides=()):
    """Load a config file, apply overrides, fill defaults and validate.
    Returns the fully resolved config dict."""
    with open(path) as f:
        raw = json.load(f)
    for ov in overrides:
        set_path(raw, *parse_override(ov) if isinstance(ov, str) else ov)
    return resolve(raw)


def resolve(raw):
    unknown = set(raw) - set(DEFAULTS)
    if unknown:
        raise ConfigError(f"unknown top-level key(s): {sorted(unknown)}")
    cfg = deep_merge(DEFAULTS, {k: v for k, v in raw.items() if k != "foils"})
    if not raw.get("foils"):
        raise ConfigError("'foils' must be a non-empty list")
    foils = []
    for i, f in enumerate(raw["foils"]):
        unknown = set(f) - set(FOIL_DEFAULTS) - {"name", "dV_V"}
        if unknown:
            raise ConfigError(f"foil {i}: unknown key(s) {sorted(unknown)}")
        rf = deep_merge(cfg["foil"], f)
        rf.setdefault("name", f"gem{i + 1}")
        foils.append(rf)
    cfg["foils"] = foils
    if cfg.get("hv") is not None:
        apply_hv(cfg)
    if cfg.get("misalignment_mc") is not None:
        cfg["misalignment_mc"] = deep_merge(MC_DEFAULTS, cfg["misalignment_mc"])
    try:
        cfg["corrections"] = corrections.resolve(raw.get("corrections"))
    except ValueError as e:
        raise ConfigError(str(e))
    validate(cfg)
    return cfg


def _default_gap_names(gaps):
    for i, gp in enumerate(gaps):
        gp.setdefault("name", ("drift" if i == 0 else "induction" if i == len(gaps) - 1
                               else f"transfer{i}"))


def apply_hv(cfg):
    """Derive foil voltages and gap fields from the 'hv' block."""
    hv = cfg["hv"]
    unknown = set(hv) - HV_KEYS
    _req(not unknown, f"hv: unknown key(s) {sorted(unknown)}")
    for k in ("V_gem", "drift_length_mm"):
        _req(k in hv, f"hv: missing '{k}'")
    _req(("divider" in hv) != ("divider_resistance_kohm" in hv),
         "hv: give exactly one of 'divider' (shares) or 'divider_resistance_kohm'")
    _req(("V_drift_gap" in hv) != ("V_drift" in hv), "hv: give exactly one of 'V_drift_gap' or 'V_drift'")
    gaps = cfg["gaps"]
    _req(isinstance(gaps, list) and len(gaps) == len(cfg["foils"]) + 1,
         f"'gaps' must list {len(cfg['foils']) + 1} gaps (drift, transfer..., induction)")
    _default_gap_names(gaps)
    Vg = float(hv["V_gem"])
    Vd = float(hv["V_drift_gap"]) if "V_drift_gap" in hv else float(hv["V_drift"]) - Vg
    _req(Vg > 0 and Vd > 0 and hv["drift_length_mm"] > 0,
         "hv: V_gem, the drift-gap voltage and drift_length_mm must be positive")
    names = [f["name"] for f in cfg["foils"]] + [g["name"] for g in gaps[1:]]
    if "divider" in hv:
        div = hv["divider"]
        _req(set(div) == set(names), f"hv.divider must give a share for exactly {names}")
        tot = sum(div.values())
        _req(abs(tot - 1) < 5e-3, f"hv.divider shares sum to {tot:.4f}, not 1")
    else:
        div = hv["divider_resistance_kohm"]
        _req(set(div) == set(names), f"hv.divider_resistance_kohm must give a resistance for exactly {names}")
        _req(all(v > 0 for v in div.values()), "hv.divider_resistance_kohm: resistances must be positive")
        tot = sum(div.values())
    shares = {k: div[k] / tot for k in names}
    hv["shares"] = shares
    for f in cfg["foils"]:
        _req("dV_V" not in f, f"{f['name']}: dV_V is derived from 'hv'; remove it")
        f["dV_V"] = shares[f["name"]] * Vg
    for g in gaps:
        _req("field_kV_cm" not in g, f"gap '{g['name']}': field_kV_cm is derived from 'hv'; remove it")
    for g in gaps[1:]:
        g["field_kV_cm"] = shares[g["name"]] * Vg / (g["thickness_mm"] / 10.0) / 1e3
    gaps[0]["field_kV_cm"] = Vd / (hv["drift_length_mm"] / 10.0) / 1e3


def _req(cond, msg):
    if not cond:
        raise ConfigError(msg)


def validate(cfg):
    g = cfg["gas"]
    comps = g["components"]
    _req(isinstance(comps, dict) and comps, "gas.components must be a non-empty "
         "mapping, e.g. {\"ar\": 98, \"ch4\": 2}")
    _req(len(comps) <= 6, "Magboltz supports at most 6 gas components")
    _req(all(v > 0 for v in comps.values()), "gas fractions must be positive")
    _req(g["pressure_bar"] > 0 and g["temperature_K"] > 0, "pressure and temperature must be positive")
    r = g["penning"]["r"]
    _req(r is None or r == "auto" or (isinstance(r, (int, float)) and 0 <= r <= 1),
         "gas.penning.r must be null (off), \"auto\" or a number in [0, 1]")
    if r == "auto":
        _req(len(comps) == 2, "Penning r = \"auto\" (Garfield's parameterisation) "
             "only exists for binary mixtures; give an explicit r")
    pg = g["penning"]["gas"]
    _req(pg is None or pg.lower() in {c.lower() for c in comps},
         f"gas.penning.gas '{pg}' is not a component of the mixture")
    pm = g["penning"]
    _req(pm["match"] in (None, "table"), "gas.penning.match must be null or \"table\"")
    if pm["match"]:
        _req(r == "auto" or (r is not None and r > 0),
             "gas.penning.match needs Penning transfer on (r = \"auto\" or > 0)")
        fl = pm["match_fields_kV_cm"]
        _req(fl is None or (isinstance(fl, list) and fl and all(v > 0 for v in fl)),
             "gas.penning.match_fields_kV_cm must be null or a list of positive fields")
        _req(pm["match_events"] >= 50, "gas.penning.match_events must be >= 50")
    corr = cfg.get("corrections")
    if corr and corr["requires_penning_match"]:
        _req(pm["match"] == "table", "these gain corrections were calibrated with Penning "
             "matching; set gas.penning.match = \"table\" (or give explicit factors)")

    names = [f["name"] for f in cfg["foils"]]
    _req(len(set(names)) == len(names), "foil names must be unique")
    for f in cfg["foils"]:
        n = f["name"]
        _req(f["pattern"] in PATTERNS, f"{n}: pattern must be one of {PATTERNS}")
        _req("dV_V" in f, f"{n}: missing dV_V (voltage across the foil)")
        _req(f["pitch_um"] > 0 and f["kapton_thickness_um"] > 0, f"{n}: pitch and kapton thickness must be positive")
        _req(f["copper_thickness_um"] >= 0, f"{n}: copper thickness must be >= 0")
        _req(len(f["offset_um"]) == 2, f"{n}: offset_um must be [dx, dy]")
        h = f["hole"]
        _req(h["profile"] in PROFILES, f"{n}: hole.profile must be one of {PROFILES}")
        _req(h["wide_side"] in ("top", "bottom"), f"{n}: hole.wide_side must be 'top' or 'bottom'")
        D = h["outer_diameter_um"]
        d = h["inner_diameter_um"] if h["profile"] != "cylindrical" else D
        _req(0 < d <= D, f"{n}: need 0 < inner_diameter <= outer_diameter")
        _req(h["profile"] == "cylindrical" or d < D,
             f"{n}: inner == outer diameter; use profile 'cylindrical'")
        Dcu = h["copper_diameter_um"] or D
        _req(Dcu < f["pitch_um"] and D < f["pitch_um"], f"{n}: hole wider than the pitch")

    gaps = cfg["gaps"]
    _req(isinstance(gaps, list) and len(gaps) == len(cfg["foils"]) + 1,
         f"'gaps' must list {len(cfg['foils']) + 1} gaps (drift, transfer..., induction) "
         f"for {len(cfg['foils'])} foils")
    for i, gp in enumerate(gaps):
        _req(gp.get("thickness_mm", 0) > 0, f"gap {i}: thickness_mm must be positive")
        _req("field_kV_cm" in gp, f"gap {i}: missing field_kV_cm")
    _default_gap_names(gaps)

    m = cfg["mesh"]
    _req(0 < m["fine_um"] <= m["coarse_um"], "mesh: need 0 < fine_um <= coarse_um")

    s = cfg["simulation"]
    _req(s["events"] > 0 and s["avalanche_size_limit"] > 0, "simulation: events and size limit must be positive")
    st = Stack(cfg)
    if st.start_depth_mm >= cfg["gaps"][0]["thickness_mm"]:
        raise ConfigError("simulation.start_depth_um puts the primaries below the drift gap")

    mc = cfg.get("misalignment_mc")
    if mc:
        _req(mc["distribution"] in ("gaussian", "uniform"),
             "misalignment_mc.distribution must be 'gaussian' or 'uniform'")
        _req(mc["distribution"] == "uniform" or mc["sigma_um"] > 0, "misalignment_mc.sigma_um must be positive")
        for n in mc["foils"] or []:
            _req(n in names, f"misalignment_mc.foils: unknown foil '{n}'")
        _req(mc["geometries"] >= 2, "misalignment_mc.geometries must be >= 2")


# ----------------------------------------------------------------- lattice
def lattice(foil):
    """(a1, a2, cell_x, cell_y) in mm for a foil's hole pattern. The
    rectangular cell (cell_x, cell_y) is a translational period of the pattern:
    1 hole for square, 2 for hex (centre + 4 quarter-holes at the corners)."""
    p = foil["pitch_um"] / 1000.0
    if foil["pattern"] == "square":
        return (p, 0.0), (0.0, p), p, p
    return (p, 0.0), (p / 2, p * math.sqrt(3) / 2), p, p * math.sqrt(3)


def common_period(lengths, max_mult=6, tol=1e-6):
    """Smallest length that is an integer multiple of every entry (supercell
    for foils with commensurate but different patterns, e.g. 140/280 um)."""
    L0 = max(lengths)
    for k in range(1, max_mult + 1):
        L = k * L0
        if all(abs(L / l - round(L / l)) < tol for l in lengths):
            return L
    raise ConfigError(f"foil patterns are not commensurate (periods {lengths} mm); "
                      "a single periodic cell cannot hold them")


def min_reciprocal(foil):
    """Length (1/mm) of the shortest reciprocal-lattice vector of the pattern."""
    p = foil["pitch_um"] / 1000.0
    return 2 * math.pi / p if foil["pattern"] == "square" else 4 * math.pi / (math.sqrt(3) * p)


# ----------------------------------------------------------------- the stack
class Stack:
    """Geometry and electrostatics derived from a resolved config (mm, V)."""

    def __init__(self, cfg, offsets_mm=None):
        self.cfg = cfg
        self.foils = cfg["foils"]
        self.gaps = cfg["gaps"]
        cells = [lattice(f)[2:] for f in self.foils]
        self.cell = (common_period([c[0] for c in cells]), common_period([c[1] for c in cells]))
        # lateral offsets (mm): nominal from the config unless overridden
        self.offsets = {f["name"]: tuple(v / 1000.0 for v in f["offset_um"]) for f in self.foils}
        if offsets_mm:
            self.offsets.update(offsets_mm)
        self._heights()
        self._potentials()

    def foil_thickness(self, f):
        return (f["kapton_thickness_um"] + 2 * f["copper_thickness_um"]) / 1000.0

    def _heights(self):
        """Electrode planes, bottom-up. For foil f:
        z[f.name+'_bot'] / z[f.name+'_top'] = outer copper surfaces;
        kapton[f.name] = (z_kapton_bottom, z_kapton_top)."""
        z = {"readout": 0.0}
        kap = {}
        cur = 0.0
        for gap, f in zip(reversed(self.gaps[1:]), reversed(self.foils)):
            cur += gap["thickness_mm"]
            z[f["name"] + "_bot"] = cur
            tcu = f["copper_thickness_um"] / 1000.0
            tk = f["kapton_thickness_um"] / 1000.0
            kap[f["name"]] = (cur + tcu, cur + tcu + tk)
            cur += self.foil_thickness(f)
            z[f["name"] + "_top"] = cur
        cur += self.gaps[0]["thickness_mm"]
        z["drift"] = cur
        self.z = z
        self.kapton = kap
        self.height = cur

    def _potentials(self):
        V = {"readout": 0.0}
        cur = 0.0
        for gap, f in zip(reversed(self.gaps[1:]), reversed(self.foils)):
            cur -= gap["field_kV_cm"] * 1e3 * gap["thickness_mm"] / 10.0
            V[f["name"] + "_bot"] = cur
            cur -= f["dV_V"]
            V[f["name"] + "_top"] = cur
        cur -= self.gaps[0]["field_kV_cm"] * 1e3 * self.gaps[0]["thickness_mm"] / 10.0
        V["drift"] = cur
        self.V = V

    def electrodes(self):
        """Electrode names top -> bottom."""
        out = ["drift"]
        for f in self.foils:
            out += [f["name"] + "_top", f["name"] + "_bot"]
        return out + ["readout"]

    @property
    def start_depth_mm(self):
        d = self.cfg["simulation"]["start_depth_um"]
        return self.gaps[0]["thickness_mm"] / 2 if d is None else d / 1000.0

    @property
    def z_start(self):
        return self.z["drift"] - self.start_depth_mm

    @property
    def z_readout(self):
        return self.cfg["simulation"]["readout_threshold_um"] / 1000.0

    def field_key(self):
        """Everything that determines the field map, canonicalised."""
        foils = []
        for f in self.foils:
            ff = {k: f[k] for k in ("pattern", "pitch_um", "kapton_thickness_um",
                                    "copper_thickness_um", "kapton_permittivity", "hole")}
            ff["name"] = f["name"]
            ff["offset_mm"] = [round(v, 9) for v in self.offsets[f["name"]]]
            foils.append(ff)
        return {"version": FIELD_MODEL_VERSION, "foils": foils,
                "z": {k: round(v, 9) for k, v in self.z.items()},
                "V": {k: round(v, 6) for k, v in self.V.items()},
                "cell": [round(c, 9) for c in self.cell], "mesh": self.cfg["mesh"]}

    def field_hash(self):
        blob = json.dumps(self.field_key(), sort_keys=True).encode()
        return hashlib.sha1(blob).hexdigest()[:12]


def gas_key(cfg):
    g = cfg["gas"]
    return {"components": {k.lower(): v for k, v in sorted(g["components"].items())},
            "pressure_bar": g["pressure_bar"], "temperature_K": g["temperature_K"]}


BAR_PER_PSI = 0.0689475729

GAS_DISPLAY = {"ar": "Ar", "ch4": "CH4", "co2": "CO2", "cf4": "CF4", "he": "He", "ne": "Ne",
               "xe": "Xe", "kr": "Kr", "n2": "N2", "h2": "H2", "o2": "O2", "h2o": "H2O",
               "ic4h10": "iC4H10", "c2h6": "C2H6", "c3h8": "C3H8", "c2h2": "C2H2",
               "dme": "DME", "tma": "TMA"}


def gas_label(cfg):
    """e.g. 'Ar:CH4 92:8, 1.01325 bar, 293.15 K' (mixtures written with ':')."""
    g = cfg["gas"]
    tot = sum(g["components"].values())
    comp = ":".join(GAS_DISPLAY.get(k.lower(), k) for k in g["components"])
    frac = ":".join(f"{100 * v / tot:g}" for v in g["components"].values())
    return f"{comp} {frac}, {g['pressure_bar']:g} bar, {g['temperature_K']:g} K"


def run_name(config_path):
    return Path(config_path).stem
