"""
Empirical gain corrections (the optional "corrections" block of a config).

The simulated avalanche has no dielectric charging-up, and even with Penning
matching Garfield++ stays below measured GEM gains by a factor per foil whose
cause is not known (CERN GDD, RD51 WG4). Both act per multiplication stage, so
a stack of n foils is corrected by

    G_corrected = G_simulated * (charging_up_per_foil * gain_scale_per_foil) ** n.

The raw simulated gains are never changed; scans and summaries report the
corrected values next to them. Corrections do not enter the run key: changing
them never re-runs a simulation.

    "corrections": "literature"          -> LITERATURE below
    "corrections": {"charging_up_per_foil": 1.2, "gain_scale_per_foil": 1.3}
"""
import copy

KEYS = {"charging_up_per_foil", "gain_scale_per_foil", "requires_penning_match", "notes"}

# Informed choices from published data (2026-10-27); see
# docs/validation.md, section 4.
LITERATURE = {
    # Surface charging-up of pre-polarised double-mask CERN foils, measured:
    # +3..+9 % per foil (Hauer et al. NIM A 976 (2020) 164205; Correia et al.
    # JINST 9 (2014) P07025, data; Chatterjee et al. arXiv:2007.11444, 2107.00890,
    # 2206.10876, ~x1.12 per pre-polarised triple). Simulations of it overshoot
    # 2-3x. Kapton polarisation after switch-on is excluded: the field map uses
    # the static permittivity, i.e. the steady polarised state.
    "charging_up_per_foil": 1.05,
    # Measured / simulated effective gain per foil left after Penning matching
    # and the charging-up above, from current-mode single-GEM data in Ar:CO2
    # 70:30 (no pulser calibration): CERN GDD (Janssens 2019) 2.09, Bachmann et
    # al. NIM A 438 (1999) Fig. 24 ~2.2, Fig. 13 1.37; log-mean 1.85, scatter
    # x/÷ 1.3. Consistent with the x4.5-6.8 BESIII CGEM tuning of Garfield++
    # triple-GEM gains (arXiv:1904.06142; RD51 WG4), i.e. 1.65-1.9 per foil.
    "gain_scale_per_foil": 1.85,
    "requires_penning_match": True,
    "notes": ["charging_up_per_foil 1.05 (1.02-1.10): surface charging, pre-polarised foils",
              "gain_scale_per_foil 1.85 x/÷ 1.3: current-mode single-GEM data vs SimGEM "
              "with Penning matching (Ar:CO2 70:30)"],
}


def resolve(spec):
    """Expanded corrections block (dict) or None."""
    if spec is None:
        return None
    if spec == "literature":
        return copy.deepcopy(LITERATURE)
    if not isinstance(spec, dict):
        raise ValueError('corrections must be null, "literature" or a mapping')
    unknown = set(spec) - KEYS
    if unknown:
        raise ValueError(f"corrections: unknown key(s) {sorted(unknown)}")
    out = {"charging_up_per_foil": 1.0, "gain_scale_per_foil": 1.0,
           "requires_penning_match": False, "notes": []}
    out.update(copy.deepcopy(spec))
    for k in ("charging_up_per_foil", "gain_scale_per_foil"):
        if not (isinstance(out[k], (int, float)) and out[k] > 0):
            raise ValueError(f"corrections.{k} must be a positive number")
    return out


def factor(cfg):
    """Multiplicative correction of the stack gain (1.0 without corrections)."""
    c = cfg.get("corrections")
    if not c:
        return 1.0
    return (c["charging_up_per_foil"] * c["gain_scale_per_foil"]) ** len(cfg["foils"])
