#!/usr/bin/env python3
"""
Gas-mixture study: the list of mixtures, as a workbook and as the CSV the sweep
driver reads.

    python studies/gas_mixture_study/build_mixtures.py   # -> gas_mixtures.xlsx, mixtures.csv

Edit the lists below (or the input cells of the workbook) and rerun. The plan
is in PLAN.md next to this file.
"""
import csv
import os

from openpyxl import Workbook
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------- components
# key: (Garfield name, formula, molar mass g/mol, ionisation energy eV,
#       LEL in air vol% (None = not flammable), nuclei, Garfield automatic Penning with Ar,
#       role, hazards / remarks)
COMPONENTS = {
    "ar": ("Ar", "Ar", 39.948, 15.76, None, "Ar", "-", "target, main gas", ""),
    "he": ("He", "He", 4.0026, 24.59, None, "He", "none", "diluent", "fast leaks through seals"),
    "ne": ("Ne", "Ne", 20.180, 21.56, None, "Ne", "none (Ne-CO2, Ne-N2 exist)", "diluent", "cost at 10 bar volumes"),
    "ch4": ("CH4", "CH4", 16.043, 12.61, 4.4, "H, C", "p- and c-dependent (Sahin 2010)",
            "quencher, Penning acceptor (Ar 4p and up)", "flammable"),
    "co2": ("CO2", "CO2", 44.010, 13.78, None, "C, O", "1 atm: Sahin 2014; other p: JINST 12 C01035",
            "quencher (cold gas, low diffusion)", "three-body attachment with O2"),
    "cf4": ("CF4", "CF4", 88.004, 16.2, None, "C, F", "none",
            "fast drift, low diffusion", "attachment at hole fields; GWP ~6500; IP: vertical 16.2 eV"),
    "c2h6": ("C2H6", "C2H6", 30.069, 11.52, 2.4, "H, C", "0.31 (10 %, 1 atm only)",
             "quencher, Penning donor (marginal for Ar 4s)", "flammable"),
    "c3h8": ("C3H8", "C3H8", 44.097, 10.94, 1.7, "H, C", "c-dependent (1 atm only)",
             "quencher, Penning donor", "flammable"),
    "ic4h10": ("iC4H10", "iC4H10", 58.122, 10.68, 1.3, "H, C", "0.40 (10 %, 1 atm only)",
               "quencher, Penning donor, UV absorber", "flammable; condenses at ~2-3 bar partial pressure at room T"),
    "dme": ("DME", "C2H6O", 46.068, 10.03, 2.7, "H, C, O", "none",
            "quencher, Penning donor, no attachment", "flammable; attacks some plastics"),
    "tma": ("TMA", "C3H9N", 59.110, 7.85, 2.0, "H, C, N", "none (Xe-TMA only)",
            "strong Penning donor, UV absorber", "flammable, toxic, reactive"),
    "c2h2": ("C2H2", "C2H2", 26.038, 11.40, 2.3, "H, C", "0.72 (cylindrical counters)",
             "Penning donor", "decomposes explosively above ~1-2 bar partial pressure"),
    "n2": ("N2", "N2", 28.014, 15.58, None, "N", "none", "reference (no Penning, poor quencher)", ""),
    "h2": ("H2", "H2", 2.016, 15.43, 4.0, "H", "none", "free-proton target, poor quencher", "flammable"),
    "o2": ("O2", "O2", 31.999, 12.07, None, "O", "-", "impurity (attachment)", ""),
    "h2o": ("H2O", "H2O", 18.015, 12.62, None, "H, O", "-", "impurity (attachment, cooling)", ""),
}

# ----------------------------------------------------------------- mixtures
# (ID, family, priority, {component: vol %} without Ar, rationale, remarks, validation data)
# Priority A: all pressures + Penning matching at 1 and 10 bar; B: 1, 4, 10 bar; C: 1, 10 bar.
MIXTURES = [
    # references
    ("R01", "Reference", "C", {}, "Pure Ar: baseline for every alpha and drift comparison", "Not operable (no quencher)", "Magboltz validation data"),
    ("R02", "Reference", "C", {"ne": 10}, "Ne effect on alpha with no quencher", "Not operable", ""),
    ("R03", "Reference", "C", {"he": 10}, "He effect on alpha with no quencher", "Not operable", ""),
    # Ar:CH4
    ("B01", "Ar:CH4", "A", {"ch4": 2}, "Lowest-quencher candidate; SimGEM's original 8 bar mixture", "Photon feedback at high gain", "Ye et al. 98:2 counter, 1-28 atm"),
    ("B02", "Ar:CH4", "B", {"ch4": 3}, "Fills the CH4 fraction scan", "", ""),
    ("B03", "Ar:CH4", "A", {"ch4": 5}, "Middle of the CH4 scan; reference for the Ne/He and donor ternaries", "", ""),
    ("B04", "Ar:CH4", "B", {"ch4": 7}, "Fills the CH4 fraction scan", "", ""),
    ("B05", "Ar:CH4", "A", {"ch4": 10}, "ND-GAr baseline (P10)", "", "Ye et al. 90:10 counter, 1-21 atm; P10 drift data"),
    # Ar:CO2
    ("B06", "Ar:CO2", "B", {"co2": 2}, "Low-diffusion family, lowest fraction", "Photon feedback", ""),
    ("B07", "Ar:CO2", "C", {"co2": 3}, "Fills the CO2 fraction scan", "", ""),
    ("B08", "Ar:CO2", "A", {"co2": 5}, "Low-diffusion candidate", "O2 three-body attachment", ""),
    ("B09", "Ar:CO2", "C", {"co2": 7}, "Fills the CO2 fraction scan", "", ""),
    ("B10", "Ar:CO2", "A", {"co2": 10}, "Low-diffusion candidate; GORG data at 1 bar", "O2 three-body attachment", "GORG Ar:CO2 90:10 (15 psi)"),
    # Ar:CF4
    ("B11", "Ar:CF4", "C", {"cf4": 2}, "Fast drift, low diffusion", "Hole-field attachment", ""),
    ("B12", "Ar:CF4", "B", {"cf4": 5}, "Fast drift, low diffusion", "Hole-field attachment; adds F nuclei", ""),
    ("B13", "Ar:CF4", "B", {"cf4": 10}, "Fast drift, low diffusion", "Hole-field attachment; adds F nuclei", ""),
    # Ar:C2H6
    ("B14", "Ar:C2H6", "C", {"c2h6": 2}, "Better UV quencher than CH4", "Flammable", ""),
    ("B15", "Ar:C2H6", "B", {"c2h6": 5}, "Better UV quencher than CH4; Penning from Ar 4s marginal", "Flammable", ""),
    ("B16", "Ar:C2H6", "C", {"c2h6": 10}, "Only point with a published r (0.31, 1 atm)", "Flammable", "Sahin 2010 (r at 1 atm)"),
    # Ar:C3H8
    ("B17", "Ar:C3H8", "C", {"c3h8": 2}, "Penning donor (IP 10.9 eV)", "Flammable", "Sahin 2010 (r at 1 atm)"),
    ("B18", "Ar:C3H8", "C", {"c3h8": 5}, "Penning donor (IP 10.9 eV)", "Flammable", ""),
    # Ar:iC4H10
    ("B19", "Ar:iC4H10", "B", {"ic4h10": 1}, "Donor alone, lowest fraction", "Flammable", ""),
    ("B20", "Ar:iC4H10", "A", {"ic4h10": 2}, "Donor alone; same isobutane as T2K gas", "Flammable", ""),
    ("B21", "Ar:iC4H10", "B", {"ic4h10": 3}, "Donor alone", "Flammable", ""),
    ("B22", "Ar:iC4H10", "C", {"ic4h10": 5}, "Donor alone", "Flammable; partial pressure 0.5 bar at 10 bar", ""),
    # Ar:DME
    ("B23", "Ar:DME", "B", {"dme": 1}, "Donor, no attachment", "Flammable; materials compatibility", ""),
    ("B24", "Ar:DME", "B", {"dme": 2}, "Donor, no attachment", "Flammable; materials compatibility", ""),
    ("B25", "Ar:DME", "C", {"dme": 5}, "Donor, no attachment", "Flammable; materials compatibility", ""),
    # Ar:TMA
    ("B26", "Ar:TMA", "B", {"tma": 0.5}, "Strongest Penning donor (IP 7.9 eV)", "Toxic, reactive", ""),
    ("B27", "Ar:TMA", "B", {"tma": 1}, "Strongest Penning donor", "Toxic, reactive", ""),
    ("B28", "Ar:TMA", "C", {"tma": 2}, "Strongest Penning donor", "Toxic, reactive", ""),
    # Ar:C2H2
    ("B29", "Ar:C2H2", "C", {"c2h2": 1}, "Published high r (0.72)", "Partial pressure 0.1 bar at 10 bar; decomposition hazard", "Sahin 2010 (r)"),
    ("B30", "Ar:C2H2", "C", {"c2h2": 2}, "Published high r (0.72)", "Partial pressure 0.2 bar at 10 bar; decomposition hazard", ""),
    # Ar:N2, Ar:H2
    ("B31", "Ar:N2", "C", {"n2": 5}, "Reference: molecular gas without Penning", "Poor quencher", ""),
    ("B32", "Ar:N2", "C", {"n2": 10}, "Reference: molecular gas without Penning", "Poor quencher", ""),
    ("B33", "Ar:H2", "C", {"h2": 2}, "Free-proton target; no Penning", "Poor quencher; flammable", ""),
    ("B34", "Ar:H2", "C", {"h2": 5}, "Free-proton target; no Penning", "Poor quencher; flammable", ""),
    # Ne / He dilution
    ("N01", "Ar:Ne:CH4", "A", {"ne": 5, "ch4": 5}, "Does Ne raise alpha? Compare B03 (95:5) and B05 (90:10)", "", ""),
    ("N02", "Ar:He:CH4", "A", {"he": 5, "ch4": 5}, "Does He raise alpha? Compare B03 and B05", "", ""),
    ("N03", "Ar:Ne:CH4", "B", {"ne": 5, "ch4": 3}, "Ne with less quencher", "", ""),
    ("N04", "Ar:He:CH4", "C", {"he": 5, "ch4": 3}, "He with less quencher", "", ""),
    ("N05", "Ar:Ne:CO2", "B", {"ne": 5, "co2": 5}, "Ne in the low-diffusion family", "", ""),
    ("N06", "Ar:He:CO2", "B", {"he": 5, "co2": 5}, "He in the low-diffusion family", "", ""),
    ("N07", "Ar:Ne:CO2", "C", {"ne": 5, "co2": 3}, "Ne with less CO2", "", ""),
    ("N08", "Ar:Ne:iC4H10", "C", {"ne": 5, "ic4h10": 2}, "Ne with a Penning donor", "Flammable", ""),
    ("N09", "Ar:Ne:CH4", "C", {"ne": 7, "ch4": 3}, "Maximum Ne within 90 % Ar", "", ""),
    # quencher + Penning donor
    ("P01", "Ar:CH4:iC4H10", "A", {"ch4": 5, "ic4h10": 1}, "CH4 drift + donor for gain", "Flammable", ""),
    ("P02", "Ar:CH4:iC4H10", "A", {"ch4": 5, "ic4h10": 2}, "CH4 drift + donor for gain", "Flammable", ""),
    ("P03", "Ar:CH4:iC4H10", "B", {"ch4": 7, "ic4h10": 2}, "Donor at higher CH4", "Flammable", ""),
    ("P04", "Ar:CH4:iC4H10", "B", {"ch4": 8, "ic4h10": 2}, "P10-like with donor", "Flammable", ""),
    ("P05", "Ar:CH4:C2H6", "B", {"ch4": 5, "c2h6": 3}, "Ethane as a milder donor", "Flammable", ""),
    ("P06", "Ar:CH4:C2H6", "C", {"ch4": 5, "c2h6": 5}, "Ethane as a milder donor", "Flammable", ""),
    ("P07", "Ar:CH4:DME", "B", {"ch4": 5, "dme": 1}, "DME donor, no attachment", "Flammable; materials", ""),
    ("P08", "Ar:CH4:DME", "C", {"ch4": 5, "dme": 2}, "DME donor", "Flammable; materials", ""),
    ("P09", "Ar:CH4:TMA", "B", {"ch4": 5, "tma": 0.5}, "Trace TMA as donor", "Toxic, reactive", ""),
    ("P10", "Ar:CH4:TMA", "B", {"ch4": 5, "tma": 1}, "TMA as donor", "Toxic, reactive", ""),
    ("P11", "Ar:CH4:C2H2", "C", {"ch4": 5, "c2h2": 1}, "Acetylene as donor", "Decomposition hazard", ""),
    ("P12", "Ar:CH4:C3H8", "C", {"ch4": 5, "c3h8": 2}, "Propane as donor", "Flammable", ""),
    ("P13", "Ar:CO2:iC4H10", "B", {"co2": 5, "ic4h10": 1}, "Low diffusion + donor", "Flammable; O2 attachment", ""),
    ("P14", "Ar:CO2:iC4H10", "A", {"co2": 5, "ic4h10": 2}, "Low diffusion + donor", "Flammable; O2 attachment", ""),
    ("P15", "Ar:CO2:DME", "B", {"co2": 5, "dme": 1}, "Low diffusion + DME donor", "Flammable; materials", ""),
    ("P16", "Ar:CO2:DME", "C", {"co2": 5, "dme": 2}, "Low diffusion + DME donor", "Flammable; materials", ""),
    ("P17", "Ar:CO2:TMA", "C", {"co2": 5, "tma": 0.5}, "Low diffusion + trace TMA", "Toxic, reactive", ""),
    ("P18", "Ar:CO2:C2H6", "C", {"co2": 5, "c2h6": 3}, "Low diffusion + ethane", "Flammable", ""),
    # CF4-based
    ("F01", "Ar:CF4:iC4H10", "A", {"cf4": 3, "ic4h10": 2}, "T2K gas: proven Ar-rich TPC mixture", "CF4 attachment at 10 bar hole fields", "T2K TPC drift velocity, diffusion (1 bar)"),
    ("F02", "Ar:CF4:iC4H10", "B", {"cf4": 3, "ic4h10": 3}, "T2K gas with more donor", "Flammable", ""),
    ("F03", "Ar:CF4:iC4H10", "B", {"cf4": 5, "ic4h10": 2}, "T2K gas with more CF4", "More attachment", ""),
    ("F04", "Ar:CF4:CH4", "B", {"cf4": 3, "ch4": 5}, "CH4 instead of isobutane", "", ""),
    ("F05", "Ar:CF4:CH4", "C", {"cf4": 5, "ch4": 5}, "CH4 instead of isobutane", "", ""),
    ("F06", "Ar:CF4:C2H6", "C", {"cf4": 3, "c2h6": 2}, "Ethane instead of isobutane", "Flammable", ""),
    ("F07", "Ar:CF4:DME", "C", {"cf4": 3, "dme": 2}, "DME instead of isobutane", "Flammable; materials", ""),
    ("F08", "Ar:CO2:CF4", "C", {"co2": 5, "cf4": 3}, "Non-flammable fast mixture", "Attachment", ""),
    # quaternaries
    ("Q01", "Ar:Ne:CH4:iC4H10", "C", {"ne": 3, "ch4": 5, "ic4h10": 2}, "Ne + donor combined", "Flammable", ""),
    ("Q02", "Ar:CF4:CH4:iC4H10", "C", {"cf4": 3, "ch4": 5, "ic4h10": 2}, "T2K gas with CH4", "Flammable", ""),
    ("Q03", "Ar:He:CF4:iC4H10", "C", {"he": 5, "cf4": 3, "ic4h10": 2}, "T2K gas with He", "Flammable", ""),
]

# impurity runs (electron lifetime): base mixtures x impurity levels, drift fields only
IMPURITY_BASES = ["B03", "B05", "B08", "B10", "P02", "P14", "F01"]
IMPURITIES = [("o2", 1), ("o2", 10), ("o2", 100), ("h2o", 10), ("h2o", 100)]   # ppm

# ----------------------------------------------------------------- styles
FONT = "Arial"
F_BODY = Font(name=FONT, size=10)
F_HEAD = Font(name=FONT, size=10, bold=True, color="FFFFFF")
F_INPUT = Font(name=FONT, size=10, color="0000FF")
F_TITLE = Font(name=FONT, size=14, bold=True)
F_BOLD = Font(name=FONT, size=10, bold=True)
FILL_HEAD = PatternFill("solid", fgColor="305496")
FILL_A = PatternFill("solid", fgColor="E2EFDA")
FILL_B = PatternFill("solid", fgColor="FFF2CC")
FILL_C = PatternFill("solid", fgColor="F2F2F2")
FILL_KEY = PatternFill("solid", fgColor="FFFF00")
THIN = Side(style="thin", color="BFBFBF")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
WRAP = Alignment(wrap_text=True, vertical="top")


def header(ws, row, names, widths=None):
    for j, n in enumerate(names, 1):
        c = ws.cell(row=row, column=j, value=n)
        c.font, c.fill, c.alignment, c.border = F_HEAD, FILL_HEAD, Alignment(wrap_text=True, vertical="center"), BOX
    if widths:
        for j, w in enumerate(widths, 1):
            ws.column_dimensions[get_column_letter(j)].width = w


def body(c, font=F_BODY, fmt=None):
    c.font, c.border, c.alignment = font, BOX, WRAP
    if fmt:
        c.number_format = fmt


def label(mix):
    names = ["Ar"] + [COMPONENTS[k][0] for k in mix]
    ar = 100 - sum(mix.values())
    vals = [ar] + list(mix.values())
    return ":".join(names) + " " + ":".join(f"{v:g}" for v in vals)


def nuclei(mix):
    out = ["Ar"]
    for k in mix:
        for n in COMPONENTS[k][5].split(", "):
            if n not in out:
                out.append(n)
    return ", ".join(out[1:]) or "-"


def penning_note(mix):
    keys = list(mix)
    if not keys:
        return "n/a"
    if len(keys) == 1:
        return "auto: " + COMPONENTS[keys[0]][6] if COMPONENTS[keys[0]][6] not in ("none", "-") \
            and not COMPONENTS[keys[0]][6].startswith("none") else "none: bracket r"
    return "none for ternaries: bracket r"


# ----------------------------------------------------------------- workbook
def build():
    wb = Workbook()

    # --- Read me
    ws = wb.active
    ws.title = "Read me"
    ws["A1"] = "Gas-mixture study: mixtures for an Ar-rich (>= 90 %) TPC at up to 10 bar"
    ws["A1"].font = F_TITLE
    rows = [
        ("Purpose", "Magboltz transport and Townsend properties for every candidate mixture at 1-10 bar, "
                    "to choose the 4-6 mixtures for the single-GEM simulations. Plan: PLAN.md."),
        ("Constraints", "Ar >= 90 % by volume; only nuclei lighter than Ar; drift (diffusion, velocity, "
                        "attachment, electron lifetime) acceptable for a long-drift TPC; gain reachable with GEMs at 10 bar."),
        ("Mixtures", "One row per mixture. Volume fractions in %, Ar = 100 - the rest. Priority sets the pressures "
                     "run (Run grid). Grey columns are formulas."),
        ("Components", "Properties used by the formulas (molar mass, ionisation energy, LEL), with sources."),
        ("Criteria", "Detector requirements (drift length, drift field, gain) and the per-foil gain targets derived from them."),
        ("Run grid", "Pressures, field grids, Penning r grid and the compute budget. Blue cells are inputs."),
        ("Impurities", "Electron-lifetime runs: O2 / H2O admixtures to selected base mixtures (drift fields only)."),
        ("Outputs", "What is extracted per mixture and pressure, and the figure of merit it feeds."),
        ("Validation", "Measured data to check Magboltz + Garfield against before trusting the sweep."),
        ("Legend", "Blue text: input, edit freely. Black: formula. Yellow fill: key assumption to revisit. "
                   "Priority fill: green A (first, all pressures), yellow B (1, 4, 10 bar), grey C (1, 10 bar)."),
        ("Regenerate", "python studies/gas_mixture_study/build_mixtures.py (also writes mixtures.csv for the driver)."),
    ]
    for i, (k, v) in enumerate(rows, 3):
        ws.cell(row=i, column=1, value=k).font = F_BOLD
        c = ws.cell(row=i, column=2, value=v)
        c.font, c.alignment = F_BODY, WRAP
    ws.column_dimensions["A"].width = 16
    ws.column_dimensions["B"].width = 110

    # --- Components
    wc = wb.create_sheet("Components")
    cols = ["Key", "Garfield name", "Formula", "Molar mass (g/mol)", "Ionisation energy (eV)",
            "LEL in air (vol %)", "Nuclei", "Garfield automatic Penning with Ar", "Role", "Hazards / remarks"]
    header(wc, 1, cols, [9, 12, 9, 12, 13, 12, 10, 34, 34, 44])
    for i, (k, v) in enumerate(COMPONENTS.items(), 2):
        vals = [k, v[0], v[1], v[2], v[3], v[4] if v[4] is not None else 0, v[5], v[6], v[7], v[8]]
        for j, x in enumerate(vals, 1):
            body(wc.cell(row=i, column=j, value=x), F_INPUT if j in (4, 5, 6) else F_BODY,
                 "0.000" if j == 4 else "0.00" if j == 5 else "0.0" if j == 6 else None)
    n_comp = len(COMPONENTS) + 1
    src = n_comp + 2
    notes = [
        "Sources: molar masses IUPAC; ionisation energies NIST Chemistry WebBook (adiabatic; CF4 vertical value).",
        "LEL: typical values of IEC 60079-20-1 (0 = not flammable); approximate, check the supplier SDS. "
        "Mixture flammability needs an ISO 10156 assessment.",
        "Ar excited levels: 4s 11.55-11.83 eV (most abundant), 4p 12.91-13.48 eV, 3d/5s ~13.85-14.3 eV. "
        "Garfield transfers only from levels above the lowest ionisation energy of the mixture.",
        "Garfield automatic Penning: MediumGas::EnablePenningTransfer() (Garfield++ 58d22389), Sahin et al. "
        "JINST 5 (2010) P05002, NIM A 768 (2014) 104, JINST 12 (2017) C01035.",
    ]
    for i, t in enumerate(notes):
        c = wc.cell(row=src + i, column=1, value=t)
        c.font = F_BODY
    wc.freeze_panes = "B2"
    rng = f"Components!$A$2:$A${n_comp}"

    def look(key_cell, col):
        return f'INDEX(Components!${col}$2:${col}${n_comp},MATCH({key_cell},{rng},0))'

    # --- Run grid (inputs referenced by Mixtures)
    wr = wb.create_sheet("Run grid")
    wr["A1"] = "Run grid and compute budget"
    wr["A1"].font = F_TITLE
    grid = [
        ("Pressures, priority A (bar)", "1, 2, 4, 6, 8, 10", None),
        ("Number of pressures, priority A", 6, None),
        ("Pressures, priority B (bar)", "1, 4, 10", None),
        ("Number of pressures, priority B", 3, None),
        ("Pressures, priority C (bar)", "1, 10", None),
        ("Number of pressures, priority C", 2, None),
        ("Drift table: E/p range (V/cm/bar)", "5 - 500", "At 10 bar: 50 V/cm - 5 kV/cm"),
        ("Drift table: field points (log)", 12, None),
        ("Avalanche table: E/p range (kV/cm/bar)", "5 - 200", "GEM holes: ~100 kV/cm at 1 bar, ~300-600 kV/cm at 10 bar"),
        ("Avalanche table: field points (log)", 15, None),
        ("Magboltz collisions (x 1e7)", 5, None),
        ("CPU s per drift point", 60, "Measured 2026-10-28: Ar:CH4 90:10 and Ar:CO2 95:5 at 10 bar, 1 core"),
        ("CPU s per avalanche point", 120, "Measured 2026-10-28: Ar:CH4 90:10 at 10 bar, 100 kV/cm - 1 MV/cm"),
        ("Penning r grid (post-hoc, from one table)", "0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, auto", "Free: rescales alpha from stored excitation rates"),
        ("Penning matching (r_micro): CPU min per run", 15, "simgem_penning: 4-25 min depending on gas"),
        ("Penning matching: pressures per priority-A mixture", 2, "1 and 10 bar"),
        ("Cores", 4, "CPU cap: at most 6 heavy processes"),
    ]
    for i, (k, v, note) in enumerate(grid, 3):
        wr.cell(row=i, column=1, value=k).font = F_BODY
        c = wr.cell(row=i, column=2, value=v)
        c.font = F_INPUT
        if note:
            wr.cell(row=i, column=3, value=note).font = F_BODY
    for r in (14, 15):                                   # measured timings: key assumptions
        wr.cell(row=r, column=2).fill = FILL_KEY
    R = {k: f"'Run grid'!$B${i}" for i, (k, _, _) in enumerate(grid, 3)}
    nP = {"A": R["Number of pressures, priority A"], "B": R["Number of pressures, priority B"],
          "C": R["Number of pressures, priority C"]}
    per_table = (f"({R['Drift table: field points (log)']}*{R['CPU s per drift point']}"
                 f"+{R['Avalanche table: field points (log)']}*{R['CPU s per avalanche point']})/3600")
    r0 = 3 + len(grid) + 1
    wr.cell(row=r0, column=1, value="Budget").font = F_BOLD
    budget = [
        ("CPU h per mixture and pressure", f"={per_table}", "0.00"),
        ("Priority A mixtures", '=COUNTIF(Mixtures!$C:$C,"A")', "0"),
        ("Priority B mixtures", '=COUNTIF(Mixtures!$C:$C,"B")', "0"),
        ("Priority C mixtures", '=COUNTIF(Mixtures!$C:$C,"C")', "0"),
        ("Tables (mixture x pressure)", f"=B{r0+2}*{nP['A']}+B{r0+3}*{nP['B']}+B{r0+4}*{nP['C']}", "0"),
        ("Tables, impurity runs (3 pressures, drift only)", '=3*COUNTIF(Impurities!$A$2:$A$500,"I*")', "0"),
        ("CPU h, mixture tables", f"=B{r0+5}*B{r0+1}", "0.0"),
        ("CPU h, impurity tables",
         f"=B{r0+6}*{R['Drift table: field points (log)']}*{R['CPU s per drift point']}/3600", "0.0"),
        ("CPU h, Penning matching",
         f"=B{r0+2}*{R['Penning matching: pressures per priority-A mixture']}*{R['Penning matching (r_micro): CPU min per run']}/60", "0.0"),
        ("CPU h, total", f"=B{r0+7}+B{r0+8}+B{r0+9}", "0.0"),
        ("Wall time, total (days)", f"=B{r0+10}/{R['Cores']}/24", "0.0"),
        ("Wall time, priority A only (h)",
         f"=(B{r0+2}*{nP['A']}*B{r0+1}+B{r0+9})/{R['Cores']}", "0.0"),
    ]
    for i, (k, f, fmt) in enumerate(budget, r0 + 1):
        wr.cell(row=i, column=1, value=k).font = F_BODY
        c = wr.cell(row=i, column=2, value=f)
        c.font, c.number_format = F_BODY, fmt
    wr.column_dimensions["A"].width = 48
    wr.column_dimensions["B"].width = 34
    wr.column_dimensions["C"].width = 70

    # --- Mixtures
    wm = wb.create_sheet("Mixtures", 1)
    cols = ["ID", "Family", "Priority", "Mixture (vol %)",
            "Comp 2", "x2 (%)", "Comp 3", "x3 (%)", "Comp 4", "x4 (%)",
            "Ar (vol %)", "Mean molar mass (g/mol)", "Ar mass fraction", "Lowest IP (eV)",
            "Ar levels eligible for Penning", "Flammable content (vol %)", "Non-Ar nuclei",
            "Penning r in Garfield", "Pressures run", "CPU h (tables)", "Rationale", "Risks / remarks",
            "Validation data"]
    header(wm, 1, cols, [6, 17, 8, 26, 8, 7, 8, 7, 8, 7, 8, 10, 9, 9, 22, 10, 11, 30, 9, 8, 46, 38, 30])
    pri_fill = {"A": FILL_A, "B": FILL_B, "C": FILL_C}
    for i, (mid, fam, pri, mix, why, risk, val) in enumerate(MIXTURES, 2):
        comps = list(mix.items()) + [("", None)] * (3 - len(mix))
        fixed = [mid, fam, pri, label(mix)]
        for j, x in enumerate(fixed, 1):
            body(wm.cell(row=i, column=j, value=x), F_INPUT if j == 3 else F_BODY)
        wm.cell(row=i, column=3).fill = pri_fill[pri]
        for n, (k, x) in enumerate(comps):
            body(wm.cell(row=i, column=5 + 2 * n, value=k or None), F_INPUT)
            body(wm.cell(row=i, column=6 + 2 * n, value=x), F_INPUT, "0.0")
        E, Fc, G, H, I, J, K = (f"{c}{i}" for c in "EFGHIJK")
        m_ar = look('"ar"', "D")
        mm = "+".join(f'IF({k}="",0,{x}*{look(k, "D")})' for k, x in ((E, Fc), (G, H), (I, J)))
        ip = ",".join(f'IF({k}="",99,{look(k, "E")})' for k in (E, G, I))
        fl = "+".join(f'IF({k}="",0,IF({look(k, "F")}>0,{x},0))' for k, x in ((E, Fc), (G, H), (I, J)))
        formulas = [
            (11, f"=100-{Fc}-{H}-{J}", "0.0"),
            (12, f"=({K}*{m_ar}+{mm})/100", "0.00"),
            (13, f"={K}*{m_ar}/({K}*{m_ar}+{mm})", "0.0%"),
            (14, f"=MIN(15.76,{ip})", "0.00"),
            (15, f'=IF(N{i}<11.55,"4s, 4p and higher",IF(N{i}<12.91,"4p and higher",'
                 f'IF(N{i}<13.85,"3d/5s and higher",IF(N{i}<15.76,"high Rydberg only","none"))))', None),
            (16, f"={fl}", "0.0"),
            (19, f'=IF(C{i}="A",{nP["A"]},IF(C{i}="B",{nP["B"]},{nP["C"]}))', "0"),
            (20, f"=S{i}*'Run grid'!$B${r0 + 1}", "0.0"),
        ]
        for col, f, fmt in formulas:
            body(wm.cell(row=i, column=col, value=f), F_BODY, fmt)
        body(wm.cell(row=i, column=17, value=nuclei(mix)))
        body(wm.cell(row=i, column=18, value=penning_note(mix)))
        for col, x in ((21, why), (22, risk), (23, val)):
            body(wm.cell(row=i, column=col, value=x or None))
    last = len(MIXTURES) + 1
    wm.freeze_panes = "E2"
    wm.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{last}"
    dv = DataValidation(type="list", formula1='"A,B,C"', allow_blank=False)
    wm.add_data_validation(dv)
    dv.add(f"C2:C{last}")
    wm.cell(row=1, column=13).comment = Comment("Mass fraction of Ar in the gas: the fraction of target "
                                                "nuclei mass that is argon.", "SimGEM")
    wm.cell(row=1, column=15).comment = Comment("Ar excitations that can Penning-ionise the lowest-IP "
                                                "component (Garfield applies r only to these).", "SimGEM")
    wm.cell(row=1, column=16).comment = Comment("Sum of the flammable components; compare with their LEL "
                                                "and do an ISO 10156 assessment for the mixture.", "SimGEM")
    t = last + 2
    for k, v in (("Count", f'=COUNTIF(A2:A{last},"?*")'), ("Priority A", f'=COUNTIF(C2:C{last},"A")'),
                 ("CPU h, all tables", f"=SUM(T2:T{last})")):
        wm.cell(row=t, column=1, value=k).font = F_BOLD
        c = wm.cell(row=t, column=4, value=v)
        c.font, c.number_format = F_BODY, "0.0"
        t += 1

    # --- Impurities
    wi = wb.create_sheet("Impurities")
    header(wi, 1, ["Run", "Base ID", "Base mixture", "Impurity", "Level (ppm)", "Pressures (bar)",
                   "Fields", "Purpose"], [9, 8, 26, 9, 10, 14, 18, 60])
    n = 2
    lookup = {m[0]: m for m in MIXTURES}
    for b in IMPURITY_BASES:
        for imp, ppm in IMPURITIES:
            vals = [f"I{n - 1:02d}", b, f'=INDEX(Mixtures!$D$2:$D${last},MATCH(B{n},Mixtures!$A$2:$A${last},0))',
                    COMPONENTS[imp][0], ppm, "1, 4, 10", "drift table only",
                    "Attachment rate and electron lifetime vs pressure (O2 three-body: ~p^2); "
                    "check Magboltz third-body treatment against literature rate constants"]
            for j, x in enumerate(vals, 1):
                body(wi.cell(row=n, column=j, value=x), F_INPUT if j in (2, 4, 5) else F_BODY)
            n += 1
        assert b in lookup
    wi.freeze_panes = "A2"

    # --- Criteria
    wq = wb.create_sheet("Criteria", 2)
    wq["A1"] = "Detector requirements (from the user, 2026-10-28) and derived targets"
    wq["A1"].font = F_TITLE
    crit = [
        ("Maximum drift length, double-sided TPC (cm)", 250, None, "0"),
        ("Maximum drift length, single-sided TPC (cm)", 500, None, "0"),
        ("Nominal drift field (kV/cm)", 0.5, None, "0.00"),
        ("Operating pressure (bar)", 10, None, "0"),
        ("Nominal drift E/p (V/cm/bar)", "=B5*1000/B6", "Lower pressures: evaluate at the same E/p and at 0.5 kV/cm", "0"),
        ("Target effective gain of the stack", 100, "'At least O(100)'; revisit", "0"),
        ("Effective/absolute gain per foil (assumed)", 0.4, "SimGEM single standard GEM: 0.3-0.45", "0.00"),
        ("Per-foil effective gain, triple GEM", "=B8^(1/3)", None, "0.0"),
        ("Per-foil effective gain, double GEM", "=B8^(1/2)", None, "0.0"),
        ("Per-foil effective gain, single GEM", "=B8", None, "0"),
        ("Per-foil absolute gain, triple GEM", "=B10/B9", "Target of the 1-D gain estimate", "0"),
        ("Per-foil absolute gain, double GEM", "=B11/B9", "Target of the 1-D gain estimate", "0"),
        ("Per-foil absolute gain, single GEM", "=B12/B9", "Target of the 1-D gain estimate", "0"),
        ("Diffusion scale sqrt(L), double-sided (sqrt(cm))", "=SQRT(B3)", "sigma = D x sqrt(L), D in um/sqrt(cm)", "0.0"),
        ("Diffusion scale sqrt(L), single-sided (sqrt(cm))", "=SQRT(B4)", None, "0.0"),
        ("Accepted charge loss over L", 0.1, "Sets the impurity limit (O2, H2O levels unknown)", "0%"),
        ("Additives allowed for the study", "all (flammable, toxic, hazardous)", "Assumed for the study; safety decided later", None),
    ]
    for i, (k, v, note, fmt) in enumerate(crit, 3):
        wq.cell(row=i, column=1, value=k).font = F_BODY
        c = wq.cell(row=i, column=2, value=v)
        c.font = F_BODY if isinstance(v, str) and v.startswith("=") else F_INPUT
        if fmt:
            c.number_format = fmt
        if note:
            wq.cell(row=i, column=3, value=note).font = F_BODY
    wq["B8"].fill = FILL_KEY
    wq["B9"].fill = FILL_KEY
    wq.column_dimensions["A"].width = 50
    wq.column_dimensions["B"].width = 30
    wq.column_dimensions["C"].width = 60

    # --- Outputs
    wo = wb.create_sheet("Outputs")
    header(wo, 1, ["Quantity", "From", "Evaluated at", "Figure of merit / use"], [36, 30, 40, 70])
    outs = [
        ("Drift velocity v_d(E)", "drift table", "E/p 5-500 V/cm/bar", "v_d at the nominal field; flatness dv/dE (stability against field/pressure drifts)"),
        ("Longitudinal diffusion D_L (um/sqrt(cm))", "drift table", "same", "sigma_z over the maximum drift length"),
        ("Transverse diffusion D_T (um/sqrt(cm))", "drift table", "same", "sigma_xy over the maximum drift length vs pad pitch"),
        ("Attachment eta(E)", "drift table", "same", "fraction lost over max drift; lifetime with impurities"),
        ("Characteristic energy eps_k = e D_T / mu", "drift table", "same", "electron temperature; cooling by the quencher"),
        ("Townsend alpha(E), no Penning", "avalanche table", "E/p 5-200 kV/cm/bar", "baseline gain potential"),
        ("Townsend alpha(E; r), r grid + auto", "avalanche table + stored rates", "same", "gain potential vs Penning efficiency"),
        ("Penning-eligible excitations / ionisation", "avalanche table (level rates)", "same", "upper bound of the Penning boost"),
        ("Radiative Ar excitations / ionisation", "avalanche table (level rates)", "same", "photon-feedback proxy: VUV photons per avalanche electron x (1 - r)"),
        ("Attachment in the holes eta(E)", "avalanche table", "same", "gain loss (CF4, O2)"),
        ("E/p needed for alpha/p = target", "derived", "alpha d = ln(per-foil absolute gain, Criteria) over d = 50 um", "ranking of gain reach per mixture and pressure"),
        ("1-D single-GEM gain estimate G(dV)", "derived: alpha integrated along the hole axis of a standard GEM field map", "dV scan", "dV for the per-foil absolute gains of a single, double and triple GEM (Criteria)"),
        ("Drift time over L", "v_d at the nominal field", "L = 250 and 500 cm", "sets the electron-lifetime requirement"),
        ("Max O2 (and H2O) level for 10 % loss over L", "impurity runs: attachment rate vs ppm", "nominal field, 1-10 bar", "purity requirement per mixture (the purity is not known yet)"),
        ("r_micro (Penning matching)", "simgem_penning", "1 and 10 bar, priority A", "input for the single-GEM simulations"),
    ]
    for i, row in enumerate(outs, 2):
        for j, x in enumerate(row, 1):
            body(wo.cell(row=i, column=j, value=x))

    # --- Validation
    wv = wb.create_sheet("Validation")
    header(wv, 1, ["Data", "Mixture / conditions", "Compare", "Status"], [40, 44, 46, 34])
    vals = [
        ("Sahin et al., NIM A 768 (2014) 104", "Ar:CO2, wire counters, ~1 atm", "gain with table + Penning", "done 2026-10-23: within a few %"),
        ("Z. Ye et al. (in Sahin 2010, Fig. 11)", "Ar:CH4 90:10, 1-21 atm; 98:2, 1-28 atm", "counter gain at 7, 14 atm (90:10) and 98:2", "3.5 and 7 atm done: sim/data 0.3-0.7; extend"),
        ("P10 drift velocity / diffusion (literature)", "Ar:CH4 90:10, 1 to several bar", "v_d, D_L, D_T", "to do"),
        ("T2K TPC gas measurements", "Ar:CF4:iC4H10 95:3:2, 1 bar", "v_d, D_T", "to do"),
        ("Three-body O2 attachment rate constants (e.g. Blum, Riegler, Rolandi)", "Ar + O2 with CH4 / CO2 / iC4H10 third bodies", "attachment vs Magboltz", "to do: Magboltz O2 has 3-body attachment; third-body species dependence unclear"),
        ("GORG", "Ar:CH4 92:8 (14.5, 22 psi), Ar:CO2 90:10 (15 psi)", "triple-GEM gain (check of the GEM simulations)", "done with corrections"),
        ("Orthen et al., physics/0302076", "Ar-based triple GEM at several pressures", "gain vs pressure, max gain (discharge limit)", "to do (literature)"),
    ]
    for i, row in enumerate(vals, 2):
        for j, x in enumerate(row, 1):
            body(wv.cell(row=i, column=j, value=x))

    out = os.path.join(HERE, "gas_mixtures.xlsx")
    wb.save(out)
    return out


def write_csv():
    out = os.path.join(HERE, "mixtures.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "family", "priority", "label", "components"])
        for mid, fam, pri, mix, *_ in MIXTURES:
            ar = 100 - sum(mix.values())
            comps = ",".join([f"ar:{ar:g}"] + [f"{k}:{v:g}" for k, v in mix.items()])
            w.writerow([mid, fam, pri, label(mix), comps])
    return out


if __name__ == "__main__":
    print(build())
    print(write_csv())
