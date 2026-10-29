// penning.cpp -- the Penning transfer probability for the MICROSCOPIC avalanche
// that reproduces the Townsend coefficient of the table formula
//     alpha_P(E) = alpha(E) * (1 + r * sum(nu_exc) / sum(nu_ion)),
// with which the published r values were defined and fitted to gain data
// (Sahin et al., JINST 5 (2010) P05002, eq. 3.1; NIM A 768 (2014) 104).
// Microscopic tracking with the same r gives a lower alpha (1-6 % in Ar:CO2 and
// Ar:CH4), which in a GEM stack compounds to tens of percent in gain.
//
// For every field E: target = alpha_P - eta from a Magboltz table with the
// requested Penning settings (--penning-r auto|<r>, --penning-gas, as for the
// gain program). Microscopic alpha - eta in a uniform field is measured from
// the growth of the ionisation density (direct + Penning) along the drift:
// counts n1, n2 in two windows of width w after a run-in, alpha = ln(n2/n1)/w,
// error from the event-to-event spread (delta method). r_micro is the root of
//     F(r) = mean_E (alpha_micro(E, r) - target(E)) / target(E),
// found from a weighted straight-line fit of F through the r values tried.
//
// Usage (normally launched by gem.py):
//   simgem_penning <gas options> --fields 30000,45000,60000 --events 400
//                  --seed 1 --collisions 5 --out result.json
// Fields in V/cm.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/ComponentConstant.hh>
#include <Garfield/Sensor.hh>

#include "common.hh"

using namespace Garfield;

namespace {

// ionisation windows of the current event, in z (cm); electrons drift to -z
double gW1Top = 0., gW1Bot = 0., gW2Bot = 0.;
long gN1 = 0, gN2 = 0;

void OnIonisation(double, double, double z, double, int, int, Medium*) {
  if (z <= gW1Top && z > gW1Bot) {
    ++gN1;
  } else if (z <= gW1Bot && z > gW2Bot) {
    ++gN2;
  }
}

struct Alpha {
  double value, err;
};

// alpha - eta (1/cm) of microscopic avalanches in a uniform field E (V/cm).
Alpha MeasureAlpha(MediumMagboltz& gas, double E, double alphaGuess, long nEvents) {
  const double run = 1.5 / alphaGuess, w = 2. / alphaGuess, buf = 0.5 / alphaGuess;
  const double L = run + 2. * w + buf;
  // Start strictly inside the area: an electron started on the boundary can
  // make AvalancheMicroscopic::Terminate loop forever.
  const double top = L + 1. / alphaGuess;
  ComponentConstant cmp;
  cmp.SetElectricField(0., 0., E);
  cmp.SetArea(-10., -10., 0., 10., 10., top);
  cmp.SetMedium(&gas);
  Sensor sensor(&cmp);
  sensor.SetArea(-10., -10., 0., 10., 10., top);
  AvalancheMicroscopic aval(&sensor);
  aval.SetUserHandleIonisation(OnIonisation);
  gW1Top = L - run;
  gW1Bot = gW1Top - w;
  gW2Bot = gW1Bot - w;
  double s1 = 0., s2 = 0., s11 = 0., s22 = 0., s12 = 0.;
  for (long k = 0; k < nEvents; ++k) {
    gN1 = gN2 = 0;
    aval.AvalancheElectron(0., 0., L, 0., 0.1, 0., 0., 0.);
    s1 += gN1;
    s2 += gN2;
    s11 += double(gN1) * gN1;
    s22 += double(gN2) * gN2;
    s12 += double(gN1) * gN2;
  }
  if (s1 <= 0. || s2 <= 0.) simgem::Args::Fail("no ionisations in the windows at E = " + std::to_string(E));
  const double n = nEvents;
  const double v1 = s11 / n - (s1 / n) * (s1 / n), v2 = s22 / n - (s2 / n) * (s2 / n);
  const double c12 = s12 / n - (s1 / n) * (s2 / n);
  const double varLn = n * (v1 / (s1 * s1) + v2 / (s2 * s2) - 2. * c12 / (s1 * s2));
  return {std::log(s2 / s1) / w, std::sqrt(std::max(varLn, 0.)) / w};
}

void SetPenning(MediumMagboltz& gas, double r, double lambda, const std::string& pgas) {
  const bool ok = pgas.empty() ? gas.EnablePenningTransfer(r, lambda)
                               : gas.EnablePenningTransfer(r, lambda, pgas);
  if (!ok) simgem::Args::Fail("could not set Penning r = " + std::to_string(r));
}

}  // namespace

int main(int argc, char* argv[]) {
  simgem::Args a(argc, argv);
  std::vector<double> fields;
  for (const auto& s : simgem::Split(a.Str("fields"), ',')) fields.push_back(std::stod(s));
  if (fields.empty()) simgem::Args::Fail("--fields needs at least one value");
  const long nEvents = a.Int("events", 400);
  const int maxEval = a.Int("max-evaluations", 4);
  RandomEngineSTL rng;
  simgem::Seed(rng, a.Int("seed", 1));

  // Target: the table formula with the requested Penning settings.
  auto tab = simgem::MakeGas(a, false);
  tab->SetMaxElectronEnergy(a.Num("max-energy", 200.));
  tab->SetFieldGrid(fields, {0.}, {HalfPi});
  tab->GenerateGasTable(a.Int("collisions", 5), false);
  std::vector<double> alpha0(fields.size()), target(fields.size());
  for (size_t i = 0; i < fields.size(); ++i) {
    tab->ElectronTownsend(0., 0., fields[i], 0., 0., 0., alpha0[i]);
  }
  simgem::ApplyPenning(*tab, a);
  for (size_t i = 0; i < fields.size(); ++i) {
    double al = 0., eta = 0.;
    tab->ElectronTownsend(0., 0., fields[i], 0., 0., 0., al);
    tab->ElectronAttachment(0., 0., fields[i], 0., 0., 0., eta);
    target[i] = al - eta;
  }
  // The r and component the table used (resolves "auto").
  double rTarget = 0., lambda = 0.;
  std::string pgas;
  if (a.Str("penning-r", "off") == "auto" || !a.Str("penning-gas", "").empty()) {
    for (const auto& c : simgem::Split(a.Str("gas"), ',')) {
      const std::string name = simgem::Split(c, ':')[0];
      double r = 0., l = 0.;
      if (tab->MediumGas::GetPenningTransfer(name, r, l) && r > 0.) {
        rTarget = r;
        lambda = l;
        pgas = name;
      }
    }
  } else {
    rTarget = a.Num("penning-r");
    lambda = a.Num("penning-lambda", 0.);
  }
  if (rTarget <= 0.) simgem::Args::Fail("Penning matching needs Penning transfer switched on (r > 0)");
  std::cout << "Target: table Townsend coefficient with r = " << rTarget
            << (pgas.empty() ? " (all components)" : " (" + pgas + ")") << "\n";

  auto micro = simgem::MakeGas(a, false);
  struct Eval {
    double r, F, Ferr;
    std::vector<Alpha> alpha;
  };
  std::vector<Eval> evals;
  auto evaluate = [&](double r) {
    SetPenning(*micro, r, lambda, pgas);
    Eval ev{r, 0., 0., {}};
    double var = 0.;
    for (size_t i = 0; i < fields.size(); ++i) {
      const Alpha al = MeasureAlpha(*micro, fields[i], target[i], nEvents);
      ev.alpha.push_back(al);
      ev.F += (al.value - target[i]) / target[i];
      var += std::pow(al.err / target[i], 2);
      std::cout << "  r = " << std::setprecision(4) << r << "  E = " << fields[i]
                << " V/cm: micro " << al.value << " +- " << al.err << " /cm, table "
                << target[i] << " /cm" << std::endl;
    }
    ev.F /= fields.size();
    ev.Ferr = std::sqrt(var) / fields.size();
    std::cout << "  r = " << r << ": mean relative mismatch " << ev.F << " +- " << ev.Ferr
              << std::endl;
    evals.push_back(ev);
  };
  // Weighted straight line F = c0 + c1 r through the evaluations; root and error.
  double rRoot = rTarget, rErr = 0.;
  auto fit = [&]() {
    double sw = 0., sx = 0., sy = 0., sxx = 0., sxy = 0.;
    for (const auto& e : evals) {
      const double wt = 1. / std::max(e.Ferr * e.Ferr, 1e-12);
      sw += wt; sx += wt * e.r; sy += wt * e.F; sxx += wt * e.r * e.r; sxy += wt * e.r * e.F;
    }
    const double det = sw * sxx - sx * sx;
    const double c1 = (sw * sxy - sx * sy) / det, c0 = (sy - c1 * sx) / sw;
    if (c1 <= 0.) simgem::Args::Fail("microscopic alpha does not rise with r; cannot match");
    rRoot = -c0 / c1;
    // var(root) from the covariance of (c0, c1)
    const double v00 = sxx / det, v11 = sw / det, v01 = -sx / det;
    rErr = std::sqrt(std::max(v00 + rRoot * rRoot * v11 + 2. * rRoot * v01, 0.)) / c1;
  };

  evaluate(rTarget);
  const double step = std::max(0.05, 0.3 * rTarget);
  evaluate(evals[0].F < 0. ? std::min(1., rTarget + step) : std::max(0., rTarget - step));
  fit();
  for (int k = 2; k < maxEval; ++k) {
    const double next = std::clamp(rRoot, 0., 1.);
    const bool close = std::any_of(evals.begin(), evals.end(), [&](const Eval& e) {
      return std::fabs(e.r - next) < std::max(0.01, rErr);
    });
    if (close) break;
    evaluate(next);
    fit();
  }
  const bool clamped = rRoot > 1. || rRoot < 0.;
  const double rMicro = std::clamp(rRoot, 0., 1.);
  std::cout << "RESULT r_micro = " << rMicro << " +- " << rErr << " (table r = " << rTarget
            << (clamped ? "; CLAMPED to [0, 1]" : "") << ")\n";

  std::ofstream out(a.Str("out"));
  if (!out) simgem::Args::Fail("cannot write " + a.Str("out"));
  out << std::setprecision(8) << "{\n  \"r_table\": " << rTarget << ",\n  \"penning_gas\": "
      << (pgas.empty() ? "null" : "\"" + pgas + "\"") << ",\n  \"lambda\": " << lambda
      << ",\n  \"r_micro\": " << rMicro << ",\n  \"r_micro_err\": " << rErr
      << ",\n  \"clamped\": " << (clamped ? "true" : "false")
      << ",\n  \"events_per_field\": " << nEvents << ",\n  \"fields_V_cm\": [";
  for (size_t i = 0; i < fields.size(); ++i) out << (i ? ", " : "") << fields[i];
  out << "],\n  \"alpha_table_no_penning\": [";
  for (size_t i = 0; i < fields.size(); ++i) out << (i ? ", " : "") << alpha0[i];
  out << "],\n  \"alpha_table_target\": [";
  for (size_t i = 0; i < fields.size(); ++i) out << (i ? ", " : "") << target[i];
  out << "],\n  \"evaluations\": [";
  for (size_t k = 0; k < evals.size(); ++k) {
    const auto& e = evals[k];
    out << (k ? ",\n" : "\n") << "    {\"r\": " << e.r << ", \"mismatch\": " << e.F
        << ", \"mismatch_err\": " << e.Ferr << ", \"alpha_micro\": [";
    for (size_t i = 0; i < e.alpha.size(); ++i) out << (i ? ", " : "") << e.alpha[i].value;
    out << "], \"alpha_micro_err\": [";
    for (size_t i = 0; i < e.alpha.size(); ++i) out << (i ? ", " : "") << e.alpha[i].err;
    out << "]}";
  }
  out << "\n  ]\n}\n";
  return 0;
}
