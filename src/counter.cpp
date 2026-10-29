// counter.cpp -- gain of a cylindrical proportional counter from a Magboltz
// table, the way published Penning transfer probabilities were fitted
// (Sahin et al., JINST 5 (2010) P05002, eq. 3.1):
//     G = exp( integral_{r_a}^{r_c} alpha(E(r)) dr ),  E(r) = V / (r ln(r_c / r_a)),
// with the Townsend coefficient of the table, without and with Penning transfer
// (--penning-*, as for the other programs). Used to check the gas model
// against measured counter gains (docs/validation.md).
//
// Usage:
//   simgem_counter <gas options> --anode-radius-cm 0.00125 --cathode-radius-cm 1.5
//                  --voltages 1690,1883,2059 [--table gas/x.gas] [--collisions 10]
//                  [--out result.csv]
// The table covers the fields of the counter at those voltages; with --table it
// is read if the file exists and written otherwise.
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "common.hh"

using namespace Garfield;

namespace {

// log G: alpha integrated over log-spaced radii from the anode to the cathode
double LogGain(MediumMagboltz& gas, double V, double ra, double rc) {
  const double L = std::log(rc / ra);
  constexpr int n = 20000;
  double s = 0.;
  for (int i = 0; i < n; ++i) {
    const double r0 = ra * std::pow(rc / ra, double(i) / n);
    const double r1 = ra * std::pow(rc / ra, double(i + 1) / n);
    double alpha = 0.;
    gas.ElectronTownsend(0., 0., V / (std::sqrt(r0 * r1) * L), 0., 0., 0., alpha);
    s += alpha * (r1 - r0);
  }
  return s;
}

}  // namespace

int main(int argc, char* argv[]) {
  simgem::Args a(argc, argv);
  const double ra = a.Num("anode-radius-cm");
  const double rc = a.Num("cathode-radius-cm");
  if (!(ra > 0. && rc > ra)) simgem::Args::Fail("need 0 < anode radius < cathode radius");
  std::vector<double> V;
  for (const auto& s : simgem::Split(a.Str("voltages"), ',')) V.push_back(std::stod(s));
  if (V.empty()) simgem::Args::Fail("--voltages needs at least one value");
  const std::string table = a.Str("table", "");

  auto gas = simgem::MakeGas(a, false);
  if (!table.empty() && std::ifstream(table)) {
    if (!gas->LoadGasFile(table)) simgem::Args::Fail("could not read " + table);
    std::cout << "Gas table: " << table << "\n";
  } else {
    const double L = std::log(rc / ra);
    double vmin = V[0], vmax = V[0];
    for (double v : V) {
      vmin = std::min(vmin, v);
      vmax = std::max(vmax, v);
    }
    // alpha is negligible below ~3 kV/cm; cover up to the field on the wire
    const double emin = std::max(3.e3, vmin / (rc * L)), emax = 1.05 * vmax / (ra * L);
    gas->SetMaxElectronEnergy(a.Num("max-energy", 200.));
    gas->SetFieldGrid(emin, emax, a.Int("field-points", 28), true);
    gas->GenerateGasTable(a.Int("collisions", 10), false);
    if (!table.empty() && !gas->WriteGasFile(table)) simgem::Args::Fail("could not write " + table);
  }
  std::vector<double> g0;
  for (double v : V) g0.push_back(std::exp(LogGain(*gas, v, ra, rc)));
  simgem::ApplyPenning(*gas, a);

  std::ofstream out;
  if (a.Has("out")) {
    out.open(a.Str("out"));
    if (!out) simgem::Args::Fail("cannot write " + a.Str("out"));
    out << "V,gain_no_penning,gain_penning\n";
  }
  for (size_t i = 0; i < V.size(); ++i) {
    const double g = std::exp(LogGain(*gas, V[i], ra, rc));
    std::cout << "RESULT V = " << V[i] << " V: gain " << g0[i] << " without Penning, " << g
              << " with Penning\n";
    if (out) out << V[i] << ',' << g0[i] << ',' << g << '\n';
  }
  return 0;
}
