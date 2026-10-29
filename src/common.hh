// common.hh -- shared pieces of the SimGEM programs: command-line parsing,
// gas set-up and field-map loading.
//
// All programs take named options (--key value). The Python driver (gem.py)
// builds these from the config, so the programs need no geometry knowledge.
//
// Gas options:
//   --gas ar:98,ch4:2       up to 6 components (any names Garfield knows)
//   --pressure-torr 6000.5  --temperature-k 293.15
//   --penning-r 0.18|auto|off   [--penning-lambda 0] [--penning-gas ar]
//
// The microscopic avalanche needs NO Magboltz gas table: MediumMagboltz
// builds its collision-rate tables from the cross sections when tracking
// starts. Tables (gastable program) are only for drift-line/Townsend methods.
#pragma once

#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Garfield/ComponentElmer.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/Random.hh>
#include <Garfield/RandomEngineSTL.hh>

namespace simgem {

class Args {
 public:
  Args(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string k = argv[i];
      if (k.rfind("--", 0) != 0) Fail("unexpected argument '" + k + "'");
      k = k.substr(2);
      if (i + 1 >= argc) Fail("option --" + k + " needs a value");
      m_kv[k] = argv[++i];
    }
  }
  bool Has(const std::string& k) const { return m_kv.count(k) > 0; }
  std::string Str(const std::string& k) const {
    auto it = m_kv.find(k);
    if (it == m_kv.end()) Fail("missing required option --" + k);
    return it->second;
  }
  std::string Str(const std::string& k, const std::string& def) const {
    return Has(k) ? Str(k) : def;
  }
  double Num(const std::string& k) const { return std::stod(Str(k)); }
  double Num(const std::string& k, double def) const { return Has(k) ? Num(k) : def; }
  long Int(const std::string& k) const { return std::stol(Str(k)); }
  long Int(const std::string& k, long def) const { return Has(k) ? Int(k) : def; }

  [[noreturn]] static void Fail(const std::string& msg) {
    std::cerr << "ERROR: " << msg << "\n";
    std::exit(2);
  }

 private:
  std::map<std::string, std::string> m_kv;
};

inline std::vector<std::string> Split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, sep)) {
    if (!item.empty()) out.push_back(item);
  }
  return out;
}

// Penning transfer from --penning-r/--penning-lambda/--penning-gas.
inline void ApplyPenning(Garfield::MediumMagboltz& gas, const Args& a) {
  const std::string r = a.Str("penning-r", "off");
  const double lambda = a.Num("penning-lambda", 0.);
  const std::string pgas = a.Str("penning-gas", "");
  bool ok = true;
  if (r == "auto") {
    ok = gas.EnablePenningTransfer();
    std::cout << "Penning: Garfield's built-in parameterisation\n";
  } else if (r != "off") {
    const double rv = std::stod(r);
    ok = pgas.empty() ? gas.EnablePenningTransfer(rv, lambda)
                      : gas.EnablePenningTransfer(rv, lambda, pgas);
    std::cout << "Penning: r = " << rv << ", lambda = " << lambda
              << (pgas.empty() ? " (all components)" : " (" + pgas + ")") << "\n";
  } else {
    std::cout << "Penning: off\n";
  }
  if (!ok) Args::Fail("could not enable Penning transfer (see Garfield message above)");
}

// Gas mixture from --gas/--pressure-torr/--temperature-k and, unless
// `penning` is false, --penning-*.
inline std::unique_ptr<Garfield::MediumMagboltz> MakeGas(const Args& a, bool penning = true) {
  std::vector<std::string> names(6, "");
  std::vector<double> fr(6, 0.);
  const auto comps = Split(a.Str("gas"), ',');
  if (comps.empty() || comps.size() > 6) Args::Fail("--gas needs 1 to 6 name:fraction entries");
  for (size_t i = 0; i < comps.size(); ++i) {
    const auto nf = Split(comps[i], ':');
    if (nf.size() != 2) Args::Fail("bad --gas entry '" + comps[i] + "'");
    names[i] = nf[0];
    fr[i] = std::stod(nf[1]);
  }
  auto gas = std::make_unique<Garfield::MediumMagboltz>();
  if (!gas->SetComposition(names[0], fr[0], names[1], fr[1], names[2], fr[2],
                           names[3], fr[3], names[4], fr[4], names[5], fr[5])) {
    Args::Fail("Garfield rejected the gas composition '" + a.Str("gas") + "'");
  }
  gas->SetPressure(a.Num("pressure-torr"));
  gas->SetTemperature(a.Num("temperature-k"));
  if (penning) {
    ApplyPenning(*gas, a);
    gas->PrintGas();
  }
  return gas;
}

// Field map from a SimGEM field directory, translational-periodic in x and y.
// Returns the cell bounding box (cm) through the reference arguments.
inline std::unique_ptr<Garfield::ComponentElmer> LoadField(
    const std::string& d, Garfield::Medium* gas, double& x0, double& y0,
    double& z0, double& x1, double& y1, double& z1) {
  auto fm = std::make_unique<Garfield::ComponentElmer>(
      d + "/stack/mesh.header", d + "/stack/mesh.elements",
      d + "/stack/mesh.nodes", d + "/dielectrics.dat",
      d + "/stack/case.result", "mm");
  if (!fm->IsReady()) Args::Fail("could not load the field map in " + d);
  // Cell extent BEFORE enabling periodicity (afterwards x/y are unbounded).
  fm->GetBoundingBox(x0, y0, z0, x1, y1, z1);
  fm->EnablePeriodicityX();
  fm->EnablePeriodicityY();
  fm->SetGas(gas);
  return fm;
}

inline void Seed(Garfield::RandomEngineSTL& rng, unsigned long seed) {
  rng.SetSeed(seed);
  Garfield::Random::SetEngine(rng);
}

}  // namespace simgem
