// transport.cpp -- microscopic electron transport in uniform fields.
//
// For each field, N electrons are drifted (no multiplication) from z = d down
// to z = 0 and the arrival spread gives drift velocity and diffusion. Used by
// gem.py to judge how much memory of one foil's hole pattern survives a
// transfer gap, i.e. whether foil misalignment can matter at all.
//
// Usage:
//   simgem_transport <gas options> --fields 278,2857 (V/cm) --distance 0.05 (cm)
//                    --electrons 300 --seed 1 --out transport.csv
#include <cmath>
#include <fstream>
#include <iostream>

#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/ComponentConstant.hh>
#include <Garfield/Sensor.hh>

#include "common.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  simgem::Args a(argc, argv);
  const double d = a.Num("distance");
  const long n = a.Int("electrons", 300);
  RandomEngineSTL rng;
  simgem::Seed(rng, a.Int("seed", 1));
  auto gas = simgem::MakeGas(a);

  std::ofstream out(a.Str("out"));
  out << "field_V_cm,distance_cm,electrons,arrived,v_drift_cm_ns,"
         "D_T_sqrt_cm,D_L_sqrt_cm,mean_energy_eV\n";
  out.precision(6);
  for (const auto& fs : simgem::Split(a.Str("fields"), ',')) {
    const double E = std::stod(fs);
    ComponentConstant cmp;
    cmp.SetElectricField(0., 0., E);  // field along +z -> electrons drift to -z
    cmp.SetArea(-10., -10., 0., 10., 10., 10. * d);
    cmp.SetMedium(gas.get());
    Sensor sensor(&cmp);
    sensor.SetArea(-10., -10., 0., 10., 10., 10. * d);
    AvalancheMicroscopic drift(&sensor);

    double sx = 0., sxx = 0., sy = 0., syy = 0., st = 0., stt = 0., se = 0.;
    long m = 0;
    for (long i = 0; i < n; ++i) {
      drift.DriftElectron(0., 0., d, 0., 0.1, 0., 0., 0.);
      const auto& el = drift.GetElectrons().front();
      const auto& p = el.path.back();
      if (p.z > 1e-6) continue;  // did not reach the plane (attached, ...)
      ++m;
      sx += p.x; sxx += p.x * p.x;
      sy += p.y; syy += p.y * p.y;
      st += p.t; stt += p.t * p.t;
      se += p.energy;
    }
    if (m < 2) {
      out << E << ',' << d << ',' << n << ',' << m << ",nan,nan,nan,nan\n";
      continue;
    }
    const double vx = sxx / m - (sx / m) * (sx / m);
    const double vy = syy / m - (sy / m) * (sy / m);
    const double vt = stt / m - (st / m) * (st / m);
    const double tm = st / m;
    const double vd = d / tm;
    const double dt = std::sqrt(0.5 * (vx + vy) / d);
    const double dl = std::sqrt(std::max(vt, 0.)) * vd / std::sqrt(d);
    out << E << ',' << d << ',' << n << ',' << m << ',' << vd << ',' << dt << ','
        << dl << ',' << se / m << '\n';
    std::cout << "E = " << E << " V/cm: v_d = " << vd * 1e3 << " um/ns, D_T = "
              << dt * 1e4 << " um/sqrt(cm), D_L = " << dl * 1e4 << " um/sqrt(cm) ("
              << m << "/" << n << " arrived)" << std::endl;
  }
  return 0;
}
