// gain.cpp -- one worker of the parallel stack-gain simulation.
//
// Single primary electrons start at --z-start (uniform over the cell in x,y)
// and are followed microscopically through the whole stack. Per event:
//   absolute_gain   all electrons produced (incl. the primary)
//   effective_gain  electrons that end below --z-readout (reach the readout)
//   n_ions          ions produced
//   truncated       1 if the avalanche-size limit was hit (gain undercounted)
//
// Usage (normally launched by gem.py):
//   simgem_gain --field fields/<hash> --gas ar:98,ch4:2 --pressure-torr 6000.5
//               --temperature-k 293.15 --penning-r 0.18 --penning-gas ar
//               --events 100 --seed 12345 --z-start 0.7 --z-readout 0.001
//               --size-limit 500000 --out events.csv [--progress file]
//               [--heartbeat file]
// Lengths in cm.
//
// Heartbeat: the number of electron collisions so far, written to the
// heartbeat file every ~2 s from the step callback. A long avalanche keeps it
// moving; a hang inside Garfield (e.g. AvalancheMicroscopic::Terminate looping
// on a point exactly on a dielectric boundary) does not, so the parent can
// tell the two apart, kill the worker and restart it.
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/Sensor.hh>

#include "common.hh"

using namespace Garfield;

namespace {

std::string gHeartbeat;
long long gSteps = 0;
auto gLastBeat = std::chrono::steady_clock::now();

void Beat() {
  std::ofstream hb(gHeartbeat, std::ios::trunc);
  hb << gSteps << '\n';
  gLastBeat = std::chrono::steady_clock::now();
}

void StepHandler(double, double, double, double, double, double, double, double, bool) {
  if ((++gSteps & 0xFFFF) != 0) return;   // look at the clock every 65536 collisions
  if (std::chrono::steady_clock::now() - gLastBeat > std::chrono::seconds(2)) Beat();
}

}  // namespace

int main(int argc, char* argv[]) {
  simgem::Args a(argc, argv);
  const long nEvents = a.Int("events");
  const double zStart = a.Num("z-start");
  const double zReadout = a.Num("z-readout");
  const long sizeLimit = a.Int("size-limit", 500000);
  const double e0 = a.Num("primary-energy", 0.1);
  const std::string progress = a.Str("progress", "");
  gHeartbeat = a.Str("heartbeat", "");

  RandomEngineSTL rng;
  simgem::Seed(rng, a.Int("seed"));
  auto gas = simgem::MakeGas(a);

  double x0, y0, z0, x1, y1, z1;
  auto fm = simgem::LoadField(a.Str("field"), gas.get(), x0, y0, z0, x1, y1, z1);
  fm->PrintRange();
  if (zStart <= z0 || zStart >= z1) simgem::Args::Fail("--z-start is outside the field map");
  const double px = x1 - x0, py = y1 - y0;

  Sensor sensor(fm.get());
  // The periodic field is defined everywhere laterally; give electrons room.
  sensor.SetArea(x0 - 1000. * px, y0 - 1000. * py, z0, x1 + 1000. * px, y1 + 1000. * py, z1);
  AvalancheMicroscopic aval(&sensor);
  aval.EnableAvalancheSizeLimit(sizeLimit);
  if (!gHeartbeat.empty()) {
    aval.SetUserHandleStep(StepHandler);
    Beat();
  }

  std::ofstream csv(a.Str("out"));
  if (!csv) simgem::Args::Fail("cannot write " + a.Str("out"));
  csv << "event,x0_cm,y0_cm,absolute_gain,effective_gain,n_ions,truncated\n";
  csv.precision(7);

  double sum = 0.;
  for (long i = 0; i < nEvents; ++i) {
    const double xs = x0 + px * RndmUniform();
    const double ys = y0 + py * RndmUniform();
    aval.AvalancheElectron(xs, ys, zStart, 0., e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    long nRead = 0;
    for (const auto& el : aval.GetElectrons()) {
      if (!el.path.empty() && el.path.back().z < zReadout) ++nRead;
    }
    const int trunc = (ne >= sizeLimit) ? 1 : 0;
    csv << i << ',' << xs << ',' << ys << ',' << ne << ',' << nRead << ',' << ni
        << ',' << trunc << '\n';
    csv.flush();
    sum += nRead;
    if (!progress.empty()) {
      std::ofstream pf(progress, std::ios::trunc);
      pf << (i + 1) << ' ' << nEvents << '\n';
    }
    if (!gHeartbeat.empty()) Beat();
  }
  std::cout << "Mean effective gain (this worker): " << sum / std::max(nEvents, 1L)
            << " over " << nEvents << " events\n";
  return 0;
}
