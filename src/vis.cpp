// vis.cpp
// -----------------------------------------------------------------------------
// Same physics as gain.cpp (the full stack, one Elmer field map),
// but instead of only counting electrons it DUMPS THE TRAJECTORIES so they can
// be plotted later. Three files per run:
//
//   vis_tracks<tag>.csv    one row per stored trajectory point
//   vis_summary<tag>.csv   one row per electron (endpoints, status, weight)
//   vis_meta<tag>.json     cell bounding box + start/readout planes, so the
//                          plotting script can draw the geometry
//
// Why a step callback instead of Garfield's stored drift lines:
// AvalancheMicroscopic::EnableDriftLines() keeps EVERY collision of EVERY
// electron in memory (SetCollisionSteps is only honoured in the B-field
// version). At 8 bar that is ~10^6-10^7 points per drift line, i.e. hundreds of
// MB for a handful of tracks and tens of GB for a real avalanche. So drift-line
// storage stays OFF (Garfield then keeps only the endpoints, which is exactly
// what the summary needs) and the trajectory is streamed to disk from
// SetUserHandleStep, thinned on the fly:
//
//   min_step    write a point only once the electron has moved this far (cm)
//               from the last point written -> file resolution set in cm, not
//               in collisions
//   max_tracks  write full paths for the first N electrons of an event only;
//               every electron still gets a summary row
//   size_limit  cap the avalanche itself
//
// The callback is a plain function pointer with no track id, so drift-line
// boundaries have to be inferred. Two signals do it: between two calls within
// one line the electron makes exactly ONE free flight, so it cannot cover more
// ground than sqrt(2E/m)*dt, nor take longer than dt_max (a very loose upper
// bound on a single flight — argon's Ramsauer minimum makes these surprisingly
// long, hence the 0.03 ns default rather than the ~1e-3 ns mean). A new line
// starts elsewhere, or earlier in time, and breaks both bounds.
//
// The detected lines are then matched back to Garfield's electrons by their
// seed point, which is bit-identical, so the `track` column joins the two CSVs.
// The split is not guaranteed exact, and the program reports both failure
// modes per event (on the baseline stack at r=0.6, ~426 electrons: 2 fragments,
// 4 merges):
//   fragments  one electron cut into several tracks. Harmless for plotting, and
//              exactly repairable: a track id that does NOT appear in the
//              summary continues the highest summary track id below it.
//   merges     two electrons drawn as one track; the second one gets
//              track = -1 in the summary. Lower dt_max to trade these for
//              fragments.
//
// Usage (normally launched by `gem.py vis`):
//   simgem_vis --field fields/<hash> <gas options> --events 3 --seed 1000
//              --z-start 0.7 --z-readout 0.001 --outdir runs/<cfg>
//              [--tag ""] [--max-tracks 200] [--size-limit 2000]
//              [--min-step 5e-4] [--dt-max 0.03]
// -----------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/FundamentalConstants.hh>
#include <Garfield/Sensor.hh>

#include "common.hh"

using namespace Garfield;

namespace {

// --- state shared with the step callback (it takes no user pointer) ----------
std::ofstream* gTrk = nullptr;
int gEvent = 0;
long gTrack = -1;          // index of the drift line being written
long gMaxTracks = 200;     // < 0 : all
double gMinStep2 = 0.;     // squared write threshold (cm^2)
double gDtMax = 0.01;      // ns; a longer gap than any single free flight
double gLastT = 0.;        // previous step (boundary detection)
double gLastE = 0.;
double gLx = 0., gLy = 0., gLz = 0.;
double gWx = 0., gWy = 0., gWz = 0.;   // last point written
long gNWritten = 0;        // points written for the current track
long gTotalWritten = 0;
std::vector<std::array<double, 4> > gStarts;   // seed of each detected track

// sqrt(2E/m) in cm/ns per sqrt(eV): the fastest an electron of energy E can go.
const double kVPerSqrtEv = SpeedOfLight * std::sqrt(2. / ElectronMass);
const double kSafety = 3.;   // slack for acceleration during the flight

void ResetTracking() {
  gTrack = -1;
  gNWritten = 0;
  gLastT = 0.;
  gLastE = 0.;
  gStarts.clear();
}

void StepHandler(double x, double y, double z, double t, double e,
                 double /*dx*/, double /*dy*/, double /*dz*/, bool hole) {
  if (hole) return;
  // One free flight separates two calls within a drift line, so the distance
  // covered is bounded by the speed at the higher of the two energies. Anything
  // beyond that (or a step back in time) means a new drift line has started.
  const double dt = t - gLastT;
  const double sx = x - gLx, sy = y - gLy, sz = z - gLz;
  const double reach = kSafety * kVPerSqrtEv * std::sqrt(std::max(e, gLastE)) * dt;
  const bool newTrack = gTrack < 0 || dt <= 0. || dt > gDtMax ||
                        (sx * sx + sy * sy + sz * sz) > reach * reach;
  if (newTrack) {
    ++gTrack;
    gNWritten = 0;
    gStarts.push_back({x, y, z, t});
  }
  gLx = x; gLy = y; gLz = z; gLastT = t; gLastE = e;

  if (gMaxTracks >= 0 && gTrack >= gMaxTracks) return;
  if (gNWritten > 0) {
    const double dx = x - gWx, dy = y - gWy, dz = z - gWz;
    if (dx * dx + dy * dy + dz * dz < gMinStep2) return;
  }
  *gTrk << gEvent << ',' << gTrack << ',' << gNWritten << ',' << x << ',' << y
        << ',' << z << ',' << t << ',' << e << '\n';
  gWx = x; gWy = y; gWz = z;
  ++gNWritten;
  ++gTotalWritten;
}

}  // namespace

int main(int argc, char* argv[]) {
  simgem::Args a(argc, argv);
  const std::string outdir = a.Str("outdir");
  const std::string tag = a.Str("tag", "");
  const long nEvents = a.Int("events", 1);
  const long maxTracks = a.Int("max-tracks", 200);
  const long sizeLimit = a.Int("size-limit", 2000);
  const double minStep = a.Num("min-step", 5e-4);  // cm
  const double dtMax = a.Num("dt-max", 0.03);      // ns
  const double zStart = a.Num("z-start");
  const double zReadout = a.Num("z-readout");
  const double e0 = a.Num("primary-energy", 0.1);

  RandomEngineSTL rng;
  simgem::Seed(rng, a.Int("seed", 1000));
  auto gas = simgem::MakeGas(a);
  double cxmin, cymin, czmin, cxmax, cymax, czmax;
  auto fm = simgem::LoadField(a.Str("field"), gas.get(), cxmin, cymin, czmin,
                              cxmax, cymax, czmax);
  fm->PrintRange();
  const double px = cxmax - cxmin;
  const double py = cymax - cymin;
  std::cout << "Starting electrons at z = " << zStart << " cm, counting those below z = "
            << zReadout << " cm\nWriting a point every " << minStep
            << " cm, full paths for up to " << maxTracks
            << " electrons/event, avalanche size limit " << sizeLimit << "\n";

  Sensor sensor(fm.get());
  sensor.SetArea(cxmin - 1000. * px, cymin - 1000. * py, czmin,
                 cxmax + 1000. * px, cymax + 1000. * py, czmax);
  AvalancheMicroscopic aval(&sensor);
  aval.EnableAvalancheSizeLimit(sizeLimit);
  aval.EnablePathLengthComputation();
  aval.SetUserHandleStep(StepHandler);
  // NOTE: EnableDriftLines() stays off on purpose -- see the header comment.

  {
    std::ofstream meta(outdir + "/vis_meta" + tag + ".json");
    meta << "{\n"
         << "  \"units\": {\"length\": \"cm\", \"time\": \"ns\", \"energy\": \"eV\"},\n"
         << "  \"cell\": {\"xmin\": " << cxmin << ", \"xmax\": " << cxmax
         << ", \"ymin\": " << cymin << ", \"ymax\": " << cymax
         << ", \"zmin\": " << czmin << ", \"zmax\": " << czmax << "},\n"
         << "  \"z_start\": " << zStart << ",\n"
         << "  \"z_readout\": " << zReadout << ",\n"
         << "  \"n_events\": " << nEvents << ",\n"
         << "  \"max_tracks\": " << maxTracks << ",\n"
         << "  \"size_limit\": " << sizeLimit << ",\n"
         << "  \"min_step_cm\": " << minStep << ",\n"
         << "  \"dt_max_ns\": " << dtMax << "\n"
         << "}\n";
  }

  std::ofstream trk(outdir + "/vis_tracks" + tag + ".csv");
  trk << "event,track,point,x_cm,y_cm,z_cm,t_ns,e_eV\n";
  std::ofstream sum(outdir + "/vis_summary" + tag + ".csv");
  sum << "event,track,status,weight,path_length_cm,"
         "x0_cm,y0_cm,z0_cm,t0_ns,x1_cm,y1_cm,z1_cm,t1_ns,reached_readout\n";
  trk.precision(7);
  sum.precision(7);

  gTrk = &trk;
  gMaxTracks = maxTracks;
  gMinStep2 = minStep * minStep;
  gDtMax = dtMax;

  long nUnmatched = 0;
  for (long i = 0; i < nEvents; ++i) {
    gEvent = i;
    ResetTracking();

    const double x0 = cxmin + px * RndmUniform();
    const double y0 = cymin + py * RndmUniform();
    aval.AvalancheElectron(x0, y0, zStart, 0., e0, 0., 0., 0.);

    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);

    // Map each detected drift line back to Garfield's electron list via its
    // seed point (the callback sees the very same doubles).
    std::map<std::array<double, 4>, long> byStart;
    for (size_t j = 0; j < gStarts.size(); ++j) byStart.emplace(gStarts[j], (long)j);

    long nReadout = 0, nMissed = 0, nElectrons = 0;
    for (const auto& el : aval.GetElectrons()) {
      if (el.path.empty()) continue;
      ++nElectrons;
      const auto& p0 = el.path.front();
      const auto& p1 = el.path.back();
      const bool reached = p1.z < zReadout;
      if (reached) ++nReadout;
      const auto it = byStart.find({p0.x, p0.y, p0.z, p0.t});
      const long track = (it == byStart.end()) ? -1 : it->second;
      if (track < 0) ++nMissed;
      sum << i << ',' << track << ',' << el.status << ',' << el.weight << ','
          << el.pathLength << ',' << p0.x << ',' << p0.y << ',' << p0.z << ',' << p0.t
          << ',' << p1.x << ',' << p1.y << ',' << p1.z << ',' << p1.t
          << ',' << (reached ? 1 : 0) << '\n';
    }
    nUnmatched += nMissed;
    // Lines with no electron of their own are fragments of the preceding one.
    const long nFragments = (long)gStarts.size() - (nElectrons - nMissed);

    std::cout << "event " << (i + 1) << "/" << nEvents << ": " << ne
              << " electrons (" << nReadout << " at readout), "
              << gStarts.size() << " tracks, "
              << gTotalWritten << " points written so far";
    if (nFragments || nMissed) {
      std::cout << "  [" << nFragments << " fragments, " << nMissed << " merged]";
    }
    std::cout << std::endl;
  }
  trk.close();
  sum.close();

  if (nUnmatched) {
    std::cout << "NOTE: " << nUnmatched << " electron(s) were merged into the"
              << " preceding track (track = -1 in the summary). Lower dt_max ("
              << dtMax << " ns) to split more aggressively.\n";
  }
  std::cout << "Trajectories: " << outdir << "/vis_tracks" << tag << ".csv"
            << "  (" << gTotalWritten << " points)\n";
  return 0;
}
