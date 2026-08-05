// vis_stack.cpp
// -----------------------------------------------------------------------------
// Same physics as gain_stack.cpp (full triple-GEM stack, one Elmer field map),
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
// Usage:
//   ./vis_stack <outdir> <gasfile> <penning_r> <n_events>
//               [seed] [tag] [max_tracks] [size_limit] [min_step_cm] [dt_max_ns]
//   e.g.  ./vis_stack runs/baseline gas/ar98_ch4_2_8bar.gas 0.18 3
//
// As in gain_stack.cpp the z-range comes from the field map's bounding box, so
// no geometry parameters are needed here.
// -----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <array>
#include <map>
#include <vector>

#include <Garfield/MediumMagboltz.hh>
#include <Garfield/ComponentElmer.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/Random.hh>
#include <Garfield/RandomEngineSTL.hh>
#include <Garfield/FundamentalConstants.hh>

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
  const std::string outdir = (argc > 1) ? argv[1] : "stack_run";
  const std::string gasfile = (argc > 2) ? argv[2] : "ar98_ch4_2_8bar.gas";
  const double rPenning = (argc > 3) ? std::atof(argv[3]) : 0.18;
  const int nEvents = (argc > 4) ? std::atoi(argv[4]) : 1;
  const unsigned int seed = (argc > 5) ? (unsigned int)std::strtoul(argv[5], nullptr, 10) : 0;
  const std::string tag = (argc > 6) ? argv[6] : "";
  const long maxTracks = (argc > 7) ? std::atol(argv[7]) : 200;
  const size_t sizeLimit = (argc > 8) ? (size_t)std::strtoul(argv[8], nullptr, 10) : 2000;
  const double minStep = (argc > 9) ? std::atof(argv[9]) : 5e-4;  // cm (5 um)
  const double dtMax = (argc > 10) ? std::atof(argv[10]) : 0.03;  // ns

  RandomEngineSTL rng;
  if (seed) {
    rng.SetSeed(seed);
    Random::SetEngine(rng);
  }

  // ---- gas -------------------------------------------------------------------
  MediumMagboltz gas;
  gas.LoadGasFile(gasfile);
  if (rPenning < 0.) {
    gas.EnablePenningTransfer();
    std::cout << "Penning: built-in parameterization (a priori, from concentration+pressure)\n";
  } else {
    gas.EnablePenningTransfer(rPenning, 0., "ar");
  }

  // ---- field map -------------------------------------------------------------
  ComponentElmer fm(outdir + "/stack/mesh.header",
                    outdir + "/stack/mesh.elements",
                    outdir + "/stack/mesh.nodes",
                    outdir + "/dielectrics.dat",
                    outdir + "/stack/case.result",
                    "mm");

  double cxmin, cymin, czmin, cxmax, cymax, czmax;
  fm.GetBoundingBox(cxmin, cymin, czmin, cxmax, cymax, czmax);

  // See gain_stack.cpp: mirror (legacy aligned) vs translational (misaligned).
  std::string pmode = "mirror";
  {
    std::ifstream pin(outdir + "/periodicity");
    if (pin) std::getline(pin, pmode);
    while (!pmode.empty() &&
           (pmode.back() == '\n' || pmode.back() == '\r' || pmode.back() == ' '))
      pmode.pop_back();
  }
  if (pmode == "translate") {
    fm.EnablePeriodicityX();
    fm.EnablePeriodicityY();
    std::cout << "Periodicity: translational (misalignment-capable)\n";
  } else {
    fm.EnableMirrorPeriodicityX();
    fm.EnableMirrorPeriodicityY();
    std::cout << "Periodicity: mirror (aligned cell)\n";
  }
  fm.SetGas(&gas);
  fm.PrintRange();

  const double px = cxmax - cxmin;
  const double py = cymax - cymin;
  const double dz = czmax - czmin;
  const double zStart = czmax - 0.02 * dz;
  const double zReadout = czmin + 0.01 * dz;
  std::cout << "Cell (cm): x[" << cxmin << "," << cxmax << "] y[" << cymin
            << "," << cymax << "] z[" << czmin << "," << czmax << "]\n"
            << "Starting electrons at z = " << zStart
            << ", counting those below z = " << zReadout << "\n"
            << "Writing a point every " << minStep << " cm, full paths for up to "
            << maxTracks << " electrons/event, avalanche size limit " << sizeLimit
            << "\n";

  Sensor sensor(&fm);
  const double pad = 40.0;
  sensor.SetArea(cxmin - pad * px, cymin - pad * py, czmin,
                 cxmax + pad * px, cymax + pad * py, czmax);

  AvalancheMicroscopic aval(&sensor);
  aval.EnableAvalancheSizeLimit(sizeLimit);
  aval.EnablePathLengthComputation();
  aval.SetUserHandleStep(StepHandler);
  // NOTE: EnableDriftLines() stays off on purpose -- see the header comment.

  // Geometry/meta for the plotting script.
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
         << "  \"penning_r\": " << rPenning << ",\n"
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
  trk.precision(6);
  sum.precision(6);

  gTrk = &trk;
  gMaxTracks = maxTracks;
  gMinStep2 = minStep * minStep;
  gDtMax = dtMax;

  long nUnmatched = 0;
  for (int i = 0; i < nEvents; ++i) {
    gEvent = i;
    ResetTracking();

    const double x0 = cxmin + px * RndmUniform();
    const double y0 = cymin + py * RndmUniform();
    aval.AvalancheElectron(x0, y0, zStart, 0., 0.1, 0., 0., 0.);

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

    // progress file, same convention as gain_stack
    std::ofstream pf(outdir + "/.progress" + tag, std::ios::trunc);
    pf << (i + 1) << ' ' << nEvents << '\n';
    pf.close();
  }
  trk.close();
  sum.close();

  if (nUnmatched) {
    std::cout << "NOTE: " << nUnmatched << " electron(s) were merged into the"
              << " preceding track (track = -1 in the summary). Lower dt_max ("
              << dtMax << " ns) to split more aggressively.\n";
  }
  std::cout << "\n=================================================\n"
            << "Trajectories: " << outdir << "/vis_tracks" << tag << ".csv"
            << "  (" << gTotalWritten << " points)\n"
            << "Per-track:    " << outdir << "/vis_summary" << tag << ".csv\n"
            << "Geometry:     " << outdir << "/vis_meta" << tag << ".json\n"
            << "=================================================\n";
  return 0;
}
