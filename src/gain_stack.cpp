// gain_stack.cpp
// -----------------------------------------------------------------------------
// Effective gain of a FULL triple-GEM stack from a single Elmer field map built
// by build_stack.py. Primaries start just below the drift plane; electrons that
// reach the readout plane are counted -> stack gain (collection + all three
// avalanches + all transfers + extraction, in one number).
//
// Usage:
//   ./gain_stack <outdir> <gasfile> <penning_r> <n_events>
//   e.g.  ./gain_stack stack_run ar98_ch4_2_8bar.gas 0.18 200
//
// The z-range is read from the field map's bounding box, so this program needs
// no geometry parameters — change the stack in config.json and it adapts.
// -----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <fstream>

#include <Garfield/MediumMagboltz.hh>
#include <Garfield/ComponentElmer.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/Random.hh>
#include <Garfield/RandomEngineSTL.hh>

using namespace Garfield;

int main(int argc, char* argv[]) {
  const std::string outdir = (argc > 1) ? argv[1] : "stack_run";
  const std::string gasfile = (argc > 2) ? argv[2] : "ar98_ch4_2_8bar.gas";
  const double rPenning = (argc > 3) ? std::atof(argv[3]) : 0.18;
  const int nEvents = (argc > 4) ? std::atoi(argv[4]) : 200;
  const unsigned int seed = (argc > 5) ? (unsigned int)std::strtoul(argv[5], nullptr, 10) : 0;
  const std::string tag = (argc > 6) ? argv[6] : "";

  // Distinct RNG per parallel worker. Recent Garfield (>=0.3) uses the Random
  // class with a pluggable engine, not a global 'randomEngine' object. The
  // engine must outlive its use, so keep it in main scope.
  RandomEngineSTL rng;
  if (seed) {
    rng.SetSeed(seed);
    Random::SetEngine(rng);
  }

  // ---- gas -------------------------------------------------------------------
  MediumMagboltz gas;
  gas.LoadGasFile(gasfile);
  // r < 0  ->  use Garfield's built-in a-priori parameterization r(c,p)
  //            (Sahin et al. 2010 for Ar/CH4); otherwise use the explicit value.
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

  // Read the true cell extent BEFORE enabling periodicity: once mirror
  // periodicity is on, GetBoundingBox returns +/- infinity in x and y.
  double cxmin, cymin, czmin, cxmax, cymax, czmax;
  fm.GetBoundingBox(cxmin, cymin, czmin, cxmax, cymax, czmax);

  // Periodicity mode is written by the field builder into <outdir>/periodicity:
  //   "mirror"    -> legacy aligned cell, mirror-symmetric (Neumann) walls
  //   "translate" -> offset-capable cell with translational periodic walls
  // Absent file -> mirror, so old run folders keep working unchanged.
  std::string pmode = "mirror";
  {
    std::ifstream pin(outdir + "/periodicity");
    if (pin) std::getline(pin, pmode);
    while (!pmode.empty() &&
           (pmode.back() == '\n' || pmode.back() == '\r' || pmode.back() == ' '))
      pmode.pop_back();
  }
  if (pmode == "translate") {
    fm.EnablePeriodicityX();        // true translational periodicity
    fm.EnablePeriodicityY();
    std::cout << "Periodicity: translational (misalignment-capable)\n";
  } else {
    fm.EnableMirrorPeriodicityX();  // matches the natural (mirror) side-wall BCs
    fm.EnableMirrorPeriodicityY();
    std::cout << "Periodicity: mirror (aligned cell)\n";
  }
  fm.SetGas(&gas);
  fm.PrintRange();

  const double px = cxmax - cxmin;                 // one pitch in x (cm)
  const double py = cymax - cymin;
  const double dz = czmax - czmin;
  const double zStart = czmax - 0.02 * dz;         // just below the drift plane
  const double zReadout = czmin + 0.01 * dz;       // "reached readout" threshold
  std::cout << "Cell (cm): x[" << cxmin << "," << cxmax << "] y[" << cymin
            << "," << cymax << "] z[" << czmin << "," << czmax << "]\n"
            << "Starting electrons at z = " << zStart
            << ", counting those below z = " << zReadout << "\n";

  Sensor sensor(&fm);
  // The mirror-periodic field is defined everywhere, so pad the active area
  // laterally and let electrons roam across cells without leaving it.
  const double pad = 40.0;
  sensor.SetArea(cxmin - pad * px, cymin - pad * py, czmin,
                 cxmax + pad * px, cymax + pad * py, czmax);

  AvalancheMicroscopic aval(&sensor);
  aval.EnableAvalancheSizeLimit(500000);

  std::ofstream csv(outdir + "/gain_stack_events" + tag + ".csv");
  csv << "event,absolute_gain,effective_gain\n";

  double sG = 0., sG2 = 0.;
  int nOk = 0;
  for (int i = 0; i < nEvents; ++i) {
    const double x0 = cxmin + px * RndmUniform();
    const double y0 = cymin + py * RndmUniform();
    aval.AvalancheElectron(x0, y0, zStart, 0., 0.1, 0., 0., 0.);

    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    int nReadout = 0;
    for (const auto& el : aval.GetElectrons()) {
      if (el.path.back().z < zReadout) ++nReadout;
    }
    csv << i << ',' << ne << ',' << nReadout << '\n';
    sG += nReadout; sG2 += double(nReadout) * nReadout;
    ++nOk;
    // progress file for the parent's live multi-bar display
    std::ofstream pf(outdir + "/.progress" + tag, std::ios::trunc);
    pf << (i + 1) << ' ' << nEvents << '\n';
    pf.close();
  }
  csv.close();

  const double mean = sG / nOk;
  const double err = std::sqrt(std::max(sG2 / nOk - mean * mean, 0.) / nOk);
  std::cout << "\n=================================================\n"
            << "Triple-GEM stack effective gain: " << mean << " +/- " << err << "\n"
            << "Penning r = " << rPenning << ",  events = " << nOk << "\n"
            << "Per-event data: " << outdir << "/gain_stack_events.csv\n"
            << "=================================================\n";
  return 0;
}