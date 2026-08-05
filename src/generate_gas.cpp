// generate_gas.cpp
// -----------------------------------------------------------------------------
// Generates a Magboltz gas table for Ar:CH4 98:2 at 8 bar.
//
// This is where the two most common "why is my gain 1000x too low" bugs live:
//   (1) pressure/density set wrong  -> gain scales enormously with N
//   (2) field grid / max energy too small for the high fields inside the holes
//
// Penning transfer is NOT applied here. The excitation and ionisation collision
// frequencies are written into the .gas file, and Penning is switched on at USE
// time (see gain.cpp) via EnablePenningTransfer, which rescales the Townsend
// coefficient as  alpha' = alpha * (1 + r * sum(f_exc)/sum(f_ion)).
// -----------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <string>
#include <Garfield/MediumMagboltz.hh>

using namespace Garfield;

// Usage: ./generate_gas [ar%] [ch4%] [pressure_bar] [temperature_K] [outfile]
// No args -> Ar/CH4 98/2 at 8 bar, 293.15 K -> ar98_ch4_2_8bar.gas
int main(int argc, char* argv[]) {

  const double arPct  = (argc > 1) ? std::atof(argv[1]) : 98.;
  const double ch4Pct = (argc > 2) ? std::atof(argv[2]) : 2.;
  const double pBar   = (argc > 3) ? std::atof(argv[3]) : 8.;
  const double tempK  = (argc > 4) ? std::atof(argv[4]) : 293.15;
  const std::string out = (argc > 5) ? argv[5] : "ar98_ch4_2_8bar.gas";

  MediumMagboltz gas("ar", arPct, "ch4", ch4Pct);

  // ---- Operating point -------------------------------------------------------
  // Pressure is in Torr.  NOTE: 1 bar != 1 atm (1 bar = 750.06 Torr).
  gas.SetPressure(pBar * 750.06);   // bar -> Torr
  gas.SetTemperature(tempK);

  // ---- Collision-rate table --------------------------------------------------
  // Inside the holes at 8 bar the field is high, so allow high electron energies.
  gas.SetMaxElectronEnergy(400.);   // eV (200 is fine at 1 bar; go higher here)

  // Field grid in V/cm. Nominal hole field ~ dV_GEM / 50um ~ 1.1e5 V/cm for 550 V,
  // with dipole peaks a few times higher. Cover well past that.
  gas.SetFieldGrid(100., 1.e6, 30, true);  // 30 points, logarithmic

  // ---- Generate --------------------------------------------------------------
  // numCollisions is in units of 1e7 collisions per grid point.
  // 10 -> 1e8 collisions/point: slow but smooth Townsend coefficients.
  const int numCollisions = 10;
  gas.GenerateGasTable(numCollisions, true);

  gas.WriteGasFile(out);
  std::cout << "Wrote " << out << " (Ar/CH4 " << arPct << "/" << ch4Pct
            << ", " << pBar << " bar, " << tempK << " K)\n";
  return 0;
}
