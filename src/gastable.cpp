// gastable.cpp -- OPTIONAL Magboltz gas table (drift velocity, diffusion,
// Townsend/attachment vs E). The microscopic gain simulation does NOT use it;
// it is only needed for drift-line (AvalancheMC / DriftLineRKF) work or as a
// cross-check of transport parameters.
//
// Usage:
//   simgem_gastable <gas options> --out gas/xyz.gas --field-min 100 --field-max 1e6
//                   --field-points 30 --max-energy 400 --collisions 10
// Pressure enters through --pressure-torr: 1 bar = 750.06 Torr (not 760).
#include <iostream>

#include "common.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  simgem::Args a(argc, argv);
  auto gas = simgem::MakeGas(a);
  gas->SetMaxElectronEnergy(a.Num("max-energy", 400.));
  gas->SetFieldGrid(a.Num("field-min", 100.), a.Num("field-max", 1.e6),
                    a.Int("field-points", 30), true);
  gas->GenerateGasTable(a.Int("collisions", 10), true);
  if (!gas->WriteGasFile(a.Str("out"))) simgem::Args::Fail("could not write " + a.Str("out"));
  std::cout << "Wrote " << a.Str("out") << "\n";
  return 0;
}
