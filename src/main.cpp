#include "current.hpp"
#include <Grid/Grid.h>

using namespace Grid;

template <typename Action>
void computePropagator(PropagatorField<Action> &prop5, PropagatorField<Action> &prop4,
                       PropagatorField<Action> &src, Action &action)
{
  ConjugateGradient<LatticeFermion> CG(1.0e-16, 100000);
  SchurRedBlackDiagTwoSolve<LatticeFermion> schur(CG);
  ZeroGuesser<LatticeFermion> zpg;

  for (int s = 0; s < Nd; s++)
  {
    for (int c = 0; c < Nc; c++)
    {
      FermionField<Action> src4(src.Grid());
      PropToFerm<Action>(src4, src, s, c);

      FermionField<Action> src5(prop5.Grid());
      action.ImportPhysicalFermionSource(src4, src5);

      FermionField<Action> result5(prop5.Grid());
      result5 = Zero();
      schur(action, src5, result5, zpg);
      std::cout << GridLogMessage << "spin " << s << " color " << c << " norm2(src5d) "
                << norm2(src5) << " norm2(result5d) " << norm2(result5) << std::endl;
      FermToProp<Action>(prop5, result5, s, c);

      FermionField<Action> result4(prop4.Grid());
      action.ExportPhysicalFermionSolution(result5, result4);
      FermToProp<Action>(prop4, result4, s, c);
    }
  }
}

template <class Action>
void TestRegression(Action &Ddwf, GaugeField<Action> &Umu, GridCartesian *FGrid,
                    GridRedBlackCartesian *FrbGrid, GridCartesian *UGrid,
                    GridRedBlackCartesian *UrbGrid, RealD mass, RealD M5,
                    GridParallelRNG *RNG4, GridParallelRNG *RNG5)
{
  PropagatorField<Action> phys_src(UGrid);
  PropagatorField<Action> seqsrc(FGrid), seqsrcalt(FGrid), diff(FGrid);
  PropagatorField<Action> prop5(FGrid);
  PropagatorField<Action> prop4(UGrid);

  random(*RNG4, phys_src);
  random(*RNG5, prop5);

  auto curr = Current::Vector;
  const int mu_J = 0;
  const int t_J = 0;
  double time;
  ComplexField<Action> ph(UGrid);
  ph = 1.0;

  time = -usecond();
  Ddwf.SeqConservedCurrent(prop5, seqsrc, phys_src, curr, mu_J, t_J, t_J, ph);
  time += usecond();
  std::cout << GridLogMessage << "Grid implementation -- norm2= " << norm2(seqsrc)
            << " -- time= " << time / 1.0e6 << "s" << std::endl;
  time = -usecond();
  seqConservedCurrent(prop5, seqsrcalt, phys_src, curr, mu_J, t_J, t_J, ph, Ddwf);
  time += usecond();
  std::cout << GridLogMessage << "this implementation -- norm2= " << norm2(seqsrcalt)
            << " -- time= " << time / 1.0e6 << "s" << std::endl;
  diff = seqsrcalt - seqsrc;
  std::cout << GridLogMessage << "rel diff norm2= " << norm2(diff) / norm2(seqsrc)
            << std::endl;
}

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads"
            << std::endl;

  const int Ls = 10;
  std::vector<ComplexD> omegas;

  omegas.push_back(std::complex<double>(1.45806438985048, -0));
  omegas.push_back(std::complex<double>(1.18231318389348, -0));
  omegas.push_back(std::complex<double>(0.830951166685955, -0));
  omegas.push_back(std::complex<double>(0.542352409156791, -0));
  omegas.push_back(std::complex<double>(0.341985020453729, -0));
  omegas.push_back(std::complex<double>(0.21137902619029, -0));
  omegas.push_back(std::complex<double>(0.126074299502912, -0));
  omegas.push_back(std::complex<double>(0.0990136651962626, -0));
  omegas.push_back(std::complex<double>(0.0686324988446592, 0.0550658530827402));
  omegas.push_back(std::complex<double>(0.0686324988446592, -0.0550658530827402));

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian *FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  std::vector<int> seeds5({5, 6, 7, 8});
  GridParallelRNG RNG5(FGrid);
  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG RNG4(UGrid);
  std::vector<int> seeds4({1, 2, 3, 4});
  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4, Umu);

  RealD mass = 0.3;
  RealD M5 = 1.0;

  std::cout << GridLogMessage << "-- DWF regression test" << std::endl;
  DomainWallFermionD Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);
  TestRegression<DomainWallFermionD>(Ddwf, Umu, FGrid, FrbGrid, UGrid, UrbGrid, mass, M5,
                                     &RNG4, &RNG5);

  std::cout << GridLogMessage << "-- Moebius DWF regression test" << std::endl;
  ScaledShamirFermionD Dsham(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, 2.0);
  TestRegression<ScaledShamirFermionD>(Dsham, Umu, FGrid, FrbGrid, UGrid, UrbGrid, mass,
                                       M5, &RNG4, &RNG5);

  std::cout << GridLogMessage << "-- z-Moebius DWF regression test" << std::endl;
  ZMobiusFermionD ZDmob(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, omegas, 1.0,
                        0.0);
  TestRegression<ZMobiusFermionD>(ZDmob, Umu, FGrid, FrbGrid, UGrid, UrbGrid, mass, M5,
                                  &RNG4, &RNG5);

  Grid_finalize();

  return EXIT_SUCCESS;
}
