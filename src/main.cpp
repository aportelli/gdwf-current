/*
 * main.cpp, part of gdwf-current (https://github.com/aportelli/gdwf-current)
 *
 * Copyright Antonin Portelli (C) 2025
 *
 * gdwf-current is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * gdwf-current is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gdwf-current.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

#include "current.hpp"
#include <Grid/Grid.h>

#define BIG_SEP "===================="
#define SEP "--------------------"

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
      FermionField<Action> src5(prop5.Grid());

      if (src.Grid()->Nd() == Nd)
      {
        FermionField<Action> src4(src.Grid());
        PropToFerm<Action>(src4, src, s, c);

        action.ImportPhysicalFermionSource(src4, src5);
      }
      else
      {
        PropToFerm<Action>(src5, src, s, c);
      }

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
void testRegression(Action &Ddwf, GaugeField<Action> &Umu, GridCartesian *FGrid,
                    GridRedBlackCartesian *FrbGrid, GridCartesian *UGrid,
                    GridRedBlackCartesian *UrbGrid, GridParallelRNG *RNG4,
                    GridParallelRNG *RNG5)
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

template <class Action>
void testDerivative(Action &Ddwf, GaugeField<Action> &Umu, GridCartesian *FGrid,
                    GridRedBlackCartesian *FrbGrid, GridCartesian *UGrid,
                    GridRedBlackCartesian *UrbGrid, GridParallelRNG *RNG4,
                    GridParallelRNG *RNG5)
{
  const int mu_J = 0;
  int L_mu = UGrid->GlobalDimensions()[mu_J];
  int nt = UGrid->GlobalDimensions()[Nd - 1];
  const double eps = 1.0e-3;
  std::vector<double> twist(Nd, 0.);
  typename Action::ImplParams implParams;
  PropagatorField<Action> phys_src(UGrid), prop5(FGrid), prop4(UGrid), prop5t(FGrid),
      prop4t(UGrid), seqsrc(FGrid), seqprop5(FGrid), seqprop4(UGrid);
  ComplexField<Action> ph(UGrid);
  auto curr = Current::Vector;

  ph = std::complex(0.0, -1.0);
  random(*RNG4, phys_src);

  std::cout << GridLogMessage << SEP << " untwisted propagator" << std::endl;
  computePropagator(prop5, prop4, phys_src, Ddwf);
  std::cout << GridLogMessage << SEP << " sequential propagator" << std::endl;
  seqConservedCurrent(prop5, seqsrc, phys_src, curr, mu_J, 0, nt - 1, ph, Ddwf);
  computePropagator(seqprop5, seqprop4, seqsrc, Ddwf);
  std::cout << GridLogMessage << SEP << " twisted propagator" << std::endl;
  implParams.twist_n_2pi_L[mu_J] = eps;
  Ddwf.Params = implParams;
  Ddwf.ImportGauge(Umu);
  computePropagator(prop5t, prop4t, phys_src, Ddwf);
  implParams.twist_n_2pi_L[mu_J] = 0.;
  Ddwf.Params = implParams;
  Ddwf.ImportGauge(Umu);
  prop5 = (1. / (eps * (2. * M_PI) / L_mu)) * (prop5t - prop5);
  std::cout << GridLogMessage << SEP << " results" << std::endl;
  std::cout << GridLogMessage << "Numerical derivative -- norm2= " << norm2(prop5)
            << " epsilon= " << eps * (2. * M_PI) / L_mu << std::endl;
  std::cout << GridLogMessage << "Sequential insertion -- norm2= " << norm2(seqprop5)
            << std::endl;
  std::cout << GridLogMessage << "Diff -- norm2_abs= " << norm2(prop5 - seqprop5)
            << " norm2_rel= " << 2.0 * norm2(prop5 - seqprop5) / norm2(prop5 + seqprop5)
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

  std::cout << GridLogMessage << BIG_SEP << " DWF regression test" << std::endl;
  DomainWallFermionD Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);
  testRegression(Ddwf, Umu, FGrid, FrbGrid, UGrid, UrbGrid, &RNG4, &RNG5);

  std::cout << GridLogMessage << BIG_SEP << " Moebius DWF regression test" << std::endl;
  ScaledShamirFermionD Dsham(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, 2.0);
  testRegression(Dsham, Umu, FGrid, FrbGrid, UGrid, UrbGrid, &RNG4, &RNG5);

  std::cout << GridLogMessage << BIG_SEP << " z-Moebius DWF regression test" << std::endl;
  ZMobiusFermionD ZDmob(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, omegas, 1.0,
                        0.0);
  testRegression(ZDmob, Umu, FGrid, FrbGrid, UGrid, UrbGrid, &RNG4, &RNG5);

  std::cout << GridLogMessage << BIG_SEP << " DWF perturbative test" << std::endl;
  testDerivative(Ddwf, Umu, FGrid, FrbGrid, UGrid, UrbGrid, &RNG4, &RNG5);
  std::cout << GridLogMessage << BIG_SEP << " Moebius DWF perturbative test" << std::endl;
  testDerivative(Dsham, Umu, FGrid, FrbGrid, UGrid, UrbGrid, &RNG4, &RNG5);
  std::cout << GridLogMessage << BIG_SEP << " z-Moebius DWF perturbative test"
            << std::endl;
  testDerivative(ZDmob, Umu, FGrid, FrbGrid, UGrid, UrbGrid, &RNG4, &RNG5);

  Grid_finalize();

  return EXIT_SUCCESS;
}
