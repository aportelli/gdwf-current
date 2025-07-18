/*
 * current.hpp, part of gdwf-current (https://github.com/aportelli/gdwf-current)
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

#include <Grid/Grid.h>

//-------------------------------------------------------------------------------------
// Macro to multiply 4D field by Wilon current
//-------------------------------------------------------------------------------------
//
// qout = JW_mu * qin
//
// qin: input quark (4D)
// qout: output quark (5D)
// sign: global sign
// csign: current sign (e.g. 1 for vector, -1 for tadpole)
// ph: current momentum phase
//
// cut in time to [tmin, tmax]
//
// qout(x) = - 0.5 * csign * sign * ph(x) * (1-gamma_mu) * U_mu(x) * qin(x+mu)
//           + 0.5 * sign * (1+gamma_mu) * U_mu^dag(x-mu) * ph(x-mu) * qin(x-mu)
//
// NB: ImportGauge multiplies links by -1/2, so multLinkField multiplies by -0.5*U_mu
//
#define JWMultiply(qout, qin, sign)                      \
  jtmp = Cshift(qin, mu, 1);                             \
  action.multLinkField(Utmp, action.Umu, jtmp, mu);      \
  jtmp = csign * (sign) * (Utmp * ph - gmu * Utmp * ph); \
  if (doCut)                                             \
  {                                                      \
    jtmp = where((lcoor >= tmin), jtmp, zz);             \
    jtmp = where((lcoor <= tmax), jtmp, zz);             \
  }                                                      \
  qout = jtmp;                                           \
  jtmp = qin * ph;                                       \
  jtmp = Cshift(jtmp, mu, -1);                           \
  action.multLinkField(Utmp, action.Umu, jtmp, mu + Nd); \
  jtmp = -(sign) * (Utmp + gmu * Utmp);                  \
  if (doCut)                                             \
  {                                                      \
    jtmp = where((lcoor >= tmin + tshift), jtmp, zz);    \
    jtmp = where((lcoor <= tmax + tshift), jtmp, zz);    \
  }                                                      \
  qout = qout + jtmp;

#define Gs(s) ((curr_type == Current::Axial) ? (((s) < Ls / 2) ? -1. : 1.) : 1.)
#define Pp(Q) (0.5 * (Q + g5 * Q))
#define Pm(Q) (0.5 * (Q - g5 * Q))

template <typename FImpl>
using FermionField = typename FImpl::FermionField;

template <typename FImpl>
using PropagatorField = typename FImpl::PropagatorField;

template <typename FImpl>
using GaugeField = typename FImpl::GaugeField;

template <typename FImpl>
using ComplexField = typename FImpl::ComplexField;

//-------------------------------------------------------------------------------------
// Conserved current multiplication
//-------------------------------------------------------------------------------------
// q_out = Gs * JW_mu * Omega * ph * q_in + C * JW_mu * ph * eta5$
// with eta5 = (P+ * phys_src, 0, ..., 0, P- * phys_src)
//
// vector current: Gs = 1, JW_mu is Wilson vector current
// axial current: Gs = diag(-1, ..., -1, 1, ..., 1), JW_mu is Wilson vector current
// tadpole current: Gs = 1, JW_mu is Wilson tadpole current
//
// cut in time to [tmin, tmax]
//
// q_in: input propagator
// q_out: output propagator
// phys_src: physical 4D source
// curr_type: current type (Current::Vector/Axial/Tadpole)
// ph: current momentum phase
//
template <typename FImpl, template <typename> class Action>
void seqConservedCurrent(PropagatorField<FImpl> &q_in, PropagatorField<FImpl> &q_out,
                         PropagatorField<FImpl> &phys_src, Grid::Current curr_type,
                         unsigned int mu, unsigned int tmin, unsigned int tmax,
                         ComplexField<FImpl> &ph,
                         Action<FImpl> &action) // Complex phase factor
{
  using namespace Grid;

  const unsigned int nd = q_in.Grid()->Nd(), nt = q_in.Grid()->GlobalDimensions()[nd - 1];
  const unsigned int Ls = action.Ls, tshift = (mu == nd - 1) ? 1 : 0;
  const int csign = (curr_type == Current::Tadpole) ? -1 : 1;
  const bool doCut = ((tmin != 0) || (tmax != nt - 1));
  auto UGrid = action.GaugeGrid();
  auto FGrid = action.FermionGrid();
  auto mass = action.Mass();

  Gamma::Algebra Gmu[] = {Gamma::Algebra::GammaX, Gamma::Algebra::GammaY,
                          Gamma::Algebra::GammaZ, Gamma::Algebra::GammaT};
  Gamma gmu = Gamma(Gmu[mu]);
  Gamma g5(Gamma::Algebra::Gamma5);

  PropagatorField<FImpl> srcTerm(UGrid);
  std::vector<PropagatorField<FImpl>> phi(Ls, UGrid);

  PropagatorField<FImpl> tmp1(UGrid), jtmp(UGrid), tmp2(UGrid);
  PropagatorField<FImpl> Utmp(UGrid);
  PropagatorField<FImpl> zz(UGrid);
  zz = 0.0;
  LatticeInteger lcoor(UGrid);
  LatticeCoordinate(lcoor, Nd - 1);

  for (int s = 0; s < Ls; ++s)
  {
    ExtractSlice(phi[s], q_in, s, 0);
  }

  // q_out = Gs * JW_mu * Omega * ph * q_in
  auto b = action.bs[0];
  auto c = action.cs[0];
  tmp1 = b * phi[0] + c * Pm(phi[1]) - mass * c * Pp(phi[Ls - 1]);
  JWMultiply(tmp2, tmp1, Gs(0));
  InsertSlice(tmp2, q_out, 0, 0);
  for (int s = 1; s < Ls - 1; ++s)
  {
    b = action.bs[s];
    c = action.cs[s];
    tmp1 = c * Pp(phi[s - 1]) + b * phi[s] + c * Pm(phi[s + 1]);
    JWMultiply(tmp2, tmp1, Gs(s));
    InsertSlice(tmp2, q_out, s, 0);
  }
  b = action.bs[Ls - 1];
  c = action.cs[Ls - 1];
  tmp1 = -mass * c * Pm(phi[0]) + c * Pp(phi[Ls - 2]) + b * phi[Ls - 1];
  JWMultiply(tmp2, tmp1, Gs(Ls - 1));
  InsertSlice(tmp2, q_out, Ls - 1, 0);

  // q_out += C * JW_mu * ph * eta5
  tmp1 = action.cs[0] * Pp(phys_src);
  JWMultiply(tmp2, tmp1, Gs(0));
  ExtractSlice(srcTerm, q_out, 0, 0);
  srcTerm += tmp2;
  InsertSlice(srcTerm, q_out, 0, 0);
  tmp1 = action.cs[Ls - 1] * Pm(phys_src);
  JWMultiply(tmp2, tmp1, Gs(Ls - 1));
  ExtractSlice(srcTerm, q_out, Ls - 1, 0);
  srcTerm += tmp2;
  InsertSlice(srcTerm, q_out, Ls - 1, 0);
}

#undef JWMultiply
#undef Gs
#undef Pp
#undef Pm
