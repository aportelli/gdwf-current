#include <Grid/Grid.h>

#define JWMultiply(qout, qin, sign)                      \
  jtmp = Cshift(qin, mu, 1);                             \
  action.multLinkField(Utmp, action.Umu, jtmp, mu);      \
  jtmp = (sign) * (Utmp * ph - gmu * Utmp * ph);         \
  jtmp = where((lcoor >= tmin), jtmp, zz);               \
  jtmp = where((lcoor <= tmax), jtmp, zz);               \
  qout = jtmp;                                           \
  jtmp = qin * ph;                                       \
  jtmp = Cshift(jtmp, mu, -1);                           \
  action.multLinkField(Utmp, action.Umu, jtmp, mu + Nd); \
  jtmp = -(sign) * (Utmp + gmu * Utmp);                  \
  jtmp = where((lcoor >= tmin + tshift), jtmp, zz);      \
  jtmp = where((lcoor <= tmax + tshift), jtmp, zz);      \
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

template <typename FImpl, template <typename> class Action>
void seqConservedCurrent(PropagatorField<FImpl> &q_in, PropagatorField<FImpl> &q_out,
                         PropagatorField<FImpl> &phys_src, Grid::Current curr_type,
                         unsigned int mu, unsigned int tmin, unsigned int tmax,
                         ComplexField<FImpl> &ph,
                         Action<FImpl> &action) // Complex phase factor
{
  using namespace Grid;

  int tshift = (mu == Grid::Nd - 1) ? 1 : 0;
  int Ls = action.Ls;
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
