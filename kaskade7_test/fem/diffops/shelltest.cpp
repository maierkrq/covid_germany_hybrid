#include <iostream>
#include <random>

#include "fem/diffops/materialLaws.hh"
#include "fem/diffops/shellKinematics.hh"

using namespace Kaskade;
using namespace Kaskade::Elastomechanics;

int main(void)
{
  using Real = double;

  Dune::FieldMatrix<Real,3,2> dphi_dx;
  Tensor<Real,3,2,2> ddphi_ddx;
  Real t = 0;
  Dune::FieldVector<Real,2>  dt_dx;

  constexpr int dim = 2;

  double h = 1e-8;
  VariationalArg<Real,dim,1> v, w;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1,1);
  v.value[0] = h*dis(gen);
  w.value[0] = h*dis(gen);
  for (int i=0; i<dim; ++i)
  {
    v.derivative[0][i] = h*dis(gen);
    w.derivative[0][i] = h*dis(gen);
    for (int j=0; j<dim; ++j)
    {
      v.hessian[0][i][j] = h*dis(gen);
      w.hessian[0][i][j] = h*dis(gen);
    }
  }



  t = dis(gen)+1.001;
  for (int j=0; j<2; ++j)
    dt_dx[j] = dis(gen);

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<2; ++j)
    {
      dphi_dx[i][j] = dis(gen);
      for (int k=0; k<2; ++k)
        ddphi_ddx[i][j][k] = dis(gen);
    }
  }

  ShellKinematicsDerivative<dim> dPhi(dphi_dx,ddphi_ddx,t,dt_dx);
  MaterialLaws::StVenantKirchhoff<3,double> energyW(ElasticModulus::material("steel"));
  ShellEnergy energy(energyW,1e-3,
                     dphi_dx,ddphi_ddx,t,dt_dx);

  double e = energy.d0();
  std::cout << "stored marginal energy density: " << e << "\n";

  std::cout << "=================================================\nChecking first derivative\n";

  std::cout << "row = 0\n";
  auto d1 = energy.d1<0>(v);
  std::cout << "d1 = " << d1 << "\ndiff ";
  auto d1H = d1;
  for (int i=0; i<dim+1; ++i)
  {
    auto dphi_dxH = dphi_dx; dphi_dxH[i] += v.derivative[0];
    auto ddphi_ddxH = ddphi_ddx; ddphi_ddxH[i] += v.hessian[0];
    ShellEnergy energyH(energyW,1e-3,dphi_dxH,ddphi_ddxH,t,dt_dx);
    d1H[i] = energyH.d0()-e;
  }
  std::cout << d1H << "\n";
  std::cout << "rel error: " << (d1-d1H).two_norm() / (d1.two_norm()+d1H.two_norm()) << "\n\n";

  std::cout << "row = 1\n";
  auto d1t = energy.d1<1>(v);
  std::cout << "d1 = " << d1t << "\ndiff ";
  ShellEnergy energyH(energyW,1e-3,dphi_dx,ddphi_ddx,t+v.value[0],dt_dx+v.derivative[0]);
  auto d1Ht = energyH.d0()-e;
  std::cout << d1Ht << "\n";
  std::cout << "rel error: " << (d1t-d1Ht).two_norm() / (d1t.two_norm()+std::abs(d1Ht)) << "\n\n";

  std::cout << "=================================================\nChecking second derivative\n";

  std::cout << "row = 0, col = 0\n";
  {
    auto d2 = energy.d2<0,0>(v,w);
    auto d2H = d2;
    for (int i=0; i<dim+1; ++i)
    {
      auto dphi_dxH = dphi_dx; dphi_dxH[i] += w.derivative[0];
      auto ddphi_ddxH = ddphi_ddx; ddphi_ddxH[i] += w.hessian[0];
      ShellEnergy energyH(energyW,1e-3,dphi_dxH,ddphi_ddxH,t,dt_dx);
      d2H[i] = energyH.d1<0>(v)-d1;
    }
    d2H = transpose(d2H);
    std::cout << "d2:\n" << d2 << "\ndiff:\n" << d2H << "\n";
    std::cout << "rel error:" << (d2-d2H).frobenius_norm() / (d2.frobenius_norm()+d2H.frobenius_norm()) << "\n\n";
  }

  std::cout << "row = 0, col = 1\n";
  {
    auto d2 = energy.d2<0,1>(v,w);
    double tH = t + w.value[0];
    auto dt_dxH = dt_dx + w.derivative[0];
    ShellEnergy energyH(energyW,1e-3,dphi_dx,ddphi_ddx,tH,dt_dxH);
    auto d2H = energyH.d1<0>(v)-d1;
    std::cout << "d2: \n " << d2 << "\ndiff: " << d2H << "\n";
    Dune::FieldMatrix<double,1,dim+1> d2Hmat; d2Hmat[0] = d2H;
    auto d2num = transpose(d2Hmat);
    std::cout << "rel error:" << (d2-d2num).frobenius_norm() / (d2.frobenius_norm()+d2num.frobenius_norm()) << "\n\n";
  }

  std::cout << "row = 1, col = 0\n";
  {
    auto d2 = energy.d2<1,0>(v,w);
    auto d2H = d2;
    for (int i=0; i<dim+1; ++i)
    {
      auto dphi_dxH = dphi_dx; dphi_dxH[i] += w.derivative[0];
      auto ddphi_ddxH = ddphi_ddx; ddphi_ddxH[i] += w.hessian[0];
      ShellEnergy energyH(energyW,1e-3,dphi_dxH,ddphi_ddxH,t,dt_dx);
      d2H[0][i] = energyH.d1<1>(v)-d1t;
    }
    std::cout << "d2:\n" << d2 << "\ndiff:\n" << d2H << "\n";
    std::cout << "rel error:" << (d2-d2H).frobenius_norm() / (d2.frobenius_norm()+d2H.frobenius_norm()) << "\n\n";
  }

  std::cout << "row = 1, col = 1\n";
  {
    auto d2 = energy.d2<1,1>(v,w);
    double tH = t + w.value[0];
    auto dt_dxH = dt_dx + w.derivative[0];
    ShellEnergy energyH(energyW,1e-3,dphi_dx,ddphi_ddx,tH,dt_dxH);
    auto d2H = energyH.d1<1>(v)-d1t;
    std::cout << "d2: \n " << d2 << "\ndiff: " << d2H << "\n";
    std::cout << "rel error:" << std::abs(d2[0][0]-d2H[0])/ (std::abs(d2[0][0])+std::abs(d2H[0])) << "\n\n";
  }

  return 0;
}
