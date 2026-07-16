/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2020 Zuse Institute Berlin                            */
/*  Copyright (C) 2013-2013 U Milan                                          */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include "dune/grid/config.h"

#include "fem/fixdune.hh"
#include "fem/diffops/membraneModels.hh"

#include "utilities/detailed_exception.hh"

namespace Kaskade 
{
    AlievPanfilov::AlievPanfilov(): MembraneModelBase<AlievPanfilov,1>("Aliev-Panfilov"),
      a(0.1), eps1(0.01), ga(8.0), gs(8.0) , mu1(0.07), mu2(0.3)
    {
    }
 
    double AlievPanfilov::current(double u, Gating const& v_) const
    {
      double v = v_[0];
      return -(ga*u*(u-a)*(u-1.0) + u*v);
    }
    
    double AlievPanfilov::current_du(double u, Gating const& v_) const
    {
      double v = v_[0];
      return std::min(-ga*( (2*u-a)*(u-1.0) + u*(u-a) ) - v, 0.0); // why the min here?
    }
    
    AlievPanfilov::Gating AlievPanfilov::current_dv(double u, Gating const& ) const
    {
      return -u;
    }
    
    AlievPanfilov::Gating AlievPanfilov::gatingRhs(double u, Gating const& v_, double /*Istim*/) const
    {
      double v = v_[0];
      return 0.25*(eps1+mu1*v/(u+mu2))*(-v-gs*u*(u-a-1.0));
    }
    
    AlievPanfilov::Gating AlievPanfilov::gatingRhs_du(double u, Gating const& v_) const
    {
      double v = v_[0];
      return 0.25* ( - mu1*v/((u+mu2)*(u+mu2))*(-v-gs*u*(u-a-1)) - ( eps1+mu1*v/(u+mu2))*(gs*(2*u-a-1)) );       
    }
    
    AlievPanfilov::GatingJacobian AlievPanfilov::gatingRhs_dv(double u, Gating const& v_) const
    {
      double v = v_[0];
      return std::min(0.25*( mu1/(u+mu2)*(-v-gs*u*(u-a-1)) - (eps1+mu1*v/(u+mu2)) ), 0.0);  // why the min here?
    }

  //------------------------------------------------------------------------------------------  
  
  namespace {
    
    TenTusscher16::Gating vRest16()
    {
      TenTusscher16::Gating v;
      
      /*
      v[0] =   0.000080;			// Cai
      v[1] =   0.56;				// CaSR
      v[2] =  11.6;					// Nai
      v[3] = 138.3;					// Ki
      v[4] =   0.0;					// M
      v[5] =   0.75;				// H
      v[6] =   0.75;				// J
      v[7] =   0.0;					// Xr1
      v[8] =   1.0;					// Xr2
      v[9] =   0.0;					// Xs
      v[10] =  0.0;					// R
      v[11] =  1.0;					// S
      v[12] =  0.0;					// D
      v[13] =  1.0;					// F
      v[14] =  1.0;					// FCa
      v[15] =  1.0;					// G
      */
      
      v[0] =   2.25557e-05;			// Cai
      v[1] =   0.0408502;			// CaSR
      v[2] =   5.8553;				// Nai
      v[3] = 140.392;				// Ki
      v[4] =   0.00131023;			// M
      v[5] =   0.777927;			// H
      v[6] =   0.777927;			// J
      v[7] =   0.000176461;			// Xr1
      v[8] =   0.484335;			// Xr2
      v[9] =   0.00295553;			// Xs
      v[10] =  1.957e-08;			// R
      v[11] =  0.999998;			// S
      v[12] =  1.90947e-05;			// D
      v[13] =  0.999925;			// F
      v[14] =  1.00804;				// FCa
      v[15] =  1.0;					// G
      
      return v;
    }
  }
  
  TenTusscher16::TenTusscher16(): MembraneModelBase<TenTusscher16,16>("ten Tusscher 16")      // -8.523e+1
  {
    //this->setRestState(-8.62e+1,vRest16());
    this->setRestState(-86.4956,vRest16());
  }
    
  double TenTusscher16::current(double u, Gating const& v) const
  {

    // State variables
    double V = u;
    double  Cai = v[0];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];

    // Constants
    double          Ko =     5.4;
    double         Cao =     2.0;
    double         Nao =   140.0;
    double           R =  8314.472;
    double           F = 96485.3415;
    double           T =   310.0;
    double         Gkr =     0.096;
    double        pKNa =     0.03;
    double         Gks =     0.245;
    double         GK1 =     5.405;
    double         Gto =     0.294;
    double         GNa =    14.838;
    double        GbNa =     0.00029;
    double         KmK =     1.0;
    double        KmNa =    40.0;
    double        knak =     1.362;
    double        GCaL =     0.000175;
    double        GbCa =     0.000592;
    double       knaca =  1000;
    double       KmNai =    87.5;
    double        KmCa =     1.38;
    double        ksat =     0.1;
    double           n =     0.35;
    double        GpCa =     0.825;
    double        KpCa =     0.0005;
    double         GpK =     0.0146;
 
    // Compute
    double        RTONF = (R*T)/F;
    
    double           Ek = RTONF*(log((Ko/Ki)));
    double          Ena = RTONF*(log((Nao/Nai)));
    double          Eks = RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    double          Eca = 0.5*RTONF*(log((Cao/Cai)));
    double          Ak1 = 0.1/(1.+exp(0.06*(V-Ek-200)));
    double          Bk1 = (3.*exp(0.0002*(V-Ek+100))+
                          exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
    double      rec_iK1 = Ak1/(Ak1+Bk1);
    double     rec_iNaK = (1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
    double      rec_ipK = 1./(1.+exp((25-V)/5.98));
    
    double  INa = GNa*sm*sm*sm*sh*sj*(V-Ena);
    double ICaL = GCaL*sd*sf*sfca*4*V*(F*F/(R*T))*
                  (exp(2*V*F/(R*T))*Cai-0.341*Cao)/(exp(2*V*F/(R*T))-1.);
    double  Ito = Gto*sr*ss*(V-Ek);
      
    double   IKr = Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(V-Ek);
    double   IKs = Gks*sxs*sxs*(V-Eks);
    double   IK1 = GK1*rec_iK1*(V-Ek);
    double INaCa = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
                  (1./(1+ksat*exp((n-1)*V*F/(R*T))))*
                  (exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-
                   exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    double INaK = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    double IpCa = GpCa*Cai/(KpCa+Cai);
    double  IpK = GpK*rec_ipK*(V-Ek);
    double IbNa = GbNa*(V-Ena);
    double IbCa = GbCa*(V-Eca);
    
    
    double        dV_dt = -(IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK);

    return dV_dt;
  }


    double TenTusscher16::current_du(double u, Gating const& v, double h) const
    {
      /*
      // numerical
      double fuv  = (*this).current(u-h,v);
      double fuiv = (*this).current(u+h,v);
      return (fuiv-fuv)/(2*h);
      */

 
    // State variables
    double    V = u;
    double  Cai = v[0];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];

    // Constants
    double          Ko =     5.4;
    double         Cao =     2.0;
    double         Nao =   140.0;
    double           R =  8314.472;
    double           F = 96485.3415;
    double           T =   310.0;
    double         Gkr =     0.096;
    double         Gks =     0.245;
    double         GK1 =     5.405;
    double         Gto =     0.294;
    double         GNa =    14.838;
    double        GbNa =     0.00029;
    double         KmK =     1.0;
    double        KmNa =    40.0;
    double        knak =     1.362;
    double        GCaL =     0.000175;
    double        GbCa =     0.000592;
    double       knaca =  1000;
    double       KmNai =    87.5;
    double        KmCa =     1.38;
    double        ksat =     0.1;
    double           n =     0.35;
    double         GpK =     0.0146;

    // Compute
    double        RTONF = (R*T)/F;
    
    double           Ek = RTONF*(log((Ko/Ki)));
    double          Ak1 = 0.1/(1.+exp(0.06*(V-Ek-200)));
    double          Bk1 = (3.*exp(0.0002*(V-Ek+100))+
                          exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
    double      rec_iK1 = Ak1/(Ak1+Bk1);
    double      rec_ipK = 1./(1.+exp((25-V)/5.98));

    double termA  = GCaL*sd*sf*sfca*4*V*(F*F/(R*T));
    double termE1 = exp(2*V*F/(R*T));
    double termE2 = exp(2*V*F/(R*T));
      
    termA         = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    termE1        = ksat*exp((n-1)*V*F/(R*T));
    termE1        = 1./(1+termE1);
    termE2        = exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao;
    double termE3 = exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*2.5;


    double IKr_u   = Gkr*sqrt(Ko/5.4)*sxr1*sxr2;
    double IKs_u   = Gks*sxs*sxs;
    double IK1_u   = GK1*rec_iK1;
    double Ito_u   = Gto*sr*ss;
    double INa_u   = GNa*sm*sm*sm*sh*sj;
    double IbNa_u  = GbNa;
    double IbCa_u  = GbCa;
    
    termA   = GCaL*sd*sf*sfca*4*(F*F/(R*T));
    double termB = 0.341*Cao;
    termE1  = 2*F/(R*T);
    double termQ = (exp(V*termE1)*Cai-termB)/(exp(V*termE1)-1.);
    double termQ_u = (Cai*termE1*exp(termE1*V) - (Cai*exp(termE1*V)-termB)*termE1*exp(termE1*V) ) / ( (exp(V*termE1)-1.)*(exp(V*termE1)-1.) );
    double ICaL_u  = termA*(termQ + V*termQ_u);
      
    termA      = -0.1*V*F/(R*T);
    termB      = -F/(R*T);
    termE1     = exp(termA*V);
    termE2     = exp(termB*V);
    termQ      = 1.+0.1245*termE1+0.0353*termE2;
    double rec_iNaK_u = - (0.1245*termA*termE1 + 0.0353*termB*termE2)/(termQ*termQ);
    double termC  = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa));
    double INaK_u = termC*rec_iNaK_u;


    termA   = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    double termC1  = (n-1)*F/(R*T);
    double termC2  = n*F/(R*T);
    double termC3  = (n-1)*F/(R*T);
    double termF1  = ksat;
    double termF2  = Nai*Nai*Nai*Cao;
    double termF3  = Nao*Nao*Nao*Cai*2.5;
    termE1   = exp(termC1*V);
    termE2   = exp(termC2*V);
    termE3   = exp(termC3*V);
    double termE1_u = termC1*termE1;
    double termE2_u = termC2*termE2;
    double termE3_u = termC3*termE3;
    double INaCa_u = termA * 
              ( termF2*termE2_u*(1+termF1*termE1) - termF1*termE1_u*termF2*termE2 ) / ( (1+termF1*termE1)*(1+termF1*termE1) );
    INaCa_u -= termA * 
                ( termF3*termE3_u*(1+termF1*termE1) - termF1*termE1_u*termF3*termE3 ) / ( (1+termF1*termE1)*(1+termF1*termE1) );

       
    double IpCa_u  = GbCa;
      
    termE1    =  exp((25-V)/5.98);
    termE1_u  = -termE1/5.98;
    double rec_ipK_u = -termE1_u / ((1.+termE1)*(1.+termE1));
    double IpK_u     = GpK*(rec_ipK_u*(V-Ek) + rec_ipK);

      
      return IKr_u + IKs_u + IK1_u + Ito_u + INa_u + IbNa_u + ICaL_u + IbCa_u + INaK_u + IpCa_u + IpK_u + INaCa_u;
    }
    

    TenTusscher16::Gating TenTusscher16::current_dv(double u, Gating const& v) const
    {
      double f = (*this).current(u,v);
      Gating fv_uv, dir;
      for (int i=0; i<nGating; ++i)
      {
        dir = 0;
        dir[i] = 1e-5;
        fv_uv[i] = ((*this).current(u,v+dir)-f)/1e-5;
      }
      return fv_uv;
    }
  
  
  TenTusscher16::Gating TenTusscher16::gatingRhs(double u, Gating const& v, double Istim) const
  {
 
    // State variables
    double    V = u;
    double  Cai = v[0];
    double CaSR = v[1];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];
    double   sg = v[15];

    // Constants
    double          Ko =     5.4;
    double         Cao =     2.0;
    double         Nao =   140.0;
    double          Vc =     0.016404;
    double         Vsr =     0.001094;
    double        Bufc =     0.15;
    double       Kbufc =     0.001;
    double       Bufsr =    10.;
    double      Kbufsr =     0.3;
    double        taug =     2.;
    double      Vmaxup =     0.000425;
    double         Kup =     0.00025;
    double           R =  8314.472;
    double           F = 96485.3415;
    double           T =   310.0;
    double CAPACITANCE =     0.185;
    double         Gkr =     0.096;
    double        pKNa =     0.03;
    double         Gks =     0.245;
    double         GK1 =     5.405;
    double         Gto =     0.294;
    double         GNa =    14.838;
    double        GbNa =     0.00029;
    double         KmK =     1.0;
    double        KmNa =    40.0;
    double        knak =     1.362;
    double        GCaL =     0.000175;
    double        GbCa =     0.000592;
    double       knaca =  1000;
    double       KmNai =    87.5;
    double        KmCa =     1.38;
    double        ksat =     0.1;
    double           n =     0.35;
    double        GpCa =     0.825;
    double        KpCa =     0.0005;
    double         GpK =     0.0146;

    // Compute
    double        RTONF = (R*T)/F;
    double  inverseVcF2 = 1/(2*Vc*F);
    double   inverseVcF = 1./(Vc*F);
    double    Kupsquare = Kup*Kup;
    double    BufcKbufc = Bufc*Kbufc;
    double  BufsrKbufsr = Bufsr*Kbufsr;
    
    double           Ek = RTONF*(log((Ko/Ki)));
    double          Ena = RTONF*(log((Nao/Nai)));
    double          Eks = RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    double          Eca = 0.5*RTONF*(log((Cao/Cai)));
    double          Ak1 = 0.1/(1.+exp(0.06*(V-Ek-200)));
    double          Bk1 = (3.*exp(0.0002*(V-Ek+100))+
                          exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
    double      rec_iK1 = Ak1/(Ak1+Bk1);
    double     rec_iNaK = (1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
    double      rec_ipK = 1./(1.+exp((25-V)/5.98));
    
    double  INa = GNa*sm*sm*sm*sh*sj*(V-Ena);
    double ICaL = GCaL*sd*sf*sfca*4*V*(F*F/(R*T))*
                  (exp(2*V*F/(R*T))*Cai-0.341*Cao)/(exp(2*V*F/(R*T))-1.);
    double  Ito = Gto*sr*ss*(V-Ek);
      
    double   IKr = Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(V-Ek);
    double   IKs = Gks*sxs*sxs*(V-Eks);
    double   IK1 = GK1*rec_iK1*(V-Ek);
    double INaCa = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
                  (1./(1+ksat*exp((n-1)*V*F/(R*T))))*
                  (exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-
                   exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    double INaK = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    double IpCa = GpCa*Cai/(KpCa+Cai);
    double  IpK = GpK*rec_ipK*(V-Ek);
    double IbNa = GbNa*(V-Ena);
    double IbCa = GbCa*(V-Eca);

    double termA  = GCaL*sd*sf*sfca*4*V*(F*F/(R*T));
    double termE1 = exp(2*V*F/(R*T));
    double termE2 = exp(2*V*F/(R*T));
      
    termA         = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    termE1        = ksat*exp((n-1)*V*F/(R*T));
    termE1        = 1./(1+termE1);
    termE2        = exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao;
    double termE3 = exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*2.5;

 
    termA   = GCaL*sd*sf*sfca*4*(F*F/(R*T));
    double termB = 0.341*Cao;
    termE1  = 2*F/(R*T);
                          
    termA      = -0.1*V*F/(R*T);
    termB      = -F/(R*T);
    termE1     = exp(termA*V);
    termE2     = exp(termB*V);

    termA   = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    double termC1  = (n-1)*F/(R*T);
    double termC2  = n*F/(R*T);
    double termC3  = (n-1)*F/(R*T);
    double termF1  = ksat;
    double termF2  = Nai*Nai*Nai*Cao;
    double termF3  = Nao*Nao*Nao*Cai*2.5;
    termE1   = exp(termC1*V);
    termE2   = exp(termC2*V);
    termE3   = exp(termC3*V);
    double termE1_u = termC1*termE1;
    double termE2_u = termC2*termE2;
    double termE3_u = termC3*termE3;
    double INaCa_u = termA * 
              ( termF2*termE2_u*(1+termF1*termE1) - termF1*termE1_u*termF2*termE2 ) / ( (1+termF1*termE1)*(1+termF1*termE1) );
    INaCa_u -= termA * 
                ( termF3*termE3_u*(1+termF1*termE1) - termF1*termE1_u*termF3*termE3 ) / ( (1+termF1*termE1)*(1+termF1*termE1) );

       
      
    termE1    =  exp((25-V)/5.98);
    termE1_u  = -termE1/5.98;

    double Caisquare  = Cai*Cai;
    double CaSRsquare = CaSR*CaSR;
    double CaCurrent=-(ICaL+IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE;
    double A = 0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
    double Irel  = A*sd*sg;
    double Ileak = 0.00008*(CaSR-Cai);
    double SERCA = Vmaxup/(1.+(Kupsquare/Caisquare));
    double CaSRCurrent = SERCA-Irel-Ileak;
    double CaiBufc = 1.0/(1.0 + (BufcKbufc)/( (Cai + Kbufc)*(Cai + Kbufc) ));

    BufsrKbufsr=Bufsr*Kbufsr;
    double CaSRbufSR = 1.0/(1.0+BufsrKbufsr/((CaSR+Kbufsr)*(CaSR+Kbufsr)));
    
    double AM = 1./(1.+exp((-60.-V)/5.));
    double BM = .1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
    double TAU_M = AM*BM;
    double M_INF = 1./((1.+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));
    
    double TAU_H = 0.0;
    if (V>=-40.)
    {
      double AH_1=0.; 
      double BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
      TAU_H= 1.0/(AH_1+BH_1);
    }
    else
    {
      double AH_2=(0.057*exp(-(V+80.)/6.8));
      double BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
      TAU_H=1.0/(AH_2+BH_2);
    }
    double H_INF=1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));

    double TAU_J = 0.0;
    if(V >= -40.)
    {
      double AJ_1=0.;      
      double BJ_1=(0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
      TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
      double AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
                     exp(-0.04391*V))*(V+37.78)/
                         (1.+exp(0.311*(V+79.23))));    
      double BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
      TAU_J = 1.0/(AJ_2+BJ_2);
    }
    double J_INF = 1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));
    
    double axr1    = 450./(1.+exp((-45.-V)/10.));
    double bxr1    = 6./(1.+exp((V-(-30.))/11.5));
    double TAU_Xr1 = axr1*bxr1;
    double Xr1_INF = 1./(1.+exp((-26.-V)/7.));

    double axr2    = 3./(1.+exp((-60.-V)/20.));
    double bxr2    = 1.12/(1.+exp((V-60.)/20.));
    double TAU_Xr2 = axr2*bxr2;
    double Xr2_INF = 1./(1.+exp((V-(-88.))/24.));

    double Axs = 1100./(sqrt(1.+exp((-10.-V)/6)));
    double Bxs = 1./(1.+exp((V-60.)/20.));
    double TAU_Xs = Axs*Bxs;
    double Xs_INF = 1./(1.+exp((-5.-V)/14.));

    double TAU_R = 9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;
    double R_INF = 1./(1.+exp((20-V)/6.));

    double TAU_S = 85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;
    double S_INF=1./(1.+exp((V+20)/5.));

    double Ad = 1.4/(1.+exp((-35-V)/13))+0.25;
    double Bd = 1.4/(1.+exp((V+5)/5));
    double Cd = 1./(1.+exp((50-V)/20));
    double TAU_D = Ad*Bd+Cd;
    double D_INF = 1./(1.+exp((-5-V)/7.5));

    double F_INF = 1./(1.+exp((V+20.)/7.));
    double TAU_F = 1125*exp(-(V+27)*(V+27)/300.)+80.0+165.0/(1.+exp((25-V)/10.));

    double FCa_INF = (1./(1.+pow((Cai/0.000325),8))+
                      0.1/(1.+exp((Cai-0.0005)/0.0001))+
                      0.20/(1.+exp((Cai-0.00075)/0.0008))+ 0.23 )/1.46;
    double tau_fCa = 2.0; 
    double d_fCa = (FCa_INF-sfca)/tau_fCa;
    if ( (FCa_INF > sfca) && (V > -37) ) d_fCa = 0.0;    
    
    double G_INF = 0.0;
    if (Cai<0.00035)
       G_INF=1./(1.+pow((Cai/0.00035),6));
    else
       G_INF=1./(1.+pow((Cai/0.00035),16));
    double tau_g = taug;                     
    double d_g = (G_INF - sg)/tau_g;
    if ( (G_INF > sg) && (V > -37) ) d_g = 0.0;
    
    
    double  Cai_dt = CaiBufc * (-CaSRCurrent + CaCurrent);
    double CaSR_dt = CaSRbufSR * Vc/Vsr * CaSRCurrent;
    double  Nai_dt = -(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
    
    double   Ki_dt = -(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
    
    double   sm_dt = (M_INF-sm)/TAU_M;; 
    double   sh_dt = (H_INF-sh)/TAU_H;
    double   sj_dt = (J_INF-sj)/TAU_J;
    double sxr1_dt = (Xr1_INF-sxr1)/TAU_Xr1;
    double sxr2_dt = (Xr2_INF-sxr2)/TAU_Xr2;
    double  sxs_dt = (Xs_INF-sxs)/TAU_Xs;
    double   sr_dt = (R_INF-sr)/TAU_R;
    double   ss_dt = (S_INF-ss)/TAU_S;
    double   sd_dt = (D_INF-sd)/TAU_D;
    double   sf_dt = (F_INF-sf)/TAU_F;
    double sfca_dt = d_fCa;
    double   sg_dt = d_g;

    

    // RHS
    Gating dv_dt;
    
    dv_dt[0]  = Cai_dt;
    dv_dt[1]  = CaSR_dt;
    dv_dt[2]  = Nai_dt;
    dv_dt[3]  = Ki_dt;
    dv_dt[4]  = sm_dt;
    dv_dt[5]  = sh_dt;
    dv_dt[6]  = sj_dt;
    dv_dt[7]  = sxr1_dt;
    dv_dt[8]  = sxr2_dt;
    dv_dt[9]  = sxs_dt;
    dv_dt[10] = sr_dt;
    dv_dt[11] = ss_dt;
    dv_dt[12] = sd_dt;
    dv_dt[13] = sf_dt;
    dv_dt[14] = sfca_dt;
    dv_dt[15] = sg_dt;
    
    for (int i=0; i<nGating; ++i)
      if (std::isnan(dv_dt[i]))
//         throw DetailedException("NaN occured in ten Tusscher computation",__FILE__,__LINE__);
        std::cerr << "dv_dt[" << i << "] = NaN\n";

    return dv_dt;
  }
  

  TenTusscher16::Gating TenTusscher16::gatingRhs_du(double u, Gating const& v) const
    {
      Gating guA = (*this).gatingRhs(u-1e-5,v);
      Gating guB = (*this).gatingRhs(u+1e-5,v);
      return (guB-guA)/2e-5;
    }
    

  TenTusscher16::GatingJacobian TenTusscher16::gatingRhs_dv(double u, Gating const& v) const
    {
      GatingJacobian gv_uv;
      //std::cout << "TenTusscher16::GatingJacobian TenTusscher16::gatingRhs_dv" << std::endl;
      
      /*
      Gating col, dir, guv = (*this).gatingRhs(u,v);
      for (int i=0; i<nGating; ++i)
      {
        dir = 0;
        dir[i]=1e-8;
        col = ((*this).gatingRhs(u,v+dir)-guv)/1e-8;
        for (int j=0; j<nGating; ++j) gv_uv[j][i] = col[j];
      }
      */
     
    // analytical evaluation in diagonal 
    // State variables
    double    V = u;
    double  Cai = v[0];
    double CaSR = v[1];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];
    double   sg = v[15];

    // Constants
    double          Ko =     5.4;
    double         Cao =     2.0;
    double         Nao =   140.0;
    double          Vc =     0.016404;
    double         Vsr =     0.001094;
    double        Bufc =     0.15;
    double       Kbufc =     0.001;
    double       Bufsr =    10.;
    double      Kbufsr =     0.3;
    double        taug =     2.;
    double      Vmaxup =     0.000425;
    double         Kup =     0.00025;
    double           R =  8314.472;
    double           F = 96485.3415;
    double           T =   310.0;
    double CAPACITANCE =     0.185;
    double         Gkr =     0.096;
    double        pKNa =     0.03;
    double         Gks =     0.245;
    double         GK1 =     5.405;
    double         Gto =     0.294;
    double         GNa =    14.838;
    double        GbNa =     0.00029;
    double         KmK =     1.0;
    double        KmNa =    40.0;
    double        knak =     1.362;
    double        GCaL =     0.000175;
    double        GbCa =     0.000592;
    double       knaca =  1000;
    double       KmNai =    87.5;
    double        KmCa =     1.38;
    double        ksat =     0.1;
    double           n =     0.35;
    double        GpCa =     0.825;
    double        KpCa =     0.0005;
    double         GpK =     0.0146;

    // Compute
    double        RTONF = (R*T)/F;
    double  inverseVcF2 = 1/(2*Vc*F);
    double   inverseVcF = 1./(Vc*F);
    double    Kupsquare = Kup*Kup;
    double    BufcKbufc = Bufc*Kbufc;
    double  BufsrKbufsr = Bufsr*Kbufsr;
    
    double           Ek = RTONF*(log((Ko/Ki)));
    double          Eca = 0.5*RTONF*(log((Cao/Cai)));
    double          Ak1 = 0.1/(1.+exp(0.06*(V-Ek-200)));
    double          Bk1 = (3.*exp(0.0002*(V-Ek+100))+
                          exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
    double      rec_iK1 = Ak1/(Ak1+Bk1);
    double     rec_iNaK = (1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
    double      rec_ipK = 1./(1.+exp((25-V)/5.98));
    
    double ICaL = GCaL*sd*sf*sfca*4*V*(F*F/(R*T))*
                  (exp(2*V*F/(R*T))*Cai-0.341*Cao)/(exp(2*V*F/(R*T))-1.);
      
    double INaCa = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
                  (1./(1+ksat*exp((n-1)*V*F/(R*T))))*
                  (exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-
                   exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    double IpCa = GpCa*Cai/(KpCa+Cai);
    double IbCa = GbCa*(V-Eca);

    double Ek_Ki      = -RTONF/Ki; 
    double Eks_Ki     = -RTONF/(Ki+pKNa*Nai); 
    double Ak1_Ki = -0.1*0.06*Ek_Ki *exp(0.06*(V-Ek-200))/( (1.+exp(0.06*(V-Ek-200))) * (1.+exp(0.06*(V-Ek-200))) );
    double vv = 3.*exp(0.0002*(V-Ek+100));
    double ww = exp(0.1*(V-Ek-10));
    double zz = 1.+exp(-0.5*(V-Ek));
    double vv_Ki = vv*(-0.0002*Ek_Ki);
    double ww_Ki = ww*(-0.1*Ek_Ki);
    double zz_K1 = 0.5*Ek_Ki*(zz-1);
    double Bk1_Ki = (zz*(vv_Ki+ww_Ki) - (vv+ww)*zz_K1) / (zz*zz);
    double rec_iK1_Ki = ( (Ak1+Bk1)*Ak1_Ki - Ak1*(Ak1_Ki+Bk1_Ki) ) /( (Ak1+Bk1)*(Ak1+Bk1) );
    double IK1_Ki  =  GK1* ( rec_iK1_Ki*(V-Ek) - rec_iK1*Ek_Ki );
    double Ito_Ki  = -Gto*sr*ss*Ek_Ki;
    double IKr_Ki  = -Gkr*sqrt(Ko/5.4)*sxr1*sxr2*Ek_Ki;
    double IKs_Ki  = -Gks*sxs*sxs*Eks_Ki;
    double INaK_Ki =  0.0;
    double IpK_Ki  = -GpK*rec_ipK*Ek_Ki;
      
    double Ena_nai  = -RTONF/Nai; 
    double INa_nai   = -GNa*sm*sm*sm*sh*sj*Ena_nai;
    double IbNa_nai  = -GbNa*Ena_nai;
    double INaK_nai  =  knak*(Ko/(Ko+KmK))*rec_iNaK* KmNa/( (Nai+KmNa)*(Nai+KmNa) );
    double fac1 = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    double fac2 = (1./(1+ksat*exp((n-1)*V*F/(R*T))));
    double fac3 = exp(n*V*F/(R*T));
    double INaCa_nai =  3*fac1*fac2*fac3*Nai*Nai*Cao;

    double Eca_cai  = -0.5*RTONF/Cai;
    double IbCa_cai = -GbCa*Eca_cai;
    double IpCa_cai =  GpCa*KpCa/( (KpCa+Cai)*(KpCa+Cai) );
      
    double termA  = GCaL*sd*sf*sfca*4*V*(F*F/(R*T));
    double termE1 = exp(2*V*F/(R*T));
    double termE2 = exp(2*V*F/(R*T));
    double ICaL_cai      = termA*termE1/(termE2-1.);
      
    termA         = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    termE1        = ksat*exp((n-1)*V*F/(R*T));
    termE1        = 1./(1+termE1);
    termE2        = exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao;
    double termE3 = exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*2.5;
    double INaCa_cai = -termA*termE1*termE3;

    
    termA   = GCaL*sd*sf*sfca*4*(F*F/(R*T));
    double termB = 0.341*Cao;
    termE1  = 2*F/(R*T);
      
    termA      = -0.1*V*F/(R*T);
    termB      = -F/(R*T);
    termE1     = exp(termA*V);
    termE2     = exp(termB*V);

    termA   = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    double termC1  = (n-1)*F/(R*T);
    double termC2  = n*F/(R*T);
    double termC3  = (n-1)*F/(R*T);
    double termF1  = ksat;
    double termF2  = Nai*Nai*Nai*Cao;
    double termF3  = Nao*Nao*Nao*Cai*2.5;
    termE1   = exp(termC1*V);
    termE2   = exp(termC2*V);
    termE3   = exp(termC3*V);
    double termE1_u = termC1*termE1;
    double termE2_u = termC2*termE2;
    double termE3_u = termC3*termE3;
    double INaCa_u = termA * 
              ( termF2*termE2_u*(1+termF1*termE1) - termF1*termE1_u*termF2*termE2 ) / ( (1+termF1*termE1)*(1+termF1*termE1) );
    INaCa_u -= termA * 
                ( termF3*termE3_u*(1+termF1*termE1) - termF1*termE1_u*termF3*termE3 ) / ( (1+termF1*termE1)*(1+termF1*termE1) );

       
      
    termE1    =  exp((25-V)/5.98);
    termE1_u  = -termE1/5.98;

    double Caisquare  = Cai*Cai;
    double CaSRsquare = CaSR*CaSR;
    double CaCurrent=-(ICaL+IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE;
    double A = 0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
    double Irel  = A*sd*sg;
    double Ileak = 0.00008*(CaSR-Cai);
    double SERCA = Vmaxup/(1.+(Kupsquare/Caisquare));
    double CaSRCurrent = SERCA-Irel-Ileak;
    double CaiBufc = 1.0/(1.0 + (BufcKbufc)/( (Cai + Kbufc)*(Cai + Kbufc) ));

    BufsrKbufsr=Bufsr*Kbufsr;
    double CaSRbufSR = 1.0/(1.0+BufsrKbufsr/((CaSR+Kbufsr)*(CaSR+Kbufsr)));
    
    double AM = 1./(1.+exp((-60.-V)/5.));
    double BM = .1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
    double TAU_M = AM*BM;
    
    double TAU_H = 0.0;
    if (V>=-40.)
    {
      double AH_1=0.; 
      double BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
      TAU_H= 1.0/(AH_1+BH_1);
    }
    else
    {
      double AH_2=(0.057*exp(-(V+80.)/6.8));
      double BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
      TAU_H=1.0/(AH_2+BH_2);
    }

    double TAU_J = 0.0;
    if(V >= -40.)
    {
      double AJ_1=0.;      
      double BJ_1=(0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
      TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
      double AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
                     exp(-0.04391*V))*(V+37.78)/
                         (1.+exp(0.311*(V+79.23))));    
      double BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
      TAU_J = 1.0/(AJ_2+BJ_2);
    }
    double axr1    = 450./(1.+exp((-45.-V)/10.));
    double bxr1    = 6./(1.+exp((V-(-30.))/11.5));
    double TAU_Xr1 = axr1*bxr1;

    double axr2    = 3./(1.+exp((-60.-V)/20.));
    double bxr2    = 1.12/(1.+exp((V-60.)/20.));
    double TAU_Xr2 = axr2*bxr2;

    double Axs = 1100./(sqrt(1.+exp((-10.-V)/6)));
    double Bxs = 1./(1.+exp((V-60.)/20.));
    double TAU_Xs = Axs*Bxs;

    double TAU_R = 9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;
    double TAU_S = 85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;

    double Ad = 1.4/(1.+exp((-35-V)/13))+0.25;
    double Bd = 1.4/(1.+exp((V+5)/5));
    double Cd = 1./(1.+exp((50-V)/20));
    double TAU_D = Ad*Bd+Cd;

    double TAU_F = 1125*exp(-(V+27)*(V+27)/300.)+80.0+165.0/(1.+exp((25-V)/10.));

    double FCa_INF = (1./(1.+pow((Cai/0.000325),8))+
                      0.1/(1.+exp((Cai-0.0005)/0.0001))+
                      0.20/(1.+exp((Cai-0.00075)/0.0008))+ 0.23 )/1.46;
    double tau_fCa = 2.0; 
    
    double G_INF = 0.0;
    if (Cai<0.00035)
       G_INF=1./(1.+pow((Cai/0.00035),6));
    else
       G_INF=1./(1.+pow((Cai/0.00035),16));
    double tau_g = taug;                     
      
      double c1   = Kbufc;
      double c2   = BufcKbufc;
      double facl = (Cai+c1)*(Cai+c1) +c2;
      double dCaiBufc     = 2.0*c2*(Cai+c1) / (facl*facl);;

      c1   = Vmaxup;
      c2   = Kupsquare;
      facl = Caisquare + c2;
      double dSERCA = 2*c1*c2*Cai /(facl*facl);
      double dIrel = 0.0;
      double dIleak = -0.00008*Cai;
      double dCaSRCurrent = dSERCA-dIrel-dIleak;
      double dCaCurrent   = -(ICaL_cai+IbCa_cai+IpCa_cai-2*INaCa_cai)*inverseVcF2*CAPACITANCE;;

      gv_uv[0][0] = dCaiBufc*(-CaSRCurrent+CaCurrent) + CaiBufc*(-dCaSRCurrent+dCaCurrent);
 
 
      dIleak = 0.00008;
      double dA = -0.016464*2*0.0625*CaSR/( (0.0625+CaSRsquare)*0.0625+CaSRsquare );
      dIrel = dA*sd*sg;
      dCaSRCurrent = -dIrel - dIleak;
      c1 = Kbufsr;
      c2 = BufsrKbufsr;
      facl = (CaSR+c1)*(CaSR+c1) +c2;
      double dCaSRbufSR   = 2.0*c2*(CaSR+c1) / (facl*facl);
      
      gv_uv[1][1] = Vc/Vsr * (CaSRbufSR * dCaSRCurrent + dCaSRbufSR * CaSRCurrent);  


      gv_uv[2][2] = -(INa_nai+IbNa_nai+3*INaK_nai+3*INaCa_nai)*inverseVcF*CAPACITANCE;;  


      gv_uv[3][3] = -(IK1_Ki+Ito_Ki+IKr_Ki+IKs_Ki-2*INaK_Ki+IpK_Ki)*inverseVcF*CAPACITANCE;

      AM    = 1./(1.+exp((-60.-V)/5.));
      BM    = 0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
      TAU_M = AM*BM;
      gv_uv[4][4] = -1.0/TAU_M; 


      TAU_H = 0.0;
      if (V>=-40.)
	  {
        double AH_1=0.; 
        double BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
        TAU_H= 1.0/(AH_1+BH_1);
      }
      else
      {
        double AH_2=(0.057*exp(-(V+80.)/6.8));
        double BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
        TAU_H=1.0/(AH_2+BH_2);
      }
      gv_uv[5][5] = -1.0/TAU_H; 


      TAU_J = 0.0;
      if(V >= -40.)
      {
        double AJ_1=0.;      
        double BJ_1=(0.6*exp((0.057)*u)/(1.+exp(-0.1*(V+32.))));
        TAU_J= 1.0/(AJ_1+BJ_1);
      }
      else
      {
        double AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
                       exp(-0.04391*V))*(V+37.78)/
                        (1.+exp(0.311*(V+79.23))));    
        double BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
        TAU_J = 1.0/(AJ_2+BJ_2);
      }
      gv_uv[6][6] = -1.0/TAU_J; 


      axr1 = 450./(1.+exp((-45.-V)/10.));
      bxr1 = 6./(1.+exp((V-(-30.))/11.5));
      TAU_Xr1 = axr1*bxr1;
      gv_uv[7][7] = -1.0/TAU_Xr1; 


      axr2 = 3./(1.+exp((-60.-V)/20.));
      bxr2 = 1.12/(1.+exp((V-60.)/20.));
      TAU_Xr2 = axr2*bxr2;
      gv_uv[8][8] = -1.0/TAU_Xr2; 


      Axs = 1100./(sqrt(1.+exp((-10.-V)/6)));
      Bxs = 1./(1.+exp((V-60.)/20.));
      TAU_Xs = Axs*Bxs;
      gv_uv[9][9] = -1.0/TAU_Xs; 


      TAU_R = 9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;
      gv_uv[10][10] = -1.0/TAU_R; 


      TAU_S = 85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;
      gv_uv[11][11] = -1.0/TAU_S; 


      Ad = 1.4/(1.+exp((-35-V)/13))+0.25;
      Bd = 1.4/(1.+exp((V+5)/5));
      Cd=1./(1.+exp((50-V)/20));
      TAU_D=Ad*Bd+Cd;
      gv_uv[12][12] = -1.0/TAU_D; 


      TAU_F = 1125*exp(-(V+27)*(V+27)/300.)+80.0+165.0/(1.+exp((25-V)/10.));
      gv_uv[13][13] = -1.0/TAU_F; 


      FCa_INF = (1./(1.+pow((Cai/0.000325),8))+
                       0.1/(1.+exp((Cai-0.0005)/0.0001))+
                       0.20/(1.+exp((Cai-0.00075)/0.0008))+
                       0.23 )/1.46;
      tau_fCa = 2.0; 
      double d_fCa_fca = -1.0/tau_fCa;
      if ( (FCa_INF > sfca) && (V > -37) ) gv_uv[14][14] = 0.0;
      else gv_uv[14][14] = d_fCa_fca; 
 
 
      G_INF = 0.0;
      if (Cai<0.00035)
        G_INF=1./(1.+pow((Cai/0.00035),6));
      else
        G_INF=1./(1.+pow((Cai/0.00035),16));
      tau_g = taug;                            
      double d_g_g = - 1.0/tau_g;
      if ( (G_INF > sg) && (V > -37) ) gv_uv[15][15] = 0.0;
      else gv_uv[15][15] = d_g_g; 

      return gv_uv;
    }

 
 
 
 
 
 
 
   //----------------------------- Ten Tusscher, 18 gating variables -------------------------------- 
  
  namespace {
    
    TenTusscher18::Gating vRest18()
    {
      TenTusscher18::Gating v;
      
      /*
      v[0] =   0.000080;			// Cai
      v[1] =   0.56;				// CaSR
      v[2] =  11.6;					// Nai
      v[3] = 138.3;					// Ki
      v[4] =   0.0;					// M
      v[5] =   0.75;				// H
      v[6] =   0.75;				// J
      v[7] =   0.0;					// Xr1
      v[8] =   1.0;					// Xr2
      v[9] =   0.0;					// Xs
      v[10] =  0.0;					// R
      v[11] =  1.0;					// S
      v[12] =  0.0;					// D
      v[13] =  1.0;					// F
      v[14] =  1.0;					// FCa
      v[15] =  1.0;					// F2
      v[16] =  0.00036;				// CaSS
      v[17] =  1.0;					// RR
      */
      
      v[0] =   0.000126;			// Cai
      v[1] =   3.64;		    	// CaSR
      v[2] =   8.604;				// Nai
      v[3] = 136.89;				// Ki
      v[4] =   0.00172;		    	// M
      v[5] =   0.7444;			    // H
      v[6] =   0.7045;	 		    // J
      v[7] =   0.00621;			    // Xr1
      v[8] =   0.4712;			    // Xr2
      v[9] =   0.0095;			    // Xs
      v[10] =  2.42e-8;			    // R
      v[11] =  0.999998;			// S
      v[12] =  3.373e-5;			// D
      v[13] =  0.7888;			    // F
      v[14] =  0.9953;				// FCa
      v[15] =  0.9755;				// F2
      v[16] =  0.00036;				// CaSS
      v[17] =  0.9073;				// RR
      
      return v;
    }
  }
  
  TenTusscher18::TenTusscher18(): MembraneModelBase<TenTusscher18,18>("ten Tusscher 18")      // -8.523e+1
  {
    //this->setRestState(-8.62e+1,vRest18());
    std::cout << "TenTusscher18::TenTusscher18(): MembraneModelBase<TenTusscher18,18>" << std::endl;
    this->setRestState(-86.4956,vRest18());
  }
    
  double TenTusscher18::current(double u, Gating const& v) const
  {

    // State variables
    double    V = u;
    double  Cai = v[0];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];
    double  sf2 = v[15];
    double CaSS = v[16];

    //External concentrations
    double Ko  = 5.4;
    double Cao = 2.0;
    double Nao = 140.0;

    //Intracellular calcium flux dynamics


    //Constants
    double R     = 8314.472;
    double F     = 96485.3415;
    double T     = 310.0;
    double RTONF = (R*T)/F;

    //Parameters for currents
    double Gkr  = 0.153;
    double pKNa = 0.03;
    double Gks  = 0.392;

    double GK1 = 5.405;
    double Gto = 0.294;

    double GNa   = 14.838;
    double GbNa  = 0.00029;
    double KmK   = 1.0;
    double KmNa  = 40.0;
    double knak  = 2.724;
    double GCaL  = 0.00003980;  
    double GbCa  = 0.000592;
    double knaca = 1000;
    double KmNai = 87.5;
    double KmCa  = 1.38;
    double ksat  = 0.1;
    double n     = 0.35;
    double GpCa  = 0.1238;
    double KpCa  = 0.0005;
    double GpK   = 0.0146;
 
        
    //Needed to compute currents
    double Ek=RTONF*(log((Ko/Ki)));
    double Ena=RTONF*(log((Nao/Nai)));
    double Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    double Eca=0.5*RTONF*(log((Cao/Cai)));
    double Ak1=0.1/(1.+exp(0.06*(V-Ek-200)));
    double Bk1=(3.*exp(0.0002*(V-Ek+100))+
	       exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
    double rec_iK1=Ak1/(Ak1+Bk1);
    double rec_iNaK=(1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
    double rec_ipK=1./(1.+exp((25-V)/5.98));


    //Compute currents
    double INa=GNa*sm*sm*sm*sh*sj*(V-Ena);
    double ICaL=GCaL*sd*sf*sf2*sfca*4*(V-15)*(F*F/(R*T))*
               (0.25*exp(2*(V-15)*F/(R*T))*CaSS-Cao)/(exp(2*(V-15)*F/(R*T))-1.);
    double Ito=Gto*sr*ss*(V-Ek);
    double IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(V-Ek);
    double IKs=Gks*sxs*sxs*(V-Eks);
    double IK1=GK1*rec_iK1*(V-Ek);
    double INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*exp((n-1)*V*F/(R*T))))*
      (exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-
       exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    double INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    double IpCa=GpCa*Cai/(KpCa+Cai);
    double IpK=GpK*rec_ipK*(V-Ek);
    double IbNa=GbNa*(V-Ena);
    double IbCa=GbCa*(V-Eca);
    
    
    double dV_dt = -(IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK);
//    std::cout << "INaK = " << INaK << " INaCa = " << INaCa << " IpCa = " << IpCa << "    IpK = " << IpK << std::endl;  
//    std::cout << "INa = " << INa << " IbNa = " << IbNa << " ICaL = " << ICaL << "    IbCa = " << IbCa << std::endl;  
//    std::cout << "IKr = " << IKr << " IKs = " << IKs << " IK1 = " << IK1 << "    Ito = " << Ito << std::endl;  
//    std::cout << " Igesamt = " << dV_dt << std::endl;

    if (std::isnan(dV_dt))  std::cout << "dV_dt = NaN" << std::endl;

    return dV_dt;
  }


  double TenTusscher18::current_du(double u, Gating const& v, double h) const
  {
    double dV_du;
      
      
    // numerical
/*
      double fuv  = (*this).current(u-h,v);
      double fuiv = (*this).current(u+h,v);
      dV_du = (fuiv-fuv)/(2*h);
*/      

    // analytical

    // State variables
    double    V = u;
    double  Cai = v[0];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];
    double  sf2 = v[15];
    double CaSS = v[16];

    //External concentrations
    double Ko  = 5.4;
    double Cao = 2.0;
    double Nao = 140.0;

    //Constants
    double R     = 8314.472;
    double F     = 96485.3415;
    double T     = 310.0;
    double RTONF = (R*T)/F;

    //Parameters for currents
    double Gkr  = 0.153;
    double Gks  = 0.392;

    double GK1 = 5.405;
    double Gto = 0.294;

    double GNa   = 14.838;
    double GbNa  = 0.00029;
    double KmK   = 1.0;
    double KmNa  = 40.0;
    double knak  = 2.724;
    double GCaL  = 0.00003980;  
    double GbCa  = 0.000592;
    double knaca = 1000;
    double KmNai = 87.5;
    double KmCa  = 1.38;
    double ksat  = 0.1;
    double n     = 0.35;
    double GpK   = 0.0146;
 
    // Compute
    
        
    //Needed to compute currents
    double c1 = 0, c2 = 0, c3 = 0;
    double termA, termB, termC;
    double termA_u, termB_u, termC_u;

    double Ek  = RTONF*(log((Ko/Ki)));
    
    double Ak1    = 0.1/(1.+exp(0.06*(V-Ek-200)));
    double Ak1_u  = -0.1*0.06*exp(0.06*(V-Ek-200)) / ( (1.+exp(0.06*(V-Ek-200))) * (1.+exp(0.06*(V-Ek-200))) );
    double Bk11   = 3.*exp(0.0002*(V-Ek+100)) + exp(0.1*(V-Ek-10));
    double Bk11_u = 0.0002*3.*exp(0.0002*(V-Ek+100)) + 0.1*exp(0.1*(V-Ek-10));
    double Bk12   = 1.+exp(-0.5*(V-Ek));
    double Bk12_u = -0.5* exp(-0.5*(V-Ek));
    double Bk1    = Bk11/Bk12;
	double Bk1_u  = (Bk11_u*Bk12 - Bk11*Bk12_u) / (Bk12*Bk12);
	              
    double rec_iK1    = Ak1/(Ak1+Bk1);
    double rec_iK1_u  = (Ak1_u*(Ak1+Bk1) - Ak1*(Ak1_u+Bk1_u)) / ( (Ak1+Bk1) * (Ak1+Bk1) ) ;
    
    c1      = F/(R*T);
    termA   = 0.1245*exp(-0.1*c1*V);
    termA_u = -0.1*c1* 0.1245*exp(-0.1*c1*V);
    termB   = 0.0353*exp(-c1*V);
    termB_u = -0.1*c1*0.1245*exp(-0.1*c1*V);
    double rec_iNaK_u = -(termA_u + termB_u)/ ( (1. + termA + termB)*(1. + termA + termB) );
    
    double rec_ipK    = 1./(1.+exp((25-V)/5.98));
    double rec_ipK_u  = 1./5.98 * exp((25-V)/5.98) / ( (1.+exp((25-V)/5.98)) * (1.+exp((25-V)/5.98)) );

    //Compute currents
    
    double INa_u = GNa*sm*sm*sm*sh*sj;
    
    
    c1      = GCaL*sd*sf*sf2*sfca*4*(F*F/(R*T));
    c2      = 2*F/(R*T);
    termA   = 0.25*exp(c2*(V-15))*CaSS-Cao;
    termA_u = 0.25*CaSS*c2*exp(c2*(V-15));
    termB   = exp(c2*(V-15))-1.;
    termB_u = c2*exp(c2*(V-15));
    double ICaL_u = c1*termA/termB + c1*(V-15) * (termA_u*termB - termA*termB_u)/(termB*termB);

               
    double Ito_u = Gto*sr*ss;


    double IKr_u = Gkr*sqrt(Ko/5.4)*sxr1*sxr2;
    
    
    double IKs_u = Gks*sxs*sxs;
    
    
    double IK1_u = GK1*( rec_iK1 + rec_iK1_u*(V-Ek) );
    
    
    c1    = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao));
    c2    = (n-1)*F/(R*T);
    c3    = n*F/(R*T);

    termA   = 1./(1+ksat*exp(c2*V));
    termB   = exp(c3*V)*Nai*Nai*Nai*Cao;
    termC   = exp(c2*V)*Nao*Nao*Nao*Cai*2.5;

    termA_u = -(c2*ksat*exp(c2*V)) / ( (1+ksat*exp(c2*V))* (1+ksat*exp(c2*V)) );
    termB_u = c3*exp(c3*V)*Nai*Nai*Nai*Cao;
    termC_u = c2*exp(c2*V)*Nao*Nao*Nao*Cai*2.5;

    double INaCa_u = c1 * ( termA_u*(termB - termC) + termA*(termB_u - termC_u));

    c1            = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa));                    
    double INaK_u = c1*rec_iNaK_u;
    
    double IpCa_u = 0;
    
    double IpK_u = GpK*( rec_ipK_u*(V-Ek) + GpK*rec_ipK );
    
    double IbNa_u = GbNa;
    
    double IbCa_u = GbCa;
    
    
    dV_du = -(IKr_u + IKs_u + IK1_u + Ito_u + INa_u + IbNa_u + ICaL_u + IbCa_u + INaK_u + INaCa_u + IpCa_u + IpK_u);
    
    //std::cout << "dV_du = " << dV_du << "   dV_du_diffQ = " << dV_du_diffQ << std::endl;



    if (std::isnan(dV_du))
      {
        std::cout << "dV_du = NaN" << std::endl;
      }



    return dV_du;

    }
    

  TenTusscher18::Gating TenTusscher18::current_dv(double u, Gating const& v) const
    {
      double f = (*this).current(u,v);
      Gating fv_uv, dir;
      for (int i=0; i<nGating; ++i)
      {
        dir = 0;
        dir[i] = 1e-5;
        fv_uv[i] = ((*this).current(u,v+dir)-f)/1e-5;
      }
      return fv_uv;
    }
  
  
  
  double TenTusscher18::gatingRhs1(double u, Gating const& v, int vId, double Istim) const
  {
    // RHS
    Gating dv_dt;
    
    // State variables
    double    V = u;
    double  Cai = v[0];
    double CaSR = v[1];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];
    double  sf2 = v[15];
    double CaSS = v[16];
    double  srr = v[17];
    
    //External concentrations
    double Ko=5.4;
    double Cao=2.0;
    double Nao=140.0;
    
    //Intracellular volumes
    double Vc=0.016404;
    double Vsr=0.001094;
    double Vss=0.00005468;
    
    //Calcium buffering dynamics
    double Bufc=0.2;
    double Kbufc=0.001;
    double Bufsr=10.;
    double Kbufsr=0.3;
    double Bufss=0.4;
    double Kbufss=0.00025;
    
    //Intracellular calcium flux dynamics
    double Vmaxup=0.006375;
    double Kup=0.00025;
    double Vrel=0.102;     //40.8;
    double k1_=0.15;
    double k2_=0.045;
    double k3=0.060;
    double k4=0.005;      //0.000015;
    double EC=1.5;
    double maxsr=2.5;
    double minsr=1.;
    double Vleak=0.00036;
    double Vxfer=0.0038;
    
    
    
    //Constants
    double R=8314.472;
    double F=96485.3415;
    double T=310.0;
    double RTONF=(R*T)/F;
    
    //Cellular capacitance         
    double CAPACITANCE=0.185;
    
    //Parameters for currents
    double Gkr=0.153;
    double pKNa=0.03;
    double Gks=0.392;
    
    double GK1=5.405;
    double Gto=0.294;
    
    double GNa=14.838;
    double GbNa=0.00029;
    double KmK=1.0;
    double KmNa=40.0;
    double knak=2.724;
    double GCaL=0.00003980;  
    double GbCa=0.000592;
    double knaca=1000;
    double KmNai=87.5;
    double KmCa=1.38;
    double ksat=0.1;
    double n=0.35;
    double GpCa=0.1238;
    double KpCa=0.0005;
    double GpK=0.0146;
    
    // Compute
    double inverseVcF2=1/(2*Vc*F);
    double inverseVcF=1./(Vc*F);
    double inversevssF2=1/(2*Vss*F);
    
    
    //Needed to compute currents
    double Ek=RTONF*(log((Ko/Ki)));
    double Ena=RTONF*(log((Nao/Nai)));
    double Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    double Eca=0.5*RTONF*(log((Cao/Cai)));
    
    double Ak1=0.1/(1.+exp(0.06*(V-Ek-200)));
    double Bk1=(3.*exp(0.0002*(V-Ek+100))+
    exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
    
    double rec_iK1=Ak1/(Ak1+Bk1);
    double rec_iNaK=(1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
    double rec_ipK=1./(1.+exp((25-V)/5.98));
    
    
    //Compute currents
    double INa=GNa*sm*sm*sm*sh*sj*(V-Ena);
    double ICaL=GCaL*sd*sf*sf2*sfca*4*(V-15)*(F*F/(R*T))*
    (0.25*exp(2*(V-15)*F/(R*T))*CaSS-Cao)/(exp(2*(V-15)*F/(R*T))-1.);
    double Ito=Gto*sr*ss*(V-Ek);
    double IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(V-Ek);
    double IKs=Gks*sxs*sxs*(V-Eks);
    double IK1=GK1*rec_iK1*(V-Ek);
    double INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
    (1./(1+ksat*exp((n-1)*V*F/(R*T))))*
    (exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-
    exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    double INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    double IpCa=GpCa*Cai/(KpCa+Cai);
    double IpK=GpK*rec_ipK*(V-Ek);
    double IbNa=GbNa*(V-Ena);
    double IbCa=GbCa*(V-Eca);
    
    
    //-------------------------------------------------------------------------------------------   
    double kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))); 
    double k1 = k1_/kCaSR;
    double k2 = k2_*kCaSR;
    
    double srr_dt;
    double CaSR_dt;
    double Cai_dt;
    double Nai_dt;
    double Ki_dt;
    double sm_dt;
    double sh_dt;
    double sj_dt;
    double sxr1_dt;
    double sxr2_dt;
    double sxs_dt; 
    double ss_dt;  
    double sr_dt;  
    double sd_dt;  
    double sf_dt;  
    double sf2_dt; 
    double sfca_dt;
    
    
    
    if (vId == 17)
    {
      srr_dt = k4*(1-srr)-k2*CaSS*srr;
      return srr_dt;
    }  
    else
    {
      
      double sOO   = k1*CaSS*CaSS*srr/(k3+k1*CaSS*CaSS);
      double Irel  = Vrel*sOO*(CaSR-CaSS);
      double Ileak = Vleak*(CaSR-Cai);
      double Iup   = Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
      double Ixfer = Vxfer*(CaSS-Cai);
      
      double CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
      if (vId == 1)
      {
        CaSR_dt = Iup-Irel-Ileak;
        CaSR_dt *= CaCSQN;
        return CaSR_dt;
      }
      else
      {
        double CaSSBuf = Bufss*CaSS/(CaSS+Kbufss);
        
        if (vId == 16)
        {
          double CaSS_dt = -Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE);
          CaSS_dt *= CaSSBuf;
          return CaSS_dt;
        }
        else
        {
          double CaBuf   = Bufc*Cai/(Cai+Kbufc);
          
          if (vId == 0)
          {
            Cai_dt  = (-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer;
            Cai_dt  *= CaBuf;
            return Cai_dt;
          }
          else
          {
            if (vId == 2)
            {
              Nai_dt  = -(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
              return Nai_dt;
            }
            else
            {
              if (vId == 3)
              {
                Ki_dt   = -(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
                return Ki_dt;
              }
              else
              {
                
                //compute steady state values and time constants 
                double AH_1, BH_1, AH_2, BH_2;
                double AJ_1, BJ_1, AJ_2, BJ_2;
                double TAU_H, TAU_J;
                
                double AM=1./(1.+exp((-60.-V)/5.));
                double BM=0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
                double TAU_M=AM*BM;
                double M_INF=1./((1.+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));
                
                if (vId == 4)
                {
                  sm_dt  = (M_INF-sm)/TAU_M;
                  return sm_dt;
                }
                else
                {
                  if (V>=-40.)
                  {
                    AH_1=0.; 
                    BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
                    TAU_H= 1.0/(AH_1+BH_1);
                  }
                  else
                  {
                    AH_2=(0.057*exp(-(V+80.)/6.8));
                    BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
                    TAU_H=1.0/(AH_2+BH_2);
                  }
                  double H_INF = 1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));
                  
                  if (vId == 5)
                  {
                    sh_dt   = (H_INF-sh)/TAU_H;
                    return sh_dt;
                  }
                  else
                  {
                    if(V>=-40.)
                    {
                      AJ_1=0.;      
                      BJ_1=(0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
                      TAU_J= 1.0/(AJ_1+BJ_1);
                    }
                    else
                    {
                      AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
                      exp(-0.04391*V))*(V+37.78)/
                      (1.+exp(0.311*(V+79.23))));    
                      BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
                      TAU_J= 1.0/(AJ_2+BJ_2);
                    }
                    double J_INF=H_INF;
                    
                    if (vId == 6)
                    {
                      sj_dt   = (J_INF-sj)/TAU_J;  
                      return sj_dt;
                    }
                    else
                    {
                      double Xr1_INF=1./(1.+exp((-26.-V)/7.));
                      double axr1=450./(1.+exp((-45.-V)/10.));
                      double bxr1=6./(1.+exp((V-(-30.))/11.5));
                      double TAU_Xr1=axr1*bxr1;
                      double Xr2_INF=1./(1.+exp((V-(-88.))/24.));
                      double axr2=3./(1.+exp((-60.-V)/20.));
                      double bxr2=1.12/(1.+exp((V-60.)/20.));
                      double TAU_Xr2=axr2*bxr2;
                      
                      if (vId == 7)
                      {
                        sxr1_dt = (Xr1_INF-sxr1)/TAU_Xr1; 
                        return sxr1_dt;
                      }
                      else
                      {
                        if (vId == 8)
                        {
                          sxr2_dt = (Xr2_INF-sxr2)/TAU_Xr2; 
                          return sxr2_dt;
                        }
                        else
                        {
                          double Xs_INF=1./(1.+exp((-5.-V)/14.));
                          double Axs=(1400./(sqrt(1.+exp((5.-V)/6))));
                          double Bxs=(1./(1.+exp((V-35.)/15.)));
                          double TAU_Xs=Axs*Bxs+80;
                          
                          if (vId == 9)
                          {
                            sxs_dt  = (Xs_INF-sxs)/TAU_Xs; 
                            return sxs_dt;
                          }
                          else
                          {
                            double R_INF=1./(1.+exp((20-V)/6.));
                            double TAU_R=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;
                            
                            if (vId == 10)
                            {
                              sr_dt   = (R_INF-sr)/TAU_R;
                              return sr_dt;
                            }
                            else
                            {
                              double S_INF=1./(1.+exp((V+20)/5.));
                              double TAU_S=85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;
                              
                              if (vId == 11)
                              {
                                ss_dt = (S_INF-ss)/TAU_S;
                                return ss_dt;
                              }
                              else
                              {
                                double D_INF=1./(1.+exp((-8-V)/7.5));
                                double Ad=1.4/(1.+exp((-35-V)/13))+0.25;
                                double Bd=1.4/(1.+exp((V+5)/5));
                                double Cd=1./(1.+exp((50-V)/20));
                                double TAU_D=Ad*Bd+Cd;
                                
                                if (vId == 12)
                                {
                                  sd_dt   = (D_INF-sd)/TAU_D;
                                  return sd_dt;
                                }
                                else  
                                {
                                  double F_INF=1./(1.+exp((V+20)/7));
                                  double Af=1102.5*exp(-(V+27)*(V+27)/225);
                                  double Bf=200./(1+exp((13-V)/10.));
                                  double Cf=(180./(1+exp((V+30)/10)))+20;
                                  double TAU_F=Af+Bf+Cf;
                                  
                                  if (vId == 13)
                                  {
                                    sf_dt   = (F_INF-sf)/TAU_F;
                                    return sf_dt;
                                  }
                                  else 
                                  {
                                    double F2_INF=0.67/(1.+exp((V+35)/7))+0.33;
                                    double Af2=600*exp(-(V+25)*(V+25)/170);
                                    double Bf2=31/(1.+exp((25-V)/10));
                                    double Cf2=16/(1.+exp((V+30)/10));
                                    double TAU_F2=Af2+Bf2+Cf2;
                                    
                                    if (vId == 15)
                                    {
                                      sf2_dt  = (F2_INF-sf2)/TAU_F2;
                                      return sf2_dt;
                                    }
                                    else
                                    {
                                      double FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
                                      double TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;
                                      
                                      if (vId == 14)
                                      {
                                        sfca_dt = (FCaSS_INF-sfca)/TAU_FCaSS;
                                        return sfca_dt;
                                      }
                                      
                                    }   // else
                                  }   // else
                                }   // else
                              }   // else
                            }   // else
                          }   // else
                        }   // else
                      }   // else
                    }   // else
                  }   // else
                }   // else
              }   // else
            }   // else
          }   // else
        }   // else
      }   // else
    }   // else
    
    // Never get here
    abort();
    return 0;
  }
  
  
  
  TenTusscher18::Gating TenTusscher18::gatingRhs(double u, Gating const& v, double Istim) const
  {
 
    // State variables
    double    V = u;
    double  Cai = v[0];
    double CaSR = v[1];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];
    double  sf2 = v[15];
    double CaSS = v[16];
    double  srr = v[17];

    //External concentrations
    double Ko=5.4;
    double Cao=2.0;
    double Nao=140.0;

    //Intracellular volumes
    double Vc=0.016404;
    double Vsr=0.001094;
    double Vss=0.00005468;

    //Calcium buffering dynamics
    double Bufc=0.2;
    double Kbufc=0.001;
    double Bufsr=10.;
    double Kbufsr=0.3;
    double Bufss=0.4;
    double Kbufss=0.00025;

    //Intracellular calcium flux dynamics
    double Vmaxup=0.006375;
    double Kup=0.00025;
    double Vrel=0.102;     //40.8;
    double k1_=0.15;
    double k2_=0.045;
    double k3=0.060;
    double k4=0.005;      //0.000015;
    double EC=1.5;
    double maxsr=2.5;
    double minsr=1.;
    double Vleak=0.00036;
    double Vxfer=0.0038;



    //Constants
    double R=8314.472;
    double F=96485.3415;
    double T=310.0;
    double RTONF=(R*T)/F;

    //Cellular capacitance         
    double CAPACITANCE=0.185;

    //Parameters for currents
    double Gkr=0.153;
    double pKNa=0.03;
    double Gks=0.392;

    double GK1=5.405;
    double Gto=0.294;

    double GNa=14.838;
    double GbNa=0.00029;
    double KmK=1.0;
    double KmNa=40.0;
    double knak=2.724;
    double GCaL=0.00003980;  
    double GbCa=0.000592;
    double knaca=1000;
    double KmNai=87.5;
    double KmCa=1.38;
    double ksat=0.1;
    double n=0.35;
    double GpCa=0.1238;
    double KpCa=0.0005;
    double GpK=0.0146;
 
    // Compute
    double inverseVcF2=1/(2*Vc*F);
    double inverseVcF=1./(Vc*F);
    double inversevssF2=1/(2*Vss*F);
    
    
    //Needed to compute currents
    double Ek=RTONF*(log((Ko/Ki)));
    double Ena=RTONF*(log((Nao/Nai)));
    double Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    double Eca=0.5*RTONF*(log((Cao/Cai)));
    if (std::isnan(Eca))
      {
        std::cout << "Eca = NaN:  Cai = " << Cai << std::endl;
      }
    
    double Ak1=0.1/(1.+exp(0.06*(V-Ek-200)));
    double Bk1=(3.*exp(0.0002*(V-Ek+100))+
	       exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
	       
    double rec_iK1=Ak1/(Ak1+Bk1);
    double rec_iNaK=(1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
    double rec_ipK=1./(1.+exp((25-V)/5.98));


    //Compute currents
    double INa=GNa*sm*sm*sm*sh*sj*(V-Ena);
    double ICaL=GCaL*sd*sf*sf2*sfca*4*(V-15)*(F*F/(R*T))*
               (0.25*exp(2*(V-15)*F/(R*T))*CaSS-Cao)/(exp(2*(V-15)*F/(R*T))-1.);
//     std::cout << "ICaL = " << ICaL << "   GCaL = " << GCaL 
//               << "   sd = " << sd << "   sf = " << sf << "   sf2 = " << sf2 << "   sfca = " << sfca 
//               << "   CaSS = " << CaSS << "   V = " << V  <<std::endl;

    double Ito=Gto*sr*ss*(V-Ek);
    double IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(V-Ek);
    double IKs=Gks*sxs*sxs*(V-Eks);
    double IK1=GK1*rec_iK1*(V-Ek);
    double INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*exp((n-1)*V*F/(R*T))))*
      (exp(n*V*F/(R*T))*Nai*Nai*Nai*Cao-
       exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5);

    double INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    double IpCa=GpCa*Cai/(KpCa+Cai);
    double IpK=GpK*rec_ipK*(V-Ek);
    double IbNa=GbNa*(V-Ena);
    double IbCa=GbCa*(V-Eca);
    if (std::isnan(IbCa))
      {
        std::cout << "IbCa" << " = NaN" << std::endl;
      }
    
    
    //-------------------------------------------------------------------------------------------   

    double kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))); 
    double k1 = k1_/kCaSR;
    double k2 = k2_*kCaSR;
    double srr_dt = k4*(1-srr)-k2*CaSS*srr;
    
    double sOO   = k1*CaSS*CaSS*srr/(k3+k1*CaSS*CaSS);
    double Irel  = Vrel*sOO*(CaSR-CaSS);
    double Ileak = Vleak*(CaSR-Cai);
    double Iup   = Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
    double Ixfer = Vxfer*(CaSS-Cai);
    
    double CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    double CaSR_dt = Iup-Irel-Ileak;
    CaSR_dt *= CaCSQN;     // korektur

    double CaSSBuf = Bufss*CaSS/(CaSS+Kbufss);
    double CaSS_dt = -Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE);
    CaSS_dt *= CaSSBuf;    // korektur

    double CaBuf   = Bufc*Cai/(Cai+Kbufc);
    double Cai_dt  = (-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer;
    Cai_dt *= CaBuf;       // korektur

    double Nai_dt  = -(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;

    double Ki_dt   = -(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;


    //compute steady state values and time constants 
    double AH_1, BH_1, AH_2, BH_2;
    double AJ_1, BJ_1, AJ_2, BJ_2;
    double TAU_H, TAU_J;
    
    double AM=1./(1.+exp((-60.-V)/5.));
    double BM=0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
    double TAU_M=AM*BM;
    double M_INF=1./((1.+exp((-56.86-V)/9.03))*(1.+exp((-56.86-V)/9.03)));
    if (V>=-40.)
    {
	  AH_1=0.; 
      BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
      TAU_H= 1.0/(AH_1+BH_1);
    }
    else
    {
	  AH_2=(0.057*exp(-(V+80.)/6.8));
	  BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
	  TAU_H=1.0/(AH_2+BH_2);
    }
    double H_INF = 1./((1.+exp((V+71.55)/7.43))*(1.+exp((V+71.55)/7.43)));
    if(V>=-40.)
    {
	  AJ_1=0.;      
	  BJ_1=(0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
	  TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
	  AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
		     exp(-0.04391*V))*(V+37.78)/
	         (1.+exp(0.311*(V+79.23))));    
	  BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
	  TAU_J= 1.0/(AJ_2+BJ_2);
    }
    double J_INF=H_INF;

    double Xr1_INF=1./(1.+exp((-26.-V)/7.));
    double axr1=450./(1.+exp((-45.-V)/10.));
    double bxr1=6./(1.+exp((V-(-30.))/11.5));
    double TAU_Xr1=axr1*bxr1;
    double Xr2_INF=1./(1.+exp((V-(-88.))/24.));
    double axr2=3./(1.+exp((-60.-V)/20.));
    double bxr2=1.12/(1.+exp((V-60.)/20.));
    double TAU_Xr2=axr2*bxr2;

    double Xs_INF=1./(1.+exp((-5.-V)/14.));
    double Axs=(1400./(sqrt(1.+exp((5.-V)/6))));
    double Bxs=(1./(1.+exp((V-35.)/15.)));
    double TAU_Xs=Axs*Bxs+80;
    
    double R_INF=1./(1.+exp((20-V)/6.));
    double S_INF=1./(1.+exp((V+20)/5.));
    double TAU_R=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;
    double TAU_S=85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;


    double D_INF=1./(1.+exp((-8-V)/7.5));
    double Ad=1.4/(1.+exp((-35-V)/13))+0.25;
    double Bd=1.4/(1.+exp((V+5)/5));
    double Cd=1./(1.+exp((50-V)/20));
    double TAU_D=Ad*Bd+Cd;
    double F_INF=1./(1.+exp((V+20)/7));
    double Af=1102.5*exp(-(V+27)*(V+27)/225);
    double Bf=200./(1+exp((13-V)/10.));
    double Cf=(180./(1+exp((V+30)/10)))+20;
    double TAU_F=Af+Bf+Cf;
    double F2_INF=0.67/(1.+exp((V+35)/7))+0.33;
    double Af2=600*exp(-(V+25)*(V+25)/170);
    double Bf2=31/(1.+exp((25-V)/10));
    double Cf2=16/(1.+exp((V+30)/10));
    double TAU_F2=Af2+Bf2+Cf2;
    double FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
    double TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;

     //Update gates
     double sm_dt   = (M_INF-sm)/TAU_M;
     double sh_dt   = (H_INF-sh)/TAU_H;
     double sj_dt   = (J_INF-sj)/TAU_J;
     double sxr1_dt = (Xr1_INF-sxr1)/TAU_Xr1; 
     double sxr2_dt = (Xr2_INF-sxr2)/TAU_Xr2;
     double sxs_dt  = (Xs_INF-sxs)/TAU_Xs;
     double ss_dt   = (S_INF-ss)/TAU_S;
     double sr_dt   = (R_INF-sr)/TAU_R;
     double sd_dt   = (D_INF-sd)/TAU_D; 
     double sf_dt   = (F_INF-sf)/TAU_F; 
     double sf2_dt  = (F2_INF-sf2)/TAU_F2; 
     double sfca_dt = (FCaSS_INF-sfca)/TAU_FCaSS;
     



    // RHS
    Gating dv_dt;
    
    dv_dt[0]  = Cai_dt;
    dv_dt[1]  = CaSR_dt;
    dv_dt[2]  = Nai_dt;
    dv_dt[3]  = Ki_dt;
    dv_dt[4]  = sm_dt;
    dv_dt[5]  = sh_dt;
    dv_dt[6]  = sj_dt;
    dv_dt[7]  = sxr1_dt;
    dv_dt[8]  = sxr2_dt;
    dv_dt[9]  = sxs_dt;
    dv_dt[10] = sr_dt;
    dv_dt[11] = ss_dt;
    dv_dt[12] = sd_dt;
    dv_dt[13] = sf_dt;
    dv_dt[14] = sfca_dt;
    dv_dt[15] = sf2_dt;
    dv_dt[16] = CaSS_dt;
    dv_dt[17] = srr_dt;
    
    
    for (int i=0; i<nGating; ++i)
    {
      if (std::isnan(dv_dt[i]))
      {
        // throw DetailedException("NaN occured in ten Tusscher computation",__FILE__,__LINE__);
        std::cout << "dv_dt[" << i << "] = NaN" << std::endl;
      }
    }
    
    return dv_dt;
  }
  

  TenTusscher18::Gating TenTusscher18::gatingRhs_du(double u, Gating const& v) const
    {
      Gating guA = (*this).gatingRhs(u-1e-5,v);
      Gating guB = (*this).gatingRhs(u+1e-5,v);
      return (guB-guA)/2e-5;
    }
    

  TenTusscher18::GatingJacobian TenTusscher18::gatingRhs_dv(double u, Gating const& v) const
    {
      GatingJacobian gv_uv = 0.;
      

/*
      // numerical differentiation in diagonal elements of the matrix
      Gating guv = (*this).gatingRhs(u,v);
      Gating dir;
      double col;
      for (int i=0; i<nGating; ++i)
      {
        dir = 0;
        dir[i]=dr;     
        col = (*this).gatingRhs1(u,v+dir,i);
        gv_uv[i][i] = (col - guv[i])/dr;
      }
*/

/*
      // numerical differentiation in all elements of the matrix
      Gating col, dir, guv = (*this).gatingRhs(u,v);
      for (int i=0; i<nGating; ++i)
      {
        dir = 0;
        dir[i]=dr;      // 1e-8;
        col = ((*this).gatingRhs(u,v+dir)-guv)/dr;
        for (int j=0; j<nGating; ++j) gv_uv[j][i] = col[j];
      }

*/


      
  
    // analytical evaluation in diagonal 
    // State variables
    double    V = u;
    double  Cai = v[0];
    double CaSR = v[1];
    double  Nai = v[2];
    double   Ki = v[3];
    double   sm = v[4]; 
    double   sh = v[5];
    double   sj = v[6];
    double sxr1 = v[7]; 
    double sxr2 = v[8];
    double  sxs = v[9];
    double   sr = v[10];
    double   ss = v[11];
    double   sd = v[12];
    double   sf = v[13];
    double sfca = v[14];
    double  sf2 = v[15];
    double CaSS = v[16];
    double  srr = v[17];

    //External concentrations
    double Ko=5.4;
    double Cao=2.0;
    double Nao=140.0;

    //Intracellular volumes
    double Vc=0.016404;
    double Vsr=0.001094;
    double Vss=0.00005468;

    //Calcium buffering dynamics
    double Bufc=0.2;
    double Kbufc=0.001;
    double Bufsr=10.;
    double Kbufsr=0.3;
    double Bufss=0.4;
    double Kbufss=0.00025;

    //Intracellular calcium flux dynamics
    double Vmaxup=0.006375;
    double Kup=0.00025;
    double Vrel=0.102;     //40.8;
    double k1_=0.15;
    double k2_=0.045;
    double k3=0.060;
    double k4=0.005;      //0.000015;
    double EC=1.5;
    double maxsr=2.5;
    double minsr=1.;
    double Vleak=0.00036;
    double Vxfer=0.0038;



    //Constants
    double R=8314.472;
    double F=96485.3415;
    double T=310.0;
    double RTONF=(R*T)/F;

    //Cellular capacitance         
    double CAPACITANCE=0.185;

    //Parameters for currents
    double Gkr=0.153;
    double pKNa=0.03;
    double Gks=0.392;

    double GK1=5.405;
    double Gto=0.294;

    double GNa=14.838;
    double GbNa=0.00029;
    double KmK=1.0;
    double KmNa=40.0;
    double knak=2.724;
    double GCaL=0.00003980;  
    double GbCa=0.000592;
    double knaca=1000;
    double KmNai=87.5;
    double KmCa=1.38;
    double ksat=0.1;
    double n=0.35;
    double GpCa=0.1238;
    double KpCa=0.0005;
    double GpK=0.0146;
 
    // Compute
    double inverseVcF2=1/(2*Vc*F);
    double inverseVcF=1./(Vc*F);
    double inversevssF2=1/(2*Vss*F);
    
    

    //Needed to compute currents
    double Ek=RTONF*(log((Ko/Ki)));
    
    double Ak1=0.1/(1.+exp(0.06*(V-Ek-200)));
    double Bk1=(3.*exp(0.0002*(V-Ek+100))+
	       exp(0.1*(V-Ek-10)))/(1.+exp(-0.5*(V-Ek)));
	       
    double rec_iK1=Ak1/(Ak1+Bk1);
    double rec_iNaK=(1./(1.+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));
    double rec_ipK=1./(1.+exp((25-V)/5.98));


    
    //-------------------------------------------------------------------------------------------   
    double kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))); 
    double k1 = k1_/kCaSR;
    double k2 = k2_*kCaSR;
    //double srr_dt = k4*(1-srr)-k2*CaSS*srr;
    
    // sRR
    gv_uv[17][17] = -k4 - k2*CaSS;
    
    
    double sOO   = k1*CaSS*CaSS*srr/(k3+k1*CaSS*CaSS);

    double CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    double CaSSBuf = Bufss*CaSS/(CaSS+Kbufss);

    double c1 = 0, c2 = 0, c3 = 0;
    
    // Cai
    c1 = 0.5*RTONF*GbCa;
    double IbCa_Cai  = c1/Cai;
    c1 = GpCa;
    c2 = KpCa;
    double IpCa_Cai  = c1*c2/( (c2+Cai)*(c2+Cai) ); 
    c1 = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*exp((n-1)*V*F/(R*T))));
    double INaCa_Cai = -c1* exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*2.5;  
    c1 = Vmaxup;
    c2 = Kup*Kup;
    double Iup_Cai   = 2*c1*c2*Cai/((Cai*Cai+c2)*(Cai*Cai+c2));    
    double Ileak_Cai = -Vleak;  
    double Ixfer_Cai = -Vxfer;  
    
    double CaBuf   = Bufc*Cai/(Cai+Kbufc);

    gv_uv[0][0] = (-(IbCa_Cai + IpCa_Cai - 2*INaCa_Cai)*inverseVcF2*CAPACITANCE) - (Iup_Cai - Ileak_Cai)*(Vsr/Vc) + Ixfer_Cai;
    gv_uv[0][0] *=  CaBuf;


    //CaSR
    double Iup_CaSR = 0.0;
    double Irel_CaSR = Vrel*sOO*CaSR;
    double Ileak_CaSR = Vleak;
    gv_uv[1][1] = Iup_CaSR - Irel_CaSR - Ileak_CaSR;
    gv_uv[1][1] *= CaCSQN;

    // Nai
    c1 = GNa*sm*sm*sm*sh*sj * RTONF;
    double INa_Nai  = c1/Nai;
    c1 = GbNa * RTONF;
    double IbNa_Nai = c1/Nai;
    c1 = knak*(Ko/(Ko+KmK))*rec_iNaK;
    c2 = KmNa;
    double INaK_Nai = c1*c2/( (Nai+c2)*(Nai+c2) );
    c1 = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao)) * (1./(1+ksat*exp((n-1)*V*F/(R*T))));
    c2 = exp((n-1)*V*F/(R*T))*Nao*Nao*Nao*Cai*2.5;
    double INaCa_Nai = c1 * exp(n*V*F/(R*T))*Cao * 3.0*Nai*Nai;
    
    gv_uv[2][2]  = -(INa_Nai+IbNa_Nai+3*INaK_Nai+3*INaCa_Nai)*inverseVcF*CAPACITANCE;


    // Ki     ------>  überprüfen
    c1 = GK1*rec_iK1 * RTONF;
    double IK1_Ki = c1/Ki;
    c1 = Gto*sr*ss * RTONF;
    double Ito_Ki = c1/Ki;
    c1 = Gkr*sqrt(Ko/5.4)*sxr1*sxr2 * RTONF;
    double IKr_Ki = c1/Ki;
    c1 = Gks*sxs*sxs * RTONF;
    c2 = Ko+pKNa*Nao;
    c3 = pKNa*Nai;
    double IKs_Ki = c1/(Ki+c3);   
    double INaK_Ki = 0.0;
    c1 = GpK*rec_ipK * RTONF;
    double IpK_Ki = c1/Ki;
    gv_uv[3][3]   = -(IK1_Ki+Ito_Ki+IKr_Ki+IKs_Ki-2*INaK_Ki+IpK_Ki)*inverseVcF*CAPACITANCE;


    // CaSS     ------>  überprüfen
    double Ixfer_CaSS = Vxfer;
    double Irel_CaSS  = -Vrel*sOO;
    c1 = GCaL*sd*sf*sf2*sfca*4*(V-15)*(F*F/(R*T));
    c2 = exp(2*(V-15)*F/(R*T))-1.;
    c3 = 0.25*exp(2*(V-15)*F/(R*T));
    double ICaL_CaSS  = (c1 * c3)/c2;
    gv_uv[16][16] = -Ixfer_CaSS*(Vc/Vss) + Irel_CaSS*(Vsr/Vss) + (-ICaL_CaSS*inversevssF2*CAPACITANCE);
    gv_uv[16][16] *= CaSSBuf;

    //compute steady state values and time constants 
    double AH_1, BH_1, AH_2, BH_2;
    double AJ_1, BJ_1, AJ_2, BJ_2;
    double TAU_H, TAU_J;
    
    double AM=1./(1.+exp((-60.-V)/5.));
    double BM=0.1/(1.+exp((V+35.)/5.))+0.10/(1.+exp((V-50.)/200.));
    double TAU_M=AM*BM;
    if (V>=-40.)
    {
	  AH_1=0.; 
      BH_1=(0.77/(0.13*(1.+exp(-(V+10.66)/11.1))));
      TAU_H= 1.0/(AH_1+BH_1);
    }
    else
    {
	  AH_2=(0.057*exp(-(V+80.)/6.8));
	  BH_2=(2.7*exp(0.079*V)+(3.1e5)*exp(0.3485*V));
	  TAU_H=1.0/(AH_2+BH_2);
    }
    if(V>=-40.)
    {
	  AJ_1=0.;      
	  BJ_1=(0.6*exp((0.057)*V)/(1.+exp(-0.1*(V+32.))));
	  TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
	  AJ_2 = (((-2.5428e4)*exp(0.2444*V)-(6.948e-6)*
		     exp(-0.04391*V))*(V+37.78)/
	         (1.+exp(0.311*(V+79.23))));    
	  BJ_2 = (0.02424*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))));
	  TAU_J= 1.0/(AJ_2+BJ_2);
    }

    double axr1=450./(1.+exp((-45.-V)/10.));
    double bxr1=6./(1.+exp((V-(-30.))/11.5));
    double TAU_Xr1=axr1*bxr1;
    double axr2=3./(1.+exp((-60.-V)/20.));
    double bxr2=1.12/(1.+exp((V-60.)/20.));
    double TAU_Xr2=axr2*bxr2;

    double Axs=(1400./(sqrt(1.+exp((5.-V)/6))));
    double Bxs=(1./(1.+exp((V-35.)/15.)));
    double TAU_Xs=Axs*Bxs+80;
    
    double TAU_R=9.5*exp(-(V+40.)*(V+40.)/1800.)+0.8;
    double TAU_S=85.*exp(-(V+45.)*(V+45.)/320.)+5./(1.+exp((V-20.)/5.))+3.;


    double Ad=1.4/(1.+exp((-35-V)/13))+0.25;
    double Bd=1.4/(1.+exp((V+5)/5));
    double Cd=1./(1.+exp((50-V)/20));
    double TAU_D=Ad*Bd+Cd;
    double Af=1102.5*exp(-(V+27)*(V+27)/225);
    double Bf=200./(1+exp((13-V)/10.));
    double Cf=(180./(1+exp((V+30)/10)))+20;
    double TAU_F=Af+Bf+Cf;
    double Af2=600*exp(-(V+25)*(V+25)/170);
    double Bf2=31/(1.+exp((25-V)/10));
    double Cf2=16/(1.+exp((V+30)/10));
    double TAU_F2=Af2+Bf2+Cf2;
    double TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;


     //Update gates
     
     //      double sm_dt   = (M_INF-sm)/TAU_M;
     //      double sh_dt   = (H_INF-sh)/TAU_H;
     //      double sj_dt   = (J_INF-sj)/TAU_J;
     //      double sxr1_dt = (Xr1_INF-sxr1)/TAU_Xr1; 
     //      double sxr2_dt = (Xr2_INF-sxr2)/TAU_Xr2;
     //      double sxs_dt  = (Xs_INF-sxs)/TAU_Xs;
     //      double sr_dt   = (R_INF-sr)/TAU_R;
     //      double ss_dt   = (S_INF-ss)/TAU_S;
     //      double sd_dt   = (D_INF-sd)/TAU_D; 
     //      double sf_dt   = (F_INF-sf)/TAU_F; 
     //      double sf2_dt  = (F2_INF-sf2)/TAU_F2; 
     //      double sfca_dt = (FCaSS_INF-sfca)/TAU_FCaSS;
     


    gv_uv[4][4]     =  -1.0/TAU_M;  
    gv_uv[5][5]     =  -1.0/TAU_H;  
    gv_uv[6][6]     =  -1.0/TAU_J;  
    gv_uv[7][7]     =  -1.0/TAU_Xr1;  
    gv_uv[8][8]     =  -1.0/TAU_Xr2;  
    gv_uv[9][9]     =  -1.0/TAU_Xs;  
    gv_uv[10][10]   =  -1.0/TAU_R;  
    gv_uv[11][11]   =  -1.0/TAU_S;  
    gv_uv[12][12]   =  -1.0/TAU_D;  
    gv_uv[13][13]   =  -1.0/TAU_F;  
    gv_uv[14][14]   =  -1.0/TAU_FCaSS;     
    gv_uv[15][15]   =  -1.0/TAU_F2;  



    return gv_uv;
  }
  
  
};
