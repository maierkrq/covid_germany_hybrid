#ifndef COVID_HH
#define COVID_HH 

#include "fem/fixdune.hh"
#include "fem/functional_aux.hh"
#include "fem/variables.hh"


int printCount = 10; 
    
template <class RType, class VarSet>
class AlievPanfilovEquation: public FunctionalBase<WeakFormulation>
{
  using Self=AlievPanfilovEquation<RType,VarSet>;
  public:
  using Scalar = RType;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  using OriginVars = VarSet;

  static constexpr int dim = AnsatzVars::Grid::dimension;

  using Cell = typename AnsatzVars::Grid::template Codim<0>::Entity;
  using Position = Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension>;

  typedef typename AnsatzVars::Grid::LeafIndexSet IndexSet;

  template <int row> using TestComponents = std::integral_constant<size_t,TestVars::template Components<row>::m>;
  template <int row> using AnsatzComponents = std::integral_constant<size_t,AnsatzVars::template Components<row>::m>;

  private:
  //PDE indices
  static int const sIdx = 0; // s
  static int const eIdx = 1; // e
  static int const iIdx = 2; // i
  static int const syIdx = 3; // sy
  static int const hIdx = 4; // h
  static int const cIdx = 5; // c
  static int const hcIdx = 6; // hc
  static int const rIdx = 7; // r

  static int const sSIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,sIdx>::type::spaceIndex; 
  static int const eSIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,eIdx>::type::spaceIndex;
  static int const iSIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,iIdx>::type::spaceIndex;
  static int const sySIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,syIdx>::type::spaceIndex; 
  static int const hSIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,hIdx>::type::spaceIndex;
  static int const cSIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,cIdx>::type::spaceIndex;
  static int const hcSIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,hcIdx>::type::spaceIndex; 
  static int const rSIdx = boost::fusion::result_of::value_at_c<typename OriginVars::Variables,rIdx>::type::spaceIndex; 

  public:

  /// \class DomainCache
  ///
  class DomainCache : public CacheBase<AlievPanfilovEquation, DomainCache>
  {
    public:
    DomainCache(Self const& F_,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      F(F_), vars(vars_) 
    {}

    void moveTo(Cell const &entity) { ent = &entity; }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
      xglob = ent->geometry().global(x);

      using namespace boost::fusion;
      
      s   = at_c<sIdx>(vars.data).value(at_c<sSIdx>(evaluators));
      ds  = at_c<sIdx>(vars.data).derivative(at_c<sSIdx>(evaluators))[0];
      
      e   = at_c<eIdx>(vars.data).value(at_c<eSIdx>(evaluators));
      de  = at_c<eIdx>(vars.data).derivative(at_c<eSIdx>(evaluators))[0];
      
      i   = at_c<iIdx>(vars.data).value(at_c<iSIdx>(evaluators));
      di  = at_c<iIdx>(vars.data).derivative(at_c<iSIdx>(evaluators))[0];

      sy   = at_c<syIdx>(vars.data).value(at_c<sySIdx>(evaluators));
      dsy  = at_c<syIdx>(vars.data).derivative(at_c<sySIdx>(evaluators))[0];
      
      h   = at_c<hIdx>(vars.data).value(at_c<hSIdx>(evaluators));
      dh  = at_c<hIdx>(vars.data).derivative(at_c<hSIdx>(evaluators))[0];
      
      c   = at_c<cIdx>(vars.data).value(at_c<cSIdx>(evaluators));
      dc  = at_c<cIdx>(vars.data).derivative(at_c<cSIdx>(evaluators))[0];

      hc   = at_c<hcIdx>(vars.data).value(at_c<hcSIdx>(evaluators));
      dhc  = at_c<hcIdx>(vars.data).derivative(at_c<hcSIdx>(evaluators))[0];

      r   = at_c<rIdx>(vars.data).value(at_c<rSIdx>(evaluators));
      dr  = at_c<rIdx>(vars.data).derivative(at_c<rSIdx>(evaluators))[0];

      // get index of global coordinates 
      double error_min = std::numeric_limits<int>::max();
      double error;
      int index_min;
      xglob = ent->geometry().global(x);
      for (int idx=0; idx<ent->subEntities(2); idx++) {
        error = (ent->geometry().corner(idx)-xglob).two_norm();
        if (error < error_min) {
          error_min = error;
          index_min = idx;
        }
      }
      IndexSet const& indexSet = F.varDesc.indexSet;
      int global_index = indexSet.subIndex(*ent,index_min,2); //global index
      grad_V = F.grad_V[global_index]; 

      t = F.time();
      sigma = F.sigma; 
      gamma = F.gamma;
      eta = F.eta;
      kappa = F.kappa;
      eta_c = F.eta_c;
      phi_i = F.phi_i;
      phi_sy = F.phi_sy;
      phi_h = F.phi_h;
      phi_hc = F.phi_hc;
      
      double act_change_notAtHome = F.activityChangePercentage[int(t)][1]/100.0;
      if (F.no_zero_covid_active()) { // zero or no covid restrictions are currently active
        act_change_notAtHome = -0.91; // approximately 91 percent activity reduction
      }

      beta = (1 + act_change_notAtHome) * F.beta_e_vec[0] * M_PI; 
      if ((1 + act_change_notAtHome) > 1) {
        beta = F.beta_e_vec[0] * M_PI; // in case of ABM: new activities won't be added! so PDE does not get a higher infection rate either
      }
      D = F.D_scale * F.D_vec[0];
    }

    template<int row> 
    Scalar d1_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    { 
      //general PDEs
      if (row==sIdx) return -( (D*ds[0] + s*grad_V[0])*arg.derivative[0][0] + (D*ds[1] + s*grad_V[1])*arg.derivative[0][1]) + 
                          (- s*beta* ((i + sy))) * arg.value[0]; 
      if (row==eIdx) return -( (D*de[0] + e*grad_V[0])*arg.derivative[0][0] + (D*de[1] + e*grad_V[1])*arg.derivative[0][1]) + 
                          (s*beta* ((i + sy) ) - sigma*e) * arg.value[0];
      if (row==iIdx) return -( (D*di[0] + i*grad_V[0])*arg.derivative[0][0] + (D*di[1] + i*grad_V[1])*arg.derivative[0][1]) + 
                          (sigma*e - phi_i*i - gamma*i) * arg.value[0]; 
      if (row==syIdx) return -( (D*dsy[0] + sy*grad_V[0])*arg.derivative[0][0] + (D*dsy[1] + sy*grad_V[1])*arg.derivative[0][1]) + 
                          (gamma*i - phi_sy*sy - eta*sy) * arg.value[0]; 
      if (row==hIdx) return -( (D*dh[0] + h*grad_V[0])*arg.derivative[0][0] + (D*dh[1] + h*grad_V[1])*arg.derivative[0][1]) + 
                          (eta*sy - phi_h*h - kappa*h) * arg.value[0]; 
      if (row==cIdx) return -( (D*dc[0] + c*grad_V[0])*arg.derivative[0][0] + (D*dc[1] + c*grad_V[1])*arg.derivative[0][1]) + 
                          (kappa*h - eta_c*c) * arg.value[0]; 
      if (row==hcIdx) return -( (D*dhc[0] + hc*grad_V[0])*arg.derivative[0][0] + (D*dhc[1] + hc*grad_V[1])*arg.derivative[0][1]) + 
                          (eta_c*c - phi_hc*hc) * arg.value[0]; 
      if (row==rIdx) return -( (D*dr[0] + r*grad_V[0])*arg.derivative[0][0] + (D*dr[1] + r*grad_V[1])*arg.derivative[0][1]) + 
                          (phi_i*i + phi_sy*sy + phi_h*h + phi_hc*hc) * arg.value[0]; 
      assert("wrong index\n"==0);
      return 0;
    }

    template<int row, int col> 
    Scalar d2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      // general PDE
      if (row==sIdx) {
        if (col==eIdx || col==hIdx || col==cIdx || col==hcIdx || col==rIdx) return 0.0;
        if (col==sIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1])
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0]
                              + (- beta* ((i + sy)) ) * arg1.value[0]*arg2.value[0];
        if (col==iIdx)  return + (- beta* s) * arg1.value[0]*arg2.value[0]; 
        if (col==syIdx) return + (- beta* s) * arg1.value[0]*arg2.value[0]; 
      }

      if (row==eIdx) {
        if (col==hIdx || col==cIdx || col==hcIdx || col==rIdx) return 0.0;
        if (col==sIdx) return beta*(i + sy)* arg1.value[0] * arg2.value[0]; 
        if (col==eIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1])
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0]
                              + (- sigma) * arg1.value[0]*arg2.value[0];
        if (col==iIdx)  return (+ (beta* s) * arg1.value[0]) * arg2.value[0]; 
        if (col==syIdx) return (+ (beta* s) * arg1.value[0]) * arg2.value[0]; 
      }

      if (row==iIdx) {
        if (col==sIdx || col==syIdx || col==hIdx || col==cIdx || col==hcIdx || col==rIdx) return 0.0;
        if (col==eIdx) return + sigma*arg1.value[0]*arg2.value[0]; 
        if (col==iIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1]) 
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0]
                              + (- phi_i - gamma) * arg1.value[0]*arg2.value[0]; 
      }

      if (row==syIdx) {
        if (col==sIdx || col==eIdx || col==hIdx || col==cIdx || col==hcIdx || col==rIdx) return 0.0;
        if (col==iIdx) return + gamma*arg1.value[0]*arg2.value[0]; 
        if (col==syIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1]) 
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0]
                              + (- phi_sy - eta) * arg1.value[0]*arg2.value[0]; 
      }

      if (row==hIdx) {
        if (col==sIdx || col==eIdx || col==iIdx || col==cIdx || col==hcIdx || col==rIdx) return 0.0;
        if (col==syIdx) return + eta*arg1.value[0]*arg2.value[0]; 
        if (col==hIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1]) 
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0]
                              + (- phi_h - kappa) * arg1.value[0]*arg2.value[0]; 
      }
      if (row==cIdx) {
        if (col==sIdx || col==eIdx || col==iIdx || col==syIdx || col==hcIdx || col==rIdx) return 0.0;
        if (col==hIdx) return + kappa*arg1.value[0]*arg2.value[0]; 
        if (col==cIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1]) 
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0]
                              + (- eta_c) * arg1.value[0]*arg2.value[0]; 
      }
      if (row==hcIdx) {
        if (col==sIdx || col==eIdx || col==iIdx || col==syIdx || col==hIdx || col==rIdx) return 0.0;
        if (col==cIdx) return + eta_c*arg1.value[0]*arg2.value[0]; 
        if (col==hcIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1]) 
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0]
                              + (- phi_hc) * arg1.value[0]*arg2.value[0]; 
      }

      if (row==rIdx) {
        if (col==sIdx || col==eIdx || col==cIdx) return 0.0; 
        if (col==iIdx) return + phi_i*arg1.value[0]*arg2.value[0]; 
        if (col==syIdx) return + phi_sy*arg1.value[0]*arg2.value[0]; 
        if (col==hIdx) return + phi_h*arg1.value[0]*arg2.value[0]; 
        if (col==hcIdx) return + phi_hc*arg1.value[0]*arg2.value[0]; 
        if (col==rIdx) return - D*(arg1.derivative[0][0]*arg2.derivative[0][0] + arg1.derivative[0][1]*arg2.derivative[0][1]) 
                              + (grad_V[0]*arg1.derivative[0][0] + grad_V[1]*arg1.derivative[0][1]) * arg2.value[0];
      }

  
      assert("wrong index\n"==0);
      return 0.0;
    }

    template<int row, int col> 
    Scalar b2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const 
    {
      if (row==col)  return arg1.value*arg2.value; // b0: 1/2 * u**2, b1: u*arg.value
      else           return 0;
    }

    private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;

    Cell const* ent;

    Scalar s, e, i, sy, h, c, hc, r; 
    Scalar t, beta, D, sigma, gamma, eta, kappa, eta_c, phi_i, phi_sy, phi_h, phi_hc; 

    std::vector<double> grad_V;
    std::vector<double> DensityAgentsInPDE;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> ds, de, di, dsy, dh, dc, dhc, dr; 
    Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension> xglob;
  };

  /// \class BoundaryCache
  ///
  class BoundaryCache : public Kaskade::CacheBase<AlievPanfilovEquation, BoundaryCache>
  {
    public:
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(Self const& F_,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
      F(F_), vars(vars_), ent(0)
    {}

    void moveTo(FaceIterator const& entity)
    {
      ent = &entity;
      penalty = 1.0e9;
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension-1> const& x, Evaluators const& evaluators) 
    {} 
    
    template<int row> 
    Scalar d1_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg) const
    {
      return 0;
    }

    template<int row, int col>  
    Scalar d2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      return 0;
    }
    
    template<int row, int col> 
    Scalar b2_impl (VariationalArg<Scalar,dim,TestComponents<row>::value> const &arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const &arg2) const
    {
      return 0;
    }

    private:
    Self const& F;
    typename AnsatzVars::VariableSet const& vars;
    FaceIterator const* ent;
    Scalar penalty;
  };

  class Scaling 
  {
    public:
    Scaling(Self const& self_): self(&self_) {}

    template <class Vals>
    void scale(Cell const&, Position const&, Vals& vals) const
    {
    }
    template <class Vals>
    void operator()(Cell const& cell, Position const& xi, Vals& vals) const
    {
      scale(cell,xi,vals);
    }

    private:
    Self const* self;
  };

  public:

  bool no_zero_covid_active() const { return zero_covid_bool; }
  Self& setNoZeroCovidBool(bool zero_covid_bool_new) { zero_covid_bool = zero_covid_bool_new; return *this; }

  Scalar time() const { return t; }
  Self& setTime(Scalar tnew)  { t=tnew; return *this; }

  Scaling scaling() const { return Scaling(*this); }

  void temporalEvaluationRange(double t0, double t1) { tau = t1-t0; }

  template <int row>
  struct D1: public FunctionalBase<VariationalFunctional>::D1<row> 
  {
    static bool const present   = true;
    static bool const constant  = false;
  };

  template <int row, int col>
  struct D2: public FunctionalBase<VariationalFunctional>::D2<row,col> 
  {
    static bool const present   = true;
    static bool const symmetric = true;
    static bool const lumped    = false;
  };

  template <int row, int col>
  struct B2 
  {
    static bool const present   = true;
    static bool const symmetric = true;
    static bool const constant  = false;
    static bool const lumped    = false;
  };

  /**
   * Given an initial value, this is transfered to a properly scaled
   * finite element iterate.
   */
  template <int row, class WeakFunctionView>
  void scaleInitialValue(WeakFunctionView const& x0, typename AnsatzVars::VariableSet& x) const 
  {
    interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<row>(x.data),
                                  ScaledFunction<WeakFunctionView>(true,x0,*this));
  }

  /// \fn integrationOrder
  ///

  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder; 
    else
      return 2*shapeFunctionOrder+1; 
  }

  private:
  Scalar t, tau;
  bool zero_covid_bool = false;
     
  template <class WeakFunctionView>
  struct ScaledFunction 
  {
    using Scalar = typename WeakFunctionView::Scalar;
    static int const components = WeakFunctionView::components;
    using ValueType = Dune::FieldVector<Scalar,components>;

    ScaledFunction(bool doScaling_,
                   WeakFunctionView const& s0_,
                   AlievPanfilovEquation<RType,AnsatzVars> const& f_): doScaling(doScaling_), s0(s0_), f(f_) {}

    template <class Cell>
    int order(Cell const&) const { return std::numeric_limits<int>::max(); }

    template <class Cell>
    ValueType value(Cell const& cell,
                    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const 
    {
      ValueType s = s0.value(cell,localCoordinate);
      return s;
    }

    private:
    bool doScaling;
    WeakFunctionView const& s0;
    AlievPanfilovEquation<RType,AnsatzVars> const& f;
  };

  private:  
  Scalar D_scale, r, xi;
  double sigma, gamma, eta, kappa, eta_c, phi_i, phi_sy, phi_h, phi_hc;
  std::vector<double> beta_e_vec;
  std::vector<double> D_vec; 
  std::vector<int> T_vec;
  std::vector<std::vector<double>> activityChangePercentage;
  std::vector<std::vector<double>> grad_V;

  VarSet varDesc;
  double reactionImplicitFactor;
 
  public: 
  AlievPanfilovEquation(double sigma_, double gamma_, double eta_, double kappa_, double eta_c_, double phi_i_, double phi_sy_, double phi_h_, double phi_hc_, 
            std::vector<double> beta_e_vec_, std::vector<std::vector<double>> activityChangePercentage_, 
            std::vector<double> D_vec_, Scalar D_scale_, 
            std::vector<std::vector<double>> grad_V_, std::vector<int> T_vec_, 
            Scalar xi_, VarSet varDesc_):
            t(0), sigma(sigma_), gamma(gamma_), eta(eta_), kappa(kappa_), eta_c(eta_c_), phi_i(phi_i_), phi_sy(phi_sy_), phi_h(phi_h_), phi_hc(phi_hc_), 
            beta_e_vec(beta_e_vec_), activityChangePercentage(activityChangePercentage_), D_vec(D_vec_), D_scale(D_scale_), grad_V(grad_V_), T_vec(T_vec_), xi(xi_), 
            varDesc(varDesc_),
            reactionImplicitFactor(1) 
  {}
};

#endif