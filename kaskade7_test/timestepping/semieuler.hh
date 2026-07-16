/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2021 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SEMIEULERINNER_HH
#define SEMIEULERINNER_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"
#include "fem/functional_aux.hh"
namespace Kaskade
{
  namespace SemiEulerDetail
  {
    template <class ParabolicEquation, class SemiEuler, bool hasInnerBoundaryCache> 
    struct SemiImplicitEulerInnerBoundaryPolicyImpl
    {
    };

    template <class ParabolicEquation, class SemiEuler>
    struct SemiImplicitEulerInnerBoundaryPolicyImpl<ParabolicEquation,SemiEuler,true>
    {
    private:
      using SemiEulerStep = SemiEuler;
      
    public:
      /**
       * \brief Evaluates jump contributions at interior faces, suitable for discontinuous Galerkin methods (optional).
       * 
       * If this structure is defined, the assembler does assemble jump terms at faces of the grid. If undefined, it does not.
       */
      class InnerBoundaryCache
      {
        using AnsatzVars = typename SemiEulerStep::AnsatzVars;
        using TestVars = typename SemiEulerStep::TestVars;
        using AnsatzVariableSet = typename AnsatzVars::VariableSet;
        using Scalar = typename ParabolicEquation::Scalar;
        
      public:
        InnerBoundaryCache(SemiEulerStep const& f_,
                           AnsatzVariableSet const& vars_,
                           AnsatzVariableSet const& varsJ_,
                           AnsatzVariableSet const& du_,
                           int flags)
      : f(f_)
      , peibc(f.parabolicEquation(),vars_,flags)
      , peibcJ(f.parabolicEquation(),varsJ_,flags)
      , du(du_)  
      {}

        
        template <class FaceIterator>
        void moveTo(FaceIterator const& face)
        {
        peibc.moveTo(face);
        peibcJ.moveTo(face);
        }
        
        /**
         * \brief Prepare the cache for evaluation at given quadrature point.
         * 
         * \param x the quadrature point in the local face coordinate system
         * \param evaluators evaluators for quantities on the current cell 
         * \param neighbourEvaluators evaluators for the cell on the opposite side of the face
         */
        template <class Position, class Evaluators>
        void evaluateAt(Position const& x, Evaluators const& evaluators, Evaluators const& neighbourEvaluators)
        {
          peibc.evaluateAt(x, evaluators, neighbourEvaluators);
          peibcJ.evaluateAt(x, evaluators,neighbourEvaluators);
        }
        
        
        /**
         * Returns the directional derivative of \f$g\f$ at \f$ 0 \f$,
         * with respect to the variable \f$u_{\mathrm{row}}\f$ in direction of the test
         * function given by arg.
         * 
         * The test variational argument belongs to the current cell, i.e. the interior side of the face.
         */
        template <int row, int dim>
        Dune::FieldVector<Scalar,TestVars::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& argT) const
        {
          return f.getTau() * peibc.template d1<row>(argT);
        }
        
        /**
         * Returns the second directional derivative of 
         * \f$ g \f$ at \f$ 0 \f$, with
         * respect to the variables \f$u_{\mathrm{row}}\f$ and
         * \f$ u_{\mathrm{col}}\f$ in direction of the test functions given
         * by argT and the ansatz functions, given by argA.
         * 
         * If centerCell is true, both test and ansatz variational arguments are associated on the current cell.
         * Otherwise, the ansatz argument argA belongs to the neighbour cell, i.e. the outside of the face.
         *
         * This method need only be defined for values of row/col for which D2<row,col>::present is true.
         */
        template <int row, int col, int dim>
        Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m>
        d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA, bool centerCell) const
        {
          Dune::FieldMatrix<Scalar,TestVars::template components<row>,
                                   AnsatzVars::template components<col>> result(0);
          
          
          if (ParabolicEquation::template B2<row,col>::present)
            result = peibcJ.template b2<row,col>(argT, argA, centerCell);
          
          if (ParabolicEquation::template D2<row,col>::present)
            result -= f.getTau()*peibcJ.template d2<row,col>(argT, argA, centerCell);
          
          return result;
        }
        

      private:
        SemiEulerStep const&                           f;
        typename ParabolicEquation::InnerBoundaryCache peibc;
        typename ParabolicEquation::InnerBoundaryCache peibcJ;
        AnsatzVariableSet const&                       du;
      };

      template <class Face>
      bool considerFace(Face const& face) const
      {
        return static_cast<SemiEuler const*>(this)->parabolicEquation().considerFace(face);
      }

    };
  }
  
  
  /**
   * \ingroup timestepping
   * \brief Linearly implicit Euler method.
   *
   * A class that defines the weak formulation of the stationary
   * elliptic problem that results from the linearly implicit Euler
   * method for
   * \f[ B \dot u  = L(u) , \f]
   * which is, doing the linearization at \f$ u_0 \f$,  \f[
   * (B - \tau L'(u_0)) \delta u_k  = \tau  L(u_k), \quad u_{k+1} = u_k + \delta u_k. \f]
   * The left hand side \f$ B - \tau L'(u_0) \f$ is assembled as matrix, and the right hand
   * side \f$ \tau  L(u_k) \f$ is assembled as right hand side by the assembler, if fed with
   * this weak problem formulation:
   * \code
   * assembler.assemble(semiLinearization(SemiImplicitEulerStep,u,uJ,du));
   * \endcode
   *
   * \f$ B \f$ may depend on \f$ u \f$, made explicit by setting \code
   * B2<row,col>::constant = false \endcode in the parabolic equation
   * class. In this case, \f$ B(u) \f$ has to be an invertable Nemyckii operator (a
   * diagonal operator, mapping \f$ u(x) \f$ to \f$ B(x,u(x)) u(x)
   * \f$. Then the linearly implicit Euler scheme is (see below for derivation)
   * \f[
   * (B(u_0) - \tau L'(u_0) + \tau_k B^{-1}(u_0)B'(u_0)\mathrm{diag}L(u_0)) \delta u_k  
   * = \tau_k  L(u)+\frac{\tau_k}{\tau_{k-1}}(B(u_0)-B(u_k))\delta u_{k-1}.
   * \f]
   * The values \f$ u_0 \f$ (linearization point for the Jacobian), \f$ u_k \f$ (evaluation point 
   * for the right hand side), and \f$ \tau_k/\tau_{k-1}\, \delta u_{k-1}  \f$ (previous increment)
   * have to be specified when constructing the DomainCache.
   * This is most conveniently done via \ref SemiLinearizationAt.
   * 
   * \tparam PE the weak formulation of the parabolic equation, adhering to the ParabolicEquation interface.
   * 
   * \em Derivation: Lubich and Roche (Rosenbrock Methods for Differential-algebraic Systems with 
   * Solution-dependent Singular Matrix Multiplying the Derivative, Computing 43:325-342, 1990) suggest
   * to reformulate the parabolic system with non-constant \f$ B \f$ above as
   * \f[ \begin{aligned} \dot u &= z \\ 0 &= L(u)-B(u)z \end{aligned} \f]
   * and applying, e.g., Rosenbrock methods. Here, simply a linearly implicit Euler. Linearizing
   * the right hand side of the extended system at \f$ u_0, z_0 \f$ yields
   * \f[ \Bigg( \begin{bmatrix} I & 0 \\ 0 & 0 \end{bmatrix}  
   *     - \tau_k \begin{bmatrix} 0 & I \\ L'(u_0)-B'(u_0)z_0 & -B(u_0) \end{bmatrix} \Bigg)
   *     \begin{bmatrix} \delta u_k \\ \delta z_k \end{bmatrix}
   *     = \tau_k \begin{bmatrix} z_k \\ L(u_k)-B(u_k)z_k \end{bmatrix}. \f]
   * The first row immediately gives \f$ \delta u_k = \tau_k (z_k + \delta z_k) \f$, and the second 
   * one 
   * \f[ - \tau_k (L'(u_0)-B'(u_0)z_0) \delta u_k + \tau (B(u_0)\delta z_k + B(u_k)z_k) = \tau_k f(u_k). \f]
   * The second term on the left hand side is just 
   * \f[ \tau_k (B(u_0)(z_k+\delta z_k) + (B(u_k)-B(u_0))z_k) = B(u_0)\delta u_k +  (B(u_k)-B(u_0)) \tau_k z_k, \f]
   * hence we obtain
   * \f[ \left(B(u_0) - \tau_k(L'(u_0)-B'(u_0)z_0)\right)\delta u_k = \tau_k L(u_k)+(B(u_0)-B(u_k)) \frac{\tau_k}{\tau_{k-1}}\delta u_{k-1}. \f]
   * Choosing a consistent initial value \f$ z_0 = B(u_0)^{-1} L(u_0) \f$ finally results in the form given above.
   */
  template <class PE>
  class SemiImplicitEulerStep
  : public SemiEulerDetail::SemiImplicitEulerInnerBoundaryPolicyImpl<PE,SemiImplicitEulerStep<PE>,
                                                                     hasInnerBoundaryCache<PE>>   
  {
    using Self = SemiImplicitEulerStep<PE>;
    using Base = SemiEulerDetail::SemiImplicitEulerInnerBoundaryPolicyImpl<PE,Self,hasInnerBoundaryCache<PE>>;

  public:
    typedef PE ParabolicEquation;

    typedef typename ParabolicEquation::Scalar  Scalar;
    typedef typename ParabolicEquation::AnsatzVars AnsatzVars;
    typedef typename ParabolicEquation::TestVars TestVars;
    typedef typename ParabolicEquation::OriginVars OriginVars;
    typedef typename AnsatzVars::Grid Grid;
    using GridView = typename AnsatzVars::GridView;

    static ProblemType const type = WeakFormulation;

    typedef PE  EvolutionEquation;

    typedef typename EvolutionEquation::AnsatzVars::template CoefficientVectorRepresentation<0,EvolutionEquation::TestVars::noOfVariables>::type CoefficientVectors;

    
    
    class DomainCache
    {
    public:
      /**
       * \brief Constructor.
       * \param f the weak formulation
       * \param vars evaluation point \f$ u_k \f$ for the right hand side
       * \param varsJ linearization point \f$ u_0 \f$ for the (extended) Jacobian
       * \param duVars previous increment \f$ \delta u_{k-1} \f$ (is irrelevant if
       *               \f$ B \f$ does not depend on \f$ u \f$, i.e. B2::constant is true)
       * \param tauRatio ratio of current and previous time step \f$ \tau_k / \tau_{k-1} \f$
       */
      DomainCache(SemiImplicitEulerStep<ParabolicEquation> const& f_, // TODO: implement tauRatio
                  typename AnsatzVars::VariableSet const& vars_,
                  typename AnsatzVars::VariableSet const& varsJ_,
                  typename AnsatzVars::VariableSet const& duVars_,
                  int flags):
        f(f_), pedc(*f_.eq,vars_,flags), pedcJ(*f_.eq,varsJ_,flags), duVars(duVars_) 
      {} 

      void moveTo(typename Grid::template Codim<0>::Entity const& entity) {
        pedc.moveTo(entity);
        pedcJ.moveTo(entity);
      }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators) 
      {
        pedc.evaluateAt(x, evaluators);
        pedcJ.evaluateAt(x, evaluators);
        du = evaluateVariables(duVars,evaluators,valueMethod);
      }

      template<int row, int dim>
      Dune::FieldVector<Scalar, TestVars::template components<row>>
      d1 (VariationalArg<Scalar,dim> const& arg) const 
      {
        Dune::FieldVector<Scalar, TestVars::template components<row>> result( f.getTau() * pedc.template d1<row>(arg) );

        boost::fusion::for_each(typename AnsatzVars::Variables(),MultiplyDifference<row,dim>(pedc,pedcJ,du,arg,result));
        return result;
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar, TestVars::template components<row>, AnsatzVars::template components<col>>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const
      {
        Dune::FieldMatrix<Scalar,TestVars::template components<row>,AnsatzVars::template components<col>> result(0);


        if (ParabolicEquation::template B2<row,col>::present)
          result = pedcJ.template b2<row,col>(arg1, arg2);

        if (ParabolicEquation::template D2<row,col>::present)
          result -= f.getTau()*pedcJ.template d2<row,col>(arg1, arg2);
        
        return result;
      }

    private:
      SemiImplicitEulerStep<ParabolicEquation> const& f;
      typename ParabolicEquation::DomainCache pedc;  // for evaluating rhs at u_k
      typename ParabolicEquation::DomainCache pedcJ; // for evaluating Jacobian at u_0
      typename AnsatzVars::VariableSet const& duVars;
      EvaluateVariables<typename AnsatzVars::VariableSet,ValueMethod> du;


      template <int row, int dim>
      class MultiplyDifference
      {
      public:
        MultiplyDifference(typename ParabolicEquation::DomainCache pedc_,
                           typename ParabolicEquation::DomainCache pedcJ_,
                           EvaluateVariables<typename AnsatzVars::VariableSet,ValueMethod> const& du_,
                           VariationalArg<Scalar,dim> const& arg_,
                           Dune::FieldVector<Scalar, TestVars::template Components<row>::m>& result_)
        : pedc(pedc_), pedcJ(pedcJ_), du(du_), arg(arg_), result(result_)
        {}

        template <class VarDescription>
        void operator()(VarDescription const& vd) const 
        {
          // compute <(B(u0)-B(uk))*du,arg>
          int const col = VarDescription::id;

          if (ParabolicEquation::template B2<row,col>::present && !ParabolicEquation::template B2<row,col>::constant) 
          {
            // @TODO: should work with vector-valued shape functions
            VariationalArg<Scalar,dim> duArg;
            duArg.value[0] = 1;
            result += (pedcJ.template b2<row,col>(arg,duArg)-pedc.template b2<row,col>(arg,duArg)) * boost::fusion::at_c<col>(du);
          }
        }

      private:
        typename ParabolicEquation::DomainCache const& pedc;
        typename ParabolicEquation::DomainCache const& pedcJ;
        EvaluateVariables<typename AnsatzVars::VariableSet,ValueMethod> const& du;
        VariationalArg<Scalar,dim> const& arg;
        Dune::FieldVector<Scalar, TestVars::template Components<row>::m>& result;
      };

    };

    // TODO: currently  a dependence of the boundary part of B on u is not implemented
    class BoundaryCache
    {
    public:
      typedef typename AnsatzVars::GridView::IntersectionIterator FaceIterator;

      BoundaryCache(SemiImplicitEulerStep<ParabolicEquation> const& f_,
                    typename AnsatzVars::VariableSet const& vars_,
                    typename AnsatzVars::VariableSet const& varsJ_,
                    typename AnsatzVars::VariableSet const& du_,
                    int flags)
      : f(f_), pedc(*f_.eq,vars_,flags), pedcJ(*f_.eq,varsJ_,flags), du(du_)  {}

      void moveTo(FaceIterator const& entity) 
      {
        pedc.moveTo(entity);
        pedcJ.moveTo(entity);
      }

      template <class Evaluators>
      void evaluateAt(Dune::FieldVector<typename Grid::ctype,Grid::dimension-1> const& x,
                      Evaluators const& evaluators) 
      {
        pedc.evaluateAt(x, evaluators);
        pedcJ.evaluateAt(x, evaluators);
      }

      template<int row, int dim>
      Dune::FieldVector<Scalar, TestVars::template components<row>> d1 (VariationalArg<Scalar,dim> const& arg) const
      {
        return f.getTau() * pedc.template d1<row>(arg);
      }

      template<int row, int col, int dim>
      Dune::FieldMatrix<Scalar, TestVars::template components<row>, AnsatzVars::template components<col>>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
      {
        // to be corrected - include b2
        if (ParabolicEquation::template D2<row,col>::present)
          return - f.getTau()*pedcJ.template d2<row,col>(arg1, arg2);

        return 0.0; // never get here
      }

    private:
      SemiImplicitEulerStep<ParabolicEquation> const& f;
      typename ParabolicEquation::BoundaryCache pedc;
      typename ParabolicEquation::BoundaryCache pedcJ;
      typename AnsatzVars::VariableSet const& du;
    };



  public:
    template <int row, int col>
    struct D2: public ParabolicEquation::template D2<row,col>
    {
      static bool const present = (ParabolicEquation::template D2<row,col>::present ||
                                   ParabolicEquation::template B2<row,col>::present);
      static bool const symmetric = (!ParabolicEquation::template D2<row,col>::present ||
                                     ParabolicEquation::template D2<row,col>::symmetric) &&
                                    (!ParabolicEquation::template B2<row,col>::present ||
                                     ParabolicEquation::template B2<row,col>::symmetric);
      //     static bool const constant = (!ParabolicEquation::template D2<row,col>::present ||
      //                                   ParabolicEquation::template D2<row,col>::constant) &&
      //                                  (!ParabolicEquation::template B2<row,col>::present ||
      //                                   ParabolicEquation::template B2<row,col>::constant);
      static bool const lumped = ((ParabolicEquation::template D2<row,col>::lumped ||
          !ParabolicEquation::template D2<row,col>::present) &&
          (ParabolicEquation::template B2<row,col>::lumped||
              !ParabolicEquation::template B2<row,col>::present));

    };

    template <int row>
    struct D1: public ParabolicEquation::template D1<row>
    {
    };

    /**
     * \brief Constructor.
     * 
     * \param eq a pointer to a parabolic equation object. This needs to exist during the life
     *           time of the SemiImplicitEulerStep object.
     * \param dt the time step. dt>=0 has to hold.
     */
    SemiImplicitEulerStep(ParabolicEquation *eq_, Scalar dt) 
    : eq(eq_)
    , tau(dt)
    {}

    template <class Cell>
    int integrationOrder(Cell const& cell, int shapeFunctionOrder, bool boundary) const
    {
      return eq->integrationOrder(cell,shapeFunctionOrder,boundary);
    }

    /**
     * \brief Sets the time step size.
     * 
     * \param dt The new time step size. dt>=0 has to hold.
     */
    Self& setTau(Scalar dt)
    { 
      tau = dt; 
      assert(tau >= 0);
      return *this;
    }

    Scalar getTau() const { return tau; }


    /**
     * \name Access to the underlying parabolic equation.
     * @{
     */
    ParabolicEquation&       parabolicEquation()       { return *eq; }
    ParabolicEquation const& parabolicEquation() const { return *eq; }
    /**
     * @}
     */

    


    private:
      ParabolicEquation *eq;
      Scalar tau;                        // the time step
      bool onlyJacobian = false;         // 
      bool onlyMassMatrix = false;       // 
  };
} // namespace Kaskade
#endif
