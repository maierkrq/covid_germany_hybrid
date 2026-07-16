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

#ifndef FUNCTIONAL_AUX_HH
#define FUNCTIONAL_AUX_HH

#include <type_traits>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include "fem/functionspace.hh"
/**
 * @file
 * @brief  Utility classes for the definition and use of variational functionals
 * @author Anton Schiela, Lars Lubkoll
 */

namespace Kaskade
{
  /**
   * \ingroup functional
   *
   * \brief Convenience base class for the easy definition of variational
   * functionals and weak formulations.
   *
   * Variational functionals and weak formulations need to provide meta
   * information template classes D1 and D2 that give the assembler
   * static information about the block structure of the problem. It is
   * easy to forget parts of these meta data, in particular if they do
   * not affect your current problem (yet they need to be defined for
   * syntactical reasons).
   *
   * Problem definitions may inherit this convenience base class, which
   * defines safe but maybe inefficient defaults. These can be partially
   * overridden as in the following example.
   *
   * \code
   * class Functional: public FunctionalBase<VariationalFunctional> {
   *   ...
   *   template <int row, int col>
   *   struct D2: public FunctionalBase<VariationalFunctional>::template D2<row,col> {
   *     static bool symmetric = true;
   *   };
   *   ...
   * };
   * \endcode
   * 
   * \seealso VariationalFunctional
   * 
   * \tparam Type the kind of variational problem, either \ref VariationalFunctional or \ref WeakFormulation
   */
  template <ProblemType Type>
  class FunctionalBase
  {
  public:

    static constexpr ProblemType type = Type;

    
    // TODO: docme
    template <class T>
    bool inDomain(T const&) const
    {
      return true;
    }
    
    /**
     * \brief Static a priori default information about the right hand side blocks.
     */
    template <int row>
    struct D1
    {
      /** 
       * \brief Whether the block is to be assembled at all.
       */
      static constexpr bool present = true;
      
      /**
       * \brief The maximum ansatz and test function derivative occuring in this block.
       * 
       * A default value of 1 is provided in the base class, which is appropriate for second order PDEs.
       */
      static constexpr int derivatives = 1;
    };

    /**
     * \brief Static a priori default information about the matrix blocks.
     * 
     * These flags declare structural properties of the individual blocks arising in the 
     * matrix of this problem and can be exploited by the assembler implementation for 
     * improved efficiency.
     */
    template <int row, int col>
    struct D2
    {
      /**
       * \brief Whether the block is to be assembled and stored or not.
       * The default value is the conservative one (but maybe inefficient) that all blocks are structurally nonzero. 
       * Note that for variational problems (which are necessarily symmetric), only the lower triangular part 
       * is accessed and hence in this case the default omits the superdiagonal blocks.
       * 
       * If present is false, the functional's d2 method will never be called for this block.
       */
      static constexpr bool present = type==WeakFormulation || row>=col;
      /**
       * \brief Whether the block is structurally symmetric (or hermitian).
       * The default is true for diagonal blocks of variational problems (which are necessarily symmetric), and
       * the conservative non-symmetric otherwise. Note that actually symmetric blocks are correctly stored if
       * this flag is false, but the storage may be less efficient. The assembler implementation may choose
       * whichever storage it likes, this flag is merely a hint that symmetry-exploiting optimizations are safe.
       */
      static constexpr bool symmetric = type==VariationalFunctional && row==col;
      
      /**
       * \brief If this flag is true (and symmetric==true), the assembler will enforce positive
       *        semidefiniteness of all local matrices by projecting them onto the cone of psd 
       *        matrices.
       * 
       * While this destroys accuracy in the Taylor approximation, it can be very convenient in 
       * algorithms for nonconvex minimization.
       */
      static constexpr bool makePositive = false;

      /**
       * \brief If this flag is true, only the diagonal of this block will be assembled and stored.
       * This is mainly useful for hierarchical error estimators. Note that setting this flag true 
       * for mass matrices does not give a traditional lumped mass matrix (as the integration weights
       * are different).
       */
      static constexpr bool lumped = false;

      static constexpr bool constant = false;

      /**
       * \brief The maximum ansatz and test function derivative occuring in this block.
       * 
       * A default value of 1 is provided in the base class, which is appropriate for second order PDEs.
       */
      static constexpr int derivatives = 1;
    };

    template <int row, int col>
    struct B2
    {
      static bool const present = true;
      static bool const symmetric = false;
      static bool const lumped = false;
      static bool const constant = false;
    };
    
    /**
     * \brief Default implementation specifying the integration order.
     * 
     * This is a conservative choice for linear problems, with an integration order of twice the shape function order.
     */
    template <class Cell>
    int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool /* boundary */) const 
    {
      return 2*shapeFunctionOrder;
    }
    
    /**
     * \brief Default implementation for which faces are to be considered for assembly in case an
     *        InnerBoundaryCache is present.
     *
     * The default implementation is conservative, i.e. it assumes all faces should be considered.
     * Redefine this in the actual implementation of the variational functional if only some faces
     * need to be assembled.
     *
     * \return true
     */
    template <class Face>
    bool considerFace(Face const& ) const
    {
      return true;
    }

    /**
     * \brief Defines the eigenvalue threshold when enforcing positivity of local stiffness matrices.
     *
     * This default implementation defines this in principle as zero (i.e. aiming at positive semidefiniteness
     * in its usual ldefinition), but accounts for some roundoff error during eigenvalue computation
     * assuming assembly in double precision.
     *
     * This method needs to be defined if `D2<,>::makePositive == true`.
     */
    template <int /*row*/, int /*col*/>
    double makePositiveThreshold() const
    {
      return -10*std::numeric_limits<double>::epsilon();
    }
  };

  // ----------------------------------------------------------------------------------------------

  namespace Functional_Aux_Detail
  {
    template <class Functional, int n>
    struct checkSymmetry
    {
      static constexpr bool ok = (!Functional::template D2<n,n>::present)
                                 || Functional::template D2<n,n>::symmetric;
      static constexpr bool pass = Functional::type != VariationalFunctional
                                   || (ok && checkSymmetry<Functional,n-1>::pass);
    };

    template <class Functional>
    struct checkSymmetry<Functional,-1>
    {
      static constexpr bool pass = true;
    };
  }

  template <class Functional, int n>
  constexpr bool symmetryCheck = Functional_Aux_Detail::checkSymmetry<Functional,n>::pass;


  // ----------------------------------------------------------------------------------------------


  /// \cond internals
  namespace Functional_Aux_Detail
  {
    template <int a, int b>
    struct LessThan
    {
      static constexpr bool value = a<b;
    };

    template <int a, int b>
    struct Equals
    {
      static constexpr bool value = a==b;
    };
  }
  /// \endcond

  /**
   * \ingroup functional
   * \brief Base class simplifying the definition of domain and boundary caches.
   * 
   * It provides do-nothing default implementations of moveTo and evaluateAt, and 
   * zero default implementations for d0, d1, and d2.
   */
  template <class AnsatzVars, class TestVars>
  class TrivialCacheBase
  {
  public:
    using Scalar = typename AnsatzVars::Scalar;
    
    template <class Entity>
    void moveTo(Entity const&) {}
    
    template <class Position, class Evaluators>
    void evaluateAt(Position const&, Evaluators const&) {}
    
    Scalar d0()
    {
      return 0;
    }
    
    template <int row, int dim>
    Dune::FieldVector<Scalar,TestVars::template components<row>> d1(VariationalArg<Scalar,dim> const& arg) const
    {
      return 0;
    }
    
    template <int row, int col, int dim>
    Dune::FieldMatrix <Scalar,TestVars::template components<row>,AnsatzVars::template components<row>>
    d2(VariationalArg<Scalar,dim> const& , VariationalArg<Scalar,dim> const& ) const
    {
      return 0;
    }
    
  };

  /**
   * \ingroup functional
   * \brief Provides implementations of d1 and d2 for (possibly mixed) use of scalar and power spaces.
   *
   * Implemented using CRTP. Just derive publicly from this functional if you want to use it.
   * Uses VariationalArg as perturbation. 
   * 
   * \todo docme better
   * 
   * \tparam Functional actual functional (only transports static information)
   * \tparam Cache Functional::DomainCache or Functional::BoundaryCache, must implement
   * Scalar d1_impl(...);
   * and
   * Scalar d2_impl(...);
   */
  template <class Functional, class Cache, class ScalarT = typename Functional::Scalar>
  class CacheBase: public TrivialCacheBase<typename Functional::AnsatzVars,typename Functional::TestVars>
  {
  public:
    using Scalar = ScalarT;
    using AnsatzVars = typename Functional::AnsatzVars;
    using TestVars = typename Functional::TestVars;

  private:
    template <int row>          using Vector = Dune::FieldVector<Scalar, TestVars::template Components<row>::m>;
    template <int row, int col> using Matrix = Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>;
    template <int row>          using TestComponents =   std::integral_constant<int,TestVars::template Components<row>::m>;
    template <int row>          using AnsatzComponents = std::integral_constant<int,AnsatzVars::template Components<row>::m>;
    
  public:
    /// forward to implementation via CRTP
    template <int row, int dimw>
    Scalar d1_local(VariationalArg<Scalar,dimw,TestComponents<row>::value> const& arg) const
    {
      return static_cast<Cache const*>(this)->template d1_impl<row>(arg);
    }

    template <int row, int dimw>
    Dune::FieldVector<Scalar,1> d1(VariationalArg<Scalar,dimw,TestComponents<row>::value> const& arg) const
    {
      return Dune::FieldVector<Scalar,1>(d1_local<row>(arg));
    }  

    /// scalar argument for multi-component variables
    template <int row, int dimw, class enable = typename std::enable_if_t<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>>
    Vector<row> d1(VariationalArg<Scalar,dimw> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dimw,TestComponents<row>::value> arg;

      for(int i=0; i<TestComponents<row>::value; ++i)
      {
        arg.value = 0; arg.value[i] = scalarArg.value[0];
        arg.derivative = 0; arg.derivative[i] = scalarArg.derivative[0];

        result[i] = d1_local<row>(arg);
      }

      return result;
    }

    /// forward to implementation via CRTP
    template <int row, int col, int dim>
    Scalar d2_local(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const 
    {
      return static_cast<Cache const*>(this)->template d2_impl<row,col>(arg1,arg2);
    }

    template <int row, int col, int dim,
    class enable1 = typename std::enable_if_t<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>,
    class enable2 = typename std::enable_if_t<Functional_Aux_Detail::LessThan<1,AnsatzComponents<row>::value>::value>
    >
    Dune::FieldMatrix<Scalar,1,1> d2(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const 
    {
      return d2_local<row,col>(arg1,arg2);
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m2,
              class enable1 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
              class enable2 = typename std::enable_if<Functional_Aux_Detail::Equals<m2,AnsatzComponents<col>::value>::value>::type,
              class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
    >
    Vector<row> d2(VariationalArg<Scalar,dim> const& scalarArg, VariationalArg<Scalar,dim,m2> const& arg2) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg.derivative[0];

        result[i] = d2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m1,
              class enable1 = typename std::enable_if<Functional_Aux_Detail::Equals<m1,TestComponents<row>::value>::value>::type,
              class enable2 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
              class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
             >
    Vector<row> d2(VariationalArg<Scalar,dim,m1> const& arg1, VariationalArg<Scalar,dim> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,AnsatzComponents<row>::value> arg2;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg2.value = 0; arg2.value[i] = scalarArg.value[0];
        arg2.derivative = 0; arg2.derivative[i] = scalarArg.derivative[0];

        result[i] = d2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim>
    Matrix<row,col> d2(VariationalArg<Scalar,dim> const& scalarArg1, VariationalArg<Scalar,dim> const& scalarArg2) const
    {
      Matrix<row,col> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      VariationalArg<Scalar,dim,AnsatzComponents<col>::value> arg2;
      for(int i=0; i<TestComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg1.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg1.derivative[0];
        for(int j=0; j<AnsatzComponents<col>::value; ++j)
        {
          arg2.value = 0; arg2.value[j] = scalarArg2.value[0];
          arg2.derivative = 0; arg2.derivative[j] = scalarArg2.derivative[0];

          result[i][j] = d2_local<row,col>(arg1,arg2);
        }
      }

      return result;
    }


    /// forward to implementation via CRTP
    template <int row, int col, int dim>
    Scalar b2_local(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const
    {
      return static_cast<Cache const*>(this)->template b2_impl<row,col>(arg1,arg2);
    }

    template <int row, int col, int dim,
    class enable1 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
    class enable2 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<row>::value>::value>::type
    >
    Dune::FieldMatrix<Scalar,1,1> b2(VariationalArg<Scalar,dim,TestComponents<row>::value> const& arg1, VariationalArg<Scalar,dim,AnsatzComponents<col>::value> const& arg2) const
    {
      return b2_local<row,col>(arg1,arg2);
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m2,
    class enable1 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
    class enable2 = typename std::enable_if<Functional_Aux_Detail::Equals<m2,AnsatzComponents<col>::value>::value>::type,
    class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
    >
    Vector<row> b2(VariationalArg<Scalar,dim> const& scalarArg, VariationalArg<Scalar,dim,m2> const& arg2) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg.derivative[0];

        result[i] = b2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim, int m1,
    class enable1 = typename std::enable_if<Functional_Aux_Detail::Equals<m1,TestComponents<row>::value>::value>::type,
    class enable2 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,TestComponents<row>::value>::value>::type,
    class enable3 = typename std::enable_if<Functional_Aux_Detail::LessThan<1,AnsatzComponents<col>::value>::value>::type
    >
    Vector<row> b2(VariationalArg<Scalar,dim,m1> const& arg1, VariationalArg<Scalar,dim> const& scalarArg) const
    {
      Vector<row> result(0);
      VariationalArg<Scalar,dim,AnsatzComponents<row>::value> arg2;
      for(int i=0; i<AnsatzComponents<row>::value; ++i)
      {
        arg2.value = 0; arg2.value[i] = scalarArg.value[0];
        arg2.derivative = 0; arg2.derivative[i] = scalarArg.derivative[0];

        result[i] = b2_local<row,col>(arg1,arg2);
      }

      return result;
    }

    /// scalar argument for multi-component variables
    template <int row, int col, int dim>
    Matrix<row,col> b2(VariationalArg<Scalar,dim> const& scalarArg1, VariationalArg<Scalar,dim> const& scalarArg2) const
    {
      Matrix<row,col> result(0);
      VariationalArg<Scalar,dim,TestComponents<row>::value> arg1;
      VariationalArg<Scalar,dim,AnsatzComponents<col>::value> arg2;
      for(int i=0; i<TestComponents<row>::value; ++i)
      {
        arg1.value = 0; arg1.value[i] = scalarArg1.value[0];
        arg1.derivative = 0; arg1.derivative[i] = scalarArg1.derivative[0];
        for(int j=0; j<AnsatzComponents<col>::value; ++j)
        {
          arg2.value = 0; arg2.value[j] = scalarArg2.value[0];
          arg2.derivative = 0; arg2.derivative[j] = scalarArg2.derivative[0];

          result[i][j] = b2_local<row,col>(arg1,arg2);
        }
      }

      return result;
    }
  };
  
  // ----------------------------------------------------------------------------

  /// \cond internals
  namespace Functional_Aux_Detail
  {
    // helper function that detects if InnerBoundaryCache is declared in Functional
    // If it exists, calling this with a nullptr argument selects this more specific template,
    // otherwise the less specific version below.
    template <class Functional>
    constexpr std::true_type hasInnerBoundaryCacheImpl(typename Functional::InnerBoundaryCache*);

    /// helper function that detects if InnerBoundaryCache is not declared in Functional
    template <class Functional>
    constexpr std::false_type hasInnerBoundaryCacheImpl(...);
  }
  /// \endcond

  /** 
   * \ingroup functional 
   * \brief Checks whether InnerBoundaryCache is declared in Functional
   */
  template <class Functional>
  constexpr bool hasInnerBoundaryCache =
  decltype(Functional_Aux_Detail::hasInnerBoundaryCacheImpl<Functional>(nullptr))::value;

  

  // ----------------------------------------------------------------------------

  /// \cond internals
  namespace Functional_Aux_Detail
  {
    template <class Functional, class Linearization, bool hasInnerBoundaryCache> 
    struct InnerBoundaryPolicyImpl
    {
    };

    template <class Functional, class Linearization>
    struct InnerBoundaryPolicyImpl<Functional,Linearization,true>
    {
      using InnerBoundaryCache = typename Functional::InnerBoundaryCache;

      /// create inner boundary cache if existent
      InnerBoundaryCache createInnerBoundaryCache(int flags) const 
      { 
        return InnerBoundaryCache(static_cast<Linearization const*>(this)->getFunctional(),
                                  static_cast<Linearization const*>(this)->getOrigin(),flags);
      }

      template <class Face>
      bool considerFace(Face const& face) const
      {
        return static_cast<Linearization const*>(this)->getFunctional().considerFace(face);
      }
    };

    template <class Functional,class Linearization> 
    using InnerBoundaryPolicy = InnerBoundaryPolicyImpl<Functional,Linearization,hasInnerBoundaryCache<Functional>>;
  }
  /// \endcond

  /** 
   * \ingroup functional
   * \brief Proxy class for the linearization of a nonlinear functional. 
   * 
   * This adaptor maps the \ref NonlinearVariationalFunctional concept to the \ref VariationalFunctional concept.
   * 
   * According to the following idea: A
   * VariationalFunctional defines a mapping \f$ F : U \to V \f$, which
   * can be evaluated and linearized in an appropriate ansatz space
   * \f$ U_a\subset U \f$, only if an argument \f$ u \in U_o \f$ in a -
   * possibly different - "origin" subspace \f$ U_o \subset U \f$ is
   * given.
   *
   * Hence, from a functional \f$ f \f$ and a variable \f$ u\in U_o \f$ one has to
   * create a linearization, which is represented by this class.
   */
  template<class Functional_>
  class LinearizationAt : public Functional_Aux_Detail::InnerBoundaryPolicy<Functional_,LinearizationAt<Functional_>>
  {
  public:
    typedef Functional_ Functional;
    static ProblemType const type = Functional::type;

    typedef typename Functional::Scalar Scalar;

    typedef typename Functional::AnsatzVars    AnsatzVars;
    typedef typename Functional::TestVars      TestVars;
    typedef typename Functional::OriginVars    OriginVars;

    typedef typename Functional::DomainCache   DomainCache;
    typedef typename Functional::BoundaryCache BoundaryCache;

    /**
     * \brief Type of variables around which the functional will be linearized.
     *
     * The type of "iterate", i.e. the point around which the functional
     * is linearized in the ansatz and test spaces. Note that this can
     * be, but need not be, the ansatz space itself. In any case, it must
     * refer to the same list of FE spaces.
     * 
     * \TODO: Why that? It must have the
     * same number of variables with the same number of components each,
     * but the underlying FE spaces can be different.
     */
    typedef typename Functional::OriginVars::VariableSet PointOfLinearization;

    /** 
     * \brief Construct a linearization from a given functional
     */
    LinearizationAt(Functional const & f_, PointOfLinearization const& u_) 
    : f(f_)
    , u(u_)
    { }
    
    /// create a domain cache
    DomainCache createDomainCache(int flags) const 
    { 
      return DomainCache(f,u,flags); 
    }
    
    /// create a boundary cache
    BoundaryCache createBoundaryCache(int flags) const 
    { 
      return BoundaryCache(f,u,flags); 
    }

    /**
     * \brief Rhs block information.
     * 
     * Inherits the provided data.
     */
    template <int row>
    struct D1: public Functional::template D1<row> {};

    /// matrix block information
    template <int row, int col>
    struct D2: public Functional::template D2<row,col> {};

    /// integration order
    template <class Cell>
    int integrationOrder(Cell const& cell, int shapeFunctionOrder, bool boundary) const
    {
      return f.integrationOrder(cell,shapeFunctionOrder,boundary);
    }

    template <int row, int col>
    auto makePositiveThreshold() const
    {
      return f.template makePositiveThreshold<row,col>();
    }

    /// Get point of linearization
    PointOfLinearization const& getOrigin() const 
    {
      return u;
    }

    /// Get functional
    Functional const& getFunctional() const 
    {
      return f;
    }
    
    /**
     * \brief Gives the highest derivative that arises in the weak formulation.
     */
    int derivatives() const 
    { 
      return f.derivatives(); 
    }

    protected:
    Functional const& f;
    PointOfLinearization const& u;
  };



  /**
   * \ingroup functional
   * \relates LinearizationAt
   * \brief Convenience routine: construct linearization without having to know the type of the functional
   */
  template<class Functional>
  LinearizationAt<Functional> linearization(Functional const& f, typename Functional::OriginVars::VariableSet const & u)
  {
    return LinearizationAt<Functional>(f,u);
  }

  /** 
   * \ingroup functional
   * \brief Proxy class for evaluating differences of linearizations of a nonlinear functional. 
   * 
   * Acceptance tests in Newton's method for minimization problems usually involves some kind 
   * of comparison \f$ f(u+\delta u) < f(u) \f$, where
   * \f[ f(u) = \sum_{i=1}^N f_i(u). \f]
   * Now taking the difference of both sums (large, because it's assembly over the grid) is 
   * numerically unstable for \f$ \delta u \f$ small. A more stable evaluation is
   * \f[ \sum_{i=1}^N (f_i(u)-f_i(u+\delta u)) > 0. \f]
   * This class supports the evaluation of the sum of differences.
   * 
   * According to the following idea: A VariationalFunctional defines a mapping \f$ F : U \to V \f$, which
   * can be evaluated and linearized in an appropriate ansatz space
   * \f$ U_a\subset U \f$, only if an argument \f$ u \in U_o \f$ in a -
   * possibly different - "origin" subspace \f$ U_o \subset U \f$ is
   * given.
   *
   * Hence, from a functional \f$ f \f$ and two variables \f$ u_1, u_2 \in U_o \f$ one has to
   * create a linearization of \f$ f(u_1)-f(u_2) \f$, which is represented by this class.
   */
  template<class Functional_>
  class LinearizationDifferenceAt : public Functional_Aux_Detail::InnerBoundaryPolicy<Functional_,LinearizationDifferenceAt<Functional_> >
  {
  public:
    typedef Functional_ Functional;
    static ProblemType const type = Functional::type;

    typedef typename Functional::Scalar Scalar;

    typedef typename Functional::AnsatzVars    AnsatzVars;
    typedef typename Functional::TestVars      TestVars;
    typedef typename Functional::OriginVars    OriginVars;


    /**
     * \brief Type of variables around which the functional will be linearized.
     *
     * The type of "iterate", i.e. the point around which the functional
     * is linearized in the ansatz and test spaces. Note that this can
     * be, but need not be, the ansatz space itself. It must have the
     * same number of variables with the same number of components each,
     * but the underlying FE spaces can be different.
     */
    typedef typename Functional::OriginVars::VariableSet PointOfLinearization;

    /** 
     * \brief Construct a difference of linearizations from a given functional
     */
    LinearizationDifferenceAt(Functional const & f_,PointOfLinearization const& u1_, PointOfLinearization const& u2_) 
    : f(f_), u1(u1_), u2(u2_) {};
    
    class DomainCache
    {
    public:
      DomainCache(Functional const& f, PointOfLinearization const& u1, PointOfLinearization const& u2, int flags)
      : dc1(f,u1,flags), dc2(f,u2,flags) {}
      
      template <class Cell>
      void moveTo(Cell const& entity) { dc1.moveTo(entity); dc2.moveTo(entity); }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators)  { dc1.evaluateAt(x,evaluators); dc2.evaluateAt(x,evaluators); }

      Scalar d0() const { return dc1.d0()-dc2.d0(); }
    
      template<int row, int dim> 
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
      d1 (VariationalArg<Scalar,dim> const& arg) const { return dc1.d1<row>(arg)-dc2.d1<row>(arg); }

      template<int row, int col, int dim> 
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const { return dc1.d2<row,col>(arg1,arg2)-dc2.d2<row,col>(arg1,arg2); }
      
      
    private:
      typename Functional::DomainCache dc1, dc2;
    };
    
    class BoundaryCache
    {
    public:
      BoundaryCache(Functional const& f, PointOfLinearization const& u1, PointOfLinearization const& u2, int flags)
      : bc1(f,u1,flags), bc2(f,u2,flags) {}
      
      template <class FaceIterator>
      void moveTo(FaceIterator const& entity) { bc1.moveTo(entity); bc2.moveTo(entity); }

      template <class Position, class Evaluators>
      void evaluateAt(Position const& x, Evaluators const& evaluators) { bc1.evaluateAt(x,evaluators); bc2.evaluateAt(x,evaluators); }

      Scalar d0() const { return bc1.d0()-bc2.d0(); }
    
      template<int row, int dim> 
      Dune::FieldVector<Scalar, TestVars::template Components<row>::m> 
      d1 (VariationalArg<Scalar,dim> const& arg) const { return bc1.d1<row>(arg)-bc2.d1<row>(arg); }

      template<int row, int col, int dim> 
      Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
      d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const { return bc1.d2<row,col>(arg1,arg2)-bc2.d2<row,col>(arg1,arg2); }
      
    private:
      typename Functional::BoundaryCache bc1, bc2;
    };

    /// create a domain cache
    DomainCache createDomainCache(int flags) const 
    {
      return DomainCache(f,u1,u2,flags);
    }
    
    /// create a boundary cache
    BoundaryCache createBoundaryCache(int flags) const 
    {
      return BoundaryCache(f,u1,u2,flags);
    }

    /**
     * \brief Rhs block information.
     */
    template <int row>
    struct D1: public Functional::template D1<row> {};

    /// matrix block information
    template <int row, int col>
    struct D2: public Functional::template D2<row,col> {};

    /// integration order
    template <class Cell>
    int integrationOrder(Cell const& cell, int shapeFunctionOrder, bool boundary) const
    {
      return f.template integrationOrder<Cell>(cell, shapeFunctionOrder, boundary);
    }

    /// Get functional
    Functional const& getFunctional() const {return f;}

  protected:
    Functional const& f;
    PointOfLinearization const& u1;
    PointOfLinearization const& u2;
  };



  /// Convenience routine: construct linearization without having to know the type of the functional
  template<class Functional>
  LinearizationDifferenceAt<Functional> linearizationDifference(Functional const& f, 
                                                                typename Functional::OriginVars::VariableSet const & u1, 
                                                                typename Functional::OriginVars::VariableSet const & u2)
  {
    return LinearizationDifferenceAt<Functional>(f,u1,u2);
  }

  //---------------------------------------------------------------------

  /**
   * \ingroup functional
   *
   * \brief Proxy class for the linearization of semi-linear time stepping schemes.
   *
   * The sole difference to LinearizationAt is the second linearization
   * point \arg uJ, which is used to compute the Jacobian \f$ J \f$
   * occuring in \f[ B \dot u - J u = f(u) - J u. \f]
   * 
   * \tparam Functional a weak formulation of a parabolic equation
   */
  template <class Functional>
  class SemiLinearizationAt: public LinearizationAt<Functional>
  {
  public:
    /**
     * \brief Constructor
     * \param time dependent functional
     * \param u_ point of linearization for the time differential operator (i.e. used for Cache::b2)
     * \param uJ_ point of linearization for spatial differential operator (i.e. used for Cache::d2)
     * \param du_ correction of last time step
     */
    SemiLinearizationAt(Functional const& f_,
                        typename Functional::AnsatzVars::VariableSet const& u_,
                        typename Functional::AnsatzVars::VariableSet const& uJ_,
                        typename Functional::AnsatzVars::VariableSet const& du_):
          LinearizationAt<Functional>(f_,u_), uJ(uJ_), du(du_) {}

    typename Functional::DomainCache createDomainCache(int flags) const 
    {
      return typename Functional::DomainCache(this->f,this->u,uJ,du,flags);
    }

    typename Functional::BoundaryCache createBoundaryCache(int flags) const 
    {
      return typename Functional::BoundaryCache(this->f,this->u,uJ,du,flags);
    }
    
  private:
    typename Functional::AnsatzVars::VariableSet const& uJ;
    typename Functional::AnsatzVars::VariableSet const& du;
  };

  template <class Functional>
  class SemiLinearizationAtInner: public LinearizationAt<Functional>
  {
  public:
    /**
     * \brief Constructor
     * \param time dependent functional
     * \param u_ point of linearization for the time differential operator (i.e. used for Cache::b2)
     * \param uJ_ point of linearization for spatial differential operator (i.e. used for Cache::d2)
     * \param du_ correction of last time step
     */
    SemiLinearizationAtInner(Functional const& f_,
                        typename Functional::AnsatzVars::VariableSet const& u_,
                        typename Functional::AnsatzVars::VariableSet const& uJ_,
                        typename Functional::AnsatzVars::VariableSet const& du_):
          LinearizationAt<Functional>(f_,u_), uJ(uJ_), du(du_) {}

    typename Functional::DomainCache createDomainCache(int flags) const 
    {
      return typename Functional::DomainCache(this->f,this->u,uJ,du,flags);
    }

    typename Functional::BoundaryCache createBoundaryCache(int flags) const 
    {
      return typename Functional::BoundaryCache(this->f,this->u,uJ,du,flags);
    }
    
    typename Functional::InnerBoundaryCache createInnerBoundaryCache(int flags) const 
    {
      return typename Functional::InnerBoundaryCache(this->f,this->u,uJ,du,flags);
    }
    
  private:
    typename Functional::AnsatzVars::VariableSet const& uJ;
    typename Functional::AnsatzVars::VariableSet const& du;
  };

  /**
   * \brief Convenience routine for constructing semi-linearizations.
   * \relates SemiLinearizationAt
   */
  template<class Functional>
  SemiLinearizationAt<Functional> semiLinearization(Functional const& f,
      typename Functional::AnsatzVars::VariableSet const& u,
      typename Functional::AnsatzVars::VariableSet const& uJ,
      typename Functional::AnsatzVars::VariableSet const& du)
  {
    return SemiLinearizationAt<Functional>(f,u,uJ,du);
  }

  //---------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------
  
  
} /* end of namespace Kaskade */

#endif
