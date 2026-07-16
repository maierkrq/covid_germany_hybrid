/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DIRECT_HH
#define DIRECT_HH

#include <cmath>
#include <memory>
#include <type_traits>

#include <boost/timer/timer.hpp>

#include "dune/istl/operators.hh"
#include <dune/istl/solvers.hh>

#include "fem/istlinterface.hh"
#include "linalg/factorization.hh"

namespace Kaskade
{
  /// \internal
  namespace DirectSolver_Detail
  {
    class NoScalar{};
    class NoFieldType{};

    template <class T>
    typename T::Scalar hasScalar(typename T::Scalar*);
    template <class T>
    NoScalar hasScalar(...);

    template <class T>
    typename T::field_type hasFieldType(typename T::field_type*);
    template <class T>
    NoFieldType hasFieldType(...);

    template <class T>
    struct HasScalar
    {
      typedef decltype(hasScalar<T>(nullptr)) type;
      static constexpr bool value = !std::is_same<type,NoScalar>::value;
    };

    template <class T>
    struct HasFieldType
    {
      typedef decltype(hasFieldType<T>(nullptr)) type;
      static constexpr bool value = !std::is_same<type,NoFieldType>::value;
    };



    template <class T>
    struct ScalarType
    {
      typedef typename std::conditional< HasScalar<T>::value, typename HasScalar<T>::type,
                                           typename std::conditional< HasFieldType<T>::value, typename HasFieldType<T>::type, void>::type >::type type;

    };

    /**
     * \brief Checks a vector for nan/inf and reports their number to stderr.
     */
    template <class Scalar>
    size_t checkNanInf(Scalar const* x, size_t n, std::string const& what)
    {
      size_t nan = 0, inf = 0;
      for (size_t i=0; i<n; ++i)
      {
        if (std::isnan(x[i])) ++nan;
        if (std::isinf(x[i])) ++inf;
      }
      if (nan>0) std::cerr << nan << " " << what << " entries are nan.\n";
      if (inf>0) std::cerr << inf << " " << what << " entries are inf.\n";

      return nan+inf;
    }
  }
  /// \endinternal

  /**
   * \ingroup direct
   * \brief Dune::InverseOperator and Dune::Preconditioner interface for direct solvers.
   *
   * This keeps a factorization during the lifetime of the object, however, due to
   * shared data, efficient copying is possible.
   */
  template <class Domain_, class Range_>
  class DirectSolver: public Dune::InverseOperator<Domain_,Range_>,
                      public Dune::Preconditioner<Domain_,Range_>
  {
  public:
    typedef typename DirectSolver_Detail::ScalarType<Domain_>::type Scalar;
    typedef Domain_ Domain;
    typedef Range_ Range;

    /**
     * \brief Default constructor. 
     * 
     * A default constructed DirectSolver implements a (pretty useless) zero operator.
     */
    DirectSolver() {}
    
    /**
     * \brief Copy constructor.
     * 
     * The copy shares the factorization with the source.
     */
    DirectSolver(DirectSolver<Domain,Range> const& ds) = default;

    /**
     * \brief Constructs a direct solver from an assembled linear operator. 
     * 
     * A copy of the linear operator is held internally, thus the linear
     * operator object is better not too large.
     */
    template <class AssembledGOP>
    explicit DirectSolver(AssembledGOP const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
    : DirectSolver(A.template get<MatrixAsTriplet<typename AssembledGOP::Scalar>>(),directType,properties)
    { }

    /**
     * \brief Constructs a direct solver from a triplet matrix. 
     */
    template <class FieldType>
    explicit DirectSolver(MatrixAsTriplet<FieldType> const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
    : solver(getFactorization(directType,A,typename Factorization<Scalar>::Options(properties)))
    , n(A.N())
    {
      assert(solver.get()); // make sure the factorization has been successful - otherwise it should have raised an exception
    }

    /**
     * \brief Constructs a direct solver from a BCRS matrix. 
     * A copy of the linear operator is held internally, thus the linear
     * operator object is better not too large.
     */
    template <class FieldType>
    explicit DirectSolver(Dune::BCRSMatrix<Dune::FieldMatrix<FieldType,1,1>> const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
    : solver(getFactorization<FieldType,int>(directType,A,typename Factorization<Scalar>::Options(properties)))
    , n(A.N())
    {
      assert(solver.get()); // make sure the factorization has been successful - otherwise it should have raised an exception
    }

    /**
     * \brief Constructs a direct solver from a NumaBCRS matrix. 
     * A copy of the linear operator is held internally, thus the linear
     * operator object is better not too large.
     *
     * Warning: Don't use this if block size is greater than 1, DirectType is MUMPS and the blocks on the diagonal are not diagonal,
     * since MUMPS requires that only the lower triangle of the matrix is given.
     * Todo: Throw exception in this case?
     */
    template <class FieldType, int n, int m, class Index>
    explicit DirectSolver(NumaBCRSMatrix<Dune::FieldMatrix<FieldType,n,m>,Index> const& A, 
                          DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
    : solver(getFactorization<FieldType,n,m,Index,int>(directType,A,typename Factorization<Scalar>::Options(properties)))
    , n(A.N()*n)
    {
      assert(solver.get()); // make sure the factorization has been successful - otherwise it should have raised an exception
    }

    /**
     * \name Solution
     * These methods apply the inverse operators, i.e. they solve the linear system. We
     * provide several different interfaces such as Dune::InverseOperator and Dune::Preconditioner
     * and some more (in particular const methods) for convenience. All the methods do
     * essentially the same - they solve the system.
     * @{
     */

    /**
     * \brief Solves the system for the given right hand side `b`.
     * 
     * The provided right hand side `b` is not modified (except possibly via aliasing
     * `x`). This non-const method is provided just for interface compatibility with
     * Dune::InverseOperator.
     */
    virtual void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res)
    {
      apply(x,b,1e-10,res);
    }

    /**
     * \brief Solves the system for the given right hand side `b`.
     *
     * The provided right hand side `b` is not modified (except possibly via aliasing
     * `x`).
     */
    void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res) const
    {
      apply(x,b,1e-10,res);
    }

    /**
     * \brief Solves the system for the given right hand side `b`.
     *
     * The provided right hand side `b` is not modified (except possibly via aliasing
     * \arg x). This non-const method is provided just for interface compatibility with
     * Dune::InverseOperator.
     */
    virtual void apply(Domain& x, Range& b, double reduction, Dune::InverseOperatorResult& res)
    {
      const_cast<DirectSolver<Domain,Range> const*>(this)->apply(x,b,res);
    }

    /**
     * \brief Solves the system for the given right hand side `b`.
     */
    void apply(Domain& x, Range const& b, double reduction, Dune::InverseOperatorResult& res) const
    {
      if constexpr (false && contiguousStorage<Domain> && contiguousStorage<Range>)
      {
        auto xp = DynamicMatrixDetail::getAddress(x);
        auto bp = DynamicMatrixDetail::getAddress(b);
        assert(DirectSolver_Detail::checkNanInf(bp,n,"rhs")==0);
        solver->solve(bp,xp);
        assert(DirectSolver_Detail::checkNanInf(xp,n,"solution")==0);
      }
      else
      {
        std::vector<Scalar> rhs(b.dim()+16);    // allocate somewhat more memory - at least UMFPACK
        vectorToSequence(b,begin(rhs));         // tends to read beyond the array boundary?
        apply(rhs,reduction,res);
        vectorFromSequence(x,begin(rhs));
      }
    }

    /**
     * \brief The input vector v contains the rhs. It will be overwritten by the solution.
     */
    void apply(std::vector<Scalar>& v, double /* reduction */, Dune::InverseOperatorResult& res) const
    {
      boost::timer::cpu_timer timer;

      apply(v);

      // Write solver statistics. Currently dummy.
      res.clear();
      res.iterations = 1;
      res.reduction = 1e-10; // dummy!
      res.converged = true;
      res.conv_rate = 1e-10;
      res.elapsed = (double)(timer.elapsed().user)/1e9;
    }

    /**
     * \brief Solves the system \f$ Ax = b \f$
     *
     * \param[inout] v on entry, contains right hand side \f$ b \f$
     *                 on exit, contains the solution \f$ x \f$
     *                 v has to be large enough, at least the matrix size
     */
    virtual void  apply(std::vector<Scalar>& v) const
    {
      // Check for solver availability. Otherwise return zero solution.
      if (solver)
      {
        assert(v.size() >= n);
        assert(DirectSolver_Detail::checkNanInf(&v[0],n,"rhs")==0);

        solver->solve(v);

        assert(DirectSolver_Detail::checkNanInf(&v[0],n,"solution")==0);
      }
      else
        std::fill(v.begin(),v.end(),0.0);
    }

    /**
     * \brief Solves the system for the given right hand side `b`.
     *
     * The provided right hand side `b` is not modified (except possibly via aliasing
     * `x`). This non-const method is provided just for interface compatibility with
     * Dune::Preconditioner (virtual function).
     */
    virtual void apply(Domain& x, Range const& b)
    {
      Dune::InverseOperatorResult dummy_res;
      apply(x,b,0,dummy_res);
    }

    /**
     * \brief Solves the system for the given right hand side `b`.
     */
    virtual void apply (Domain& x, Range const& b) const
    {
      Dune::InverseOperatorResult dummy_res;
      apply(x,b,0,dummy_res);
    }

    /// @}

    /**
     * \brief unused
     *
     * This method does exactly nothing, and is just present for interface conformity.
     */
    virtual void pre (Domain &x, Range &b) {}

    /**
     * \brief unused
     *
     * This method does exactly nothing, and is just present for interface conformity.
     */
    virtual void post (Domain &x) {}

    /**
     * \brief returns the category of the operator
     * 
     * From the Dune doxygen documentation it is unclear what this is supposed to mean. We return a dummy here.
     */
    virtual Dune::SolverCategory::Category category() const override
    {
      return Dune::SolverCategory::sequential;
    }

  private:
    std::shared_ptr<Factorization<Scalar>> solver;
    size_t n = 0;   // size of the square matrix (with scalar entries)
  };

  //---------------------------------------------------------------------
  
  // partial specialization for tiny fixed-size matrices
  template <class S, int n>
  class DirectSolver<Dune::FieldVector<S,n>,Dune::FieldVector<S,n>>
  : public Dune::InverseOperator<Dune::FieldVector<S,n>,Dune::FieldVector<S,n>>
  , public Dune::Preconditioner<Dune::FieldVector<S,n>,Dune::FieldVector<S,n>>
  {
  public:
    typedef S Scalar;
    typedef Dune::FieldVector<S,n> Domain;
    typedef Dune::FieldVector<S,n> Range;

    /**
     * \brief Default constructor. 
     * 
     * A default constructed DirectSolver implements a (pretty useless) zero operator.
     */
    DirectSolver() {}

    /**
     * \brief Constructs a direct solver from a FieldMatrix matrix. 
     */
    template <class FieldType>
    explicit DirectSolver(Dune::FieldMatrix<FieldType,n,n> const& A):
        inverse(A)
    {
      inverse.invert();
    }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * The provided right hand side \arg b is not modified (except possibly via aliasing
     * \arg x).
     */
    virtual void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res) { apply(x,b,1e-10,res); }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * As an extension to the Dune::InverseOperator interface, we provide const versions of
     * the apply methods.
     */
    void apply(Domain& x, Range& b, Dune::InverseOperatorResult& res) const { apply(x,b,1e-10,res); }

    /**
     * Solves the system for the given right hand side \arg b, which is
     * guaranteed not to be overwritten (except possibly via aliasing
     * \arg x).
     */
    virtual void apply(Domain& x, Range& b, double reduction, Dune::InverseOperatorResult& res)
    {
      const_cast<DirectSolver<Domain,Range> const*>(this)->apply(x,b,res);
    }

    /**
     * \brief Solves the system for the given right hand side \arg b.
     * 
     * As an addition to the Dune::InverseOperator interface, we provide const versions of
     * the apply methods, since the factorization is not modified during the solution phase.
     */
    void apply(Domain& x, Range const& b, double reduction, Dune::InverseOperatorResult& res) const
    {
      x = inverse*b;
      
      // Write solver statistics. Currently dummy.
      res.clear();
      res.iterations = 1;
      res.reduction = 1e-10; // dummy!
      res.converged = true;
      res.conv_rate = 1e-10;
      res.elapsed = 0;
    }


    virtual void        pre (Domain &x, Range &b) {}

    virtual void        apply (Domain &v, const Range &d)
    {
      v = inverse*d;
    }

    virtual void        post (Domain &x) {}

  private:
    Dune::FieldMatrix<Scalar,n,n> inverse;
  };

  //---------------------------------------------------------------------

  /**
   * \ingroup linalgsolution
   * \brief Dune::LinearOperator interface for inverse operators.
   */
  template <class InverseOperator>
  class InverseLinearOperator: public Dune::LinearOperator<typename InverseOperator::Range, typename InverseOperator::Domain>
  {
  public:
    typedef typename InverseOperator::Domain Domain;
    typedef typename InverseOperator::Range Range;
    typedef typename InverseOperator::Scalar Scalar;

    InverseLinearOperator() = default;

    InverseLinearOperator(InverseOperator const& op_) : op(op_)
    {}

    virtual void apply(Domain const& x, Range& y) const
    {
      Dune::InverseOperatorResult result;
      Domain rhs(x);
      op.apply(y,rhs,result);
    }

    virtual void applyscaleadd(Scalar alpha, Domain const& x, Range& y) const
    {
      Range ynew(y);
      apply(x,ynew);
      y.axpy(alpha,ynew);
    }

    /**
     * \brief returns the category of the operator
     * 
     * From the Dune doxygen documentation it is unclear what this is supposed to mean. We return a dummy here.
     */
    virtual Dune::SolverCategory::Category category() const override
    {
      return Dune::SolverCategory::sequential;
    }

  private:
    InverseOperator op;
  };

  //---------------------------------------------------------------------

  /**
   * \ingroup direct
   * \brief convenience function for constructing a DirectInverseOperator
   * \relates InverseLinearOperator
   */
  template <class GOP, int firstRow, int lastRow, int firstCol, int lastCol>
  InverseLinearOperator<DirectSolver<typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Domain,
                                     typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Range> >
  directInverseOperator(AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol> const& A,
                        DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
  {
    typedef typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Domain Domain;
    typedef typename AssembledGalerkinOperator<GOP,firstRow,lastRow,firstCol,lastCol>::Range Range;
    return InverseLinearOperator<DirectSolver<Domain,Range> >(DirectSolver<Domain,Range>(A,directType,properties));
  }

  template <class Matrix, class Domain, class Range>
  InverseLinearOperator<DirectSolver<Domain,Range> >
  directInverseOperator(MatrixRepresentedOperator<Matrix,Domain,Range> const& A,
                        DirectType directType=DirectType::UMFPACK, MatrixProperties properties=MatrixProperties::GENERAL)
  {
    return InverseLinearOperator<DirectSolver<Domain,Range> >(DirectSolver<Domain,Range>(A,directType,properties));
  }

} // namespace Kaskade

//---------------------------------------------------------------------

#endif
