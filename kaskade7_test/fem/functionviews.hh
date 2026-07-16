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

#ifndef FUNCTIONVIEWS_HH
#define FUNCTIONVIEWS_HH

#include <type_traits>

namespace Kaskade
{
  /**
   * @file
   * @brief  Some useful views on FunctionSpaceElement s
   * @author Anton Schiela
   *
   */

  /**
   * \ingroup functionviews
   * \brief Namespace of views on a FunctionSpaceElement
   * Views on Function SpaceElements can be used in a restricted sense
   * like FunctionSpaceElements. The idea is to transform a given
   * function to another function while evaluating it at a certain
   * point. Such a view has provide certain functionality:
   *
   * Public Types::
   *
   * - Space         (which is a FEFunctionSpace)
   * - Scalar        (return type, such as double)
   * - ValueType     (value, usually a Dune::FieldVector<Scalar,...>)
   *
   * Usually, views are templates over a Function, and they can inherit their types
   *
   * Member Functions
   *
   * - Space const& space() const  (returns the space, the function lives on)
   * - template<class SFS> int order(SFS const& sfs) const (order of integration)
   *
   * Usually these mamber functions pass over the values from the function that is viewed on
   * the order of integration may have to be changed, according to the transformation, done.
   *
   * - ValueType value(typename Function::Space::Evaluator const& evaluator) const
   *
   * return the transformed value. Usually value(evaluator) is called for the function
   * that is viewed on. Then its value is transformed.
   *
   * Examples for such views are given by AbsSquare, Gradient, Difference
   */
  namespace FunctionViews
  {
    /**
     * \ingroup functionviews
     * View on FunctionSpaceElements. If .value() is called, then the
     * result is the square norm of the value of the FunctionSpaceElement
     * this view was constructed with.
     *
     * For more information on the construction of such views cf. FunctionViews
     */
    template <typename Function>
    class AbsSquare
    {
    public:
      typedef typename Function::Space Space;
      typedef typename Function::Scalar Scalar;
      typedef typename Dune::FieldVector<Scalar,1> ValueType;

      AbsSquare(Function const& f_) : f(f_) {}

      Space const& space() const {return f.space();}

      template<class SFS>
      int order(SFS const& sfs) const {return 2*f.order(sfs); }

      template<class EPtr, class V>
      ValueType value(EPtr const& ci, V const& v) const
      {
        return f.value(ci,v).two_norm2();
      }

      ValueType value(typename Function::Space::Evaluator const& evaluator) const
      {
        return f.value(evaluator).two_norm2();
      }

    private:
      Function const& f;
    };

    // --------------------------------------------------------------------------------------------

    /** 
     * \ingroup functionviews
     * \brief Derivative of a given finite element function.
     *
     * For a finite element function \f$ \Omega \subset \mathbb{R}^n \to \mathbb{R}^m \f$, this
     * provides a function view with values in
     * - \f$ \mathbb{R}^{m} \f$ if \f$ m == 1 \f$
     * - \f$ \mathbb{R}^{n\times m} \f$ if \f$ m > 1 \f$ and `asVector==false`
     * - \f$ \mathbb{R}^{nm} \f$ if \f$ m > 1 \f$ and `asVector==true`,
     *   with columns of \f$ f_x \f$ concatenated.
     * 
     * For more information on the construction of such views cf. FunctionViews
     */
    template <typename Function, bool asVector=Function::components==1>
    class Gradient
    {
      static int const m = Function::components;

    public:
      typedef typename Function::Space Space;
      static int const dim = Space::dim;

      typedef typename Function::Scalar Scalar;


      using ValueType = std::conditional_t<asVector,       Dune::FieldVector<Scalar,dim*m>,
                                         /*otherwise*/     Dune::FieldMatrix<Scalar,m,dim>>;

      Gradient(Function const& f_) : f(f_)
      {}

      Space const& space() const 
      {
        return f.space();
      }

      template<class SFS>
      int order(SFS const& sfs) const
      {
        return std::max(f.order(sfs)-1,0);
      }

      /**
       * \brief evaluates the derivative
       * \return for scalar functions, a derivative vector, for vectorial components, the Jacobian
       */
      ValueType value(typename Function::Space::Evaluator const& evaluator) const
      {
        return convertToValueType(f.derivative(evaluator));
      }

      /**
       * \brief evaluates the derivative
       * \return for scalar functions, a derivative vector, for vectorial components, the Jacobian
       */
      template<class Cell, class Position>
      ValueType value(Cell const& cell, Position const& xi) const
      {
        return convertToValueType(f.derivative(cell,xi));
      }


    private:
      Function const& f;

      ValueType convertToValueType(Dune::FieldMatrix<Scalar,m,dim> const& fx) const
      {
        if constexpr (m==1)
          return fx[0];
        else if constexpr (asVector)
          return Dune::asVector(fx);
        else
          return fx;
      }
    };


    /**
     * \ingroup functionviews
     * \brief creates a function view of the gradient
     * \relates Gradient
     */
    template <class Function>
    Gradient<Function,false> gradient(Function const& f)
    {
      return Gradient<Function,false>(f);
    }

    /**
     * \ingroup functionviews
     * \brief creates a function view of the gradient
     * \relates Gradient
     */
    template <class Function>
    Gradient<Function,true> gradientAsVector(Function const& f)
    {
      return Gradient<Function,true>(f);
    }

    // --------------------------------------------------------------------------------------------

    /**
     * \ingroup functionviews
     * \brief A function view returning the square of the Frobenius norm of a function's derivative.
     * 
     * For more information on the construction of such views cf. FunctionViews
     */
    template <typename Function>
    class GradientAbsSquare
    {
    public:
      typedef typename Function::Space Space;
      typedef typename Function::Scalar Scalar;
      typedef Scalar  ValueType;

      GradientAbsSquare(Function const& f_) : f(f_)
      { }

      Space const& space() const 
      {
        return f.space();
      }

      template <class SFS>
      int order(SFS const& sfs) const 
      { 
        return 2*std::max(f.order(sfs)-1,0); 
      }

      ValueType value(typename Function::Space::Evaluator const& evaluator) const
      {
        return f.derivative(evaluator).frobenius_norm2();
      }

      template <class Cell, class Position>
      ValueType value(Cell const& cell, Position const& xi) const
      {
        f.derivative(cell,xi).frobenius_norm2();
      }

    private:
      Function const& f;
    };


  }

  /**
   * \ingroup functionviews
   * \brief Construct a View on two Functions without having to state the type of these explicitely
   */
  template<template<class Fct1,class Fct2> class View, class Fct1, class Fct2>
  View<Fct1, Fct2> makeView(Fct1 const& f1, Fct2 const& f2)
  {
    return View<Fct1,Fct2>(f1,f2);
  }

  /**
   * \ingroup functionviews
   * \brief Construct a View on a Function without having to state the type of it explicitely
   */
  template<template<typename Function> class View, typename Function>
  View<Function> makeView(Function const& f1)
  {
    return View<Function>(f1);
  }

  /**
   * \ingroup functionviews
   * \brief Namespace of classes that mimick a FunctionSpaceElement in a very weak sense
   *
   * Views on Function SpaceElements can be used in a very
   * restricted sense like FunctionSpaceElements. The weak views are useful
   * to scan certain properties of the grid, such as the level of the grid.
   * Such weak views can publicly inherit from ZeroOrderView
   *
   * Member types
   * 
   * - ValueType a Dune::FieldVector<Scalar,n> vector type
   * 
   * Member functions
   *
   * - ValueType value(Cell const& ci, V const&) const
   *
   */
  namespace WeakFunctionViews
  {
    /// Base class providing int order()
    struct ZeroOrderView
    {
      template<class SFS>
      int order(SFS const& sfs) const { return 0; }

      virtual ~ZeroOrderView() {};
    };

    struct ConstantGradient
    {
      typedef Dune::FieldVector<double,2> ValueType;

      template<class SFS>
      int order(SFS const& sfs) const { return 0; }

      template<class EPtr, class V>
      ValueType value(EPtr const& ci, V const&) const
      {
        ValueType v;
        v[0]=0;
        v[1]=1;
        return v;
      }
    };

    /// Return 1, if entity mightbecoarsened
    struct MightBeCoarsened : public ZeroOrderView
    {
      template<class EPtr, class V>
      double value(EPtr const& ci, V const&) const
      {
        return ci.mightBeCoarsened()/ci.geometry().volume();
      }
    };

    /// Return 1, if entity was refined
    struct WasRefined : public ZeroOrderView
    {
      template<class EPtr, class V>
      double value(EPtr const& ci, V const&) const
      {
        return ci.wasRefined()/ci.geometry().volume();
      }
    };

    /// Get level() of entity in Grid
    struct GridLevel : public ZeroOrderView
    {
      template<class EPtr, class V>
      double value(EPtr const& ci, V const&) const
      {
        return ci.level()/ci.geometry().volume();
      }
    };

    /// Return 1, if entity is constructed from red refinement
    struct IsRegular : public ZeroOrderView
    {
      template<class EPtr, class V>
      double value(EPtr const& ci, V const&) const
      {
        return ci.isRegular()/ci.geometry().volume();
      }
    };

  }

  template<class Grid, class T> class CellData;
  template<class Space> class LocalIntegral;

  ///Evaluate WeakFunctionViews and construct CellData
  template<class View, class Space>
  typename CellData<typename Space::Grid, Dune::FieldVector<double,1> >::CellDataVector evalCellProperties(Space const& space)
  {
    LocalIntegral<Space> localIntegral;
    typename CellData<typename Space::Grid, Dune::FieldVector<double,1> >::CellDataVector
    errorIndicator(localIntegral(View(),space));
    return errorIndicator;
  }

  /**
   * \ingroup functionviews
   * \brief (Scalar) product of two functions.
   *
   * This implements the Function concept. Given two Functions \f$
   * f_1,f2:\Omega\to K^m \f$, it represents the scalar product Function
   * \f$ g(x)=f_1(x)^H f_2(x) \f$.
   */
  template <class Function>
  class ScalarProductView
  {
  public:
    /**
     * The functions f1 and f2 have to exist during the lifetime of this view.
     */
    ScalarProductView(Function const& f1_, Function const& f2_): f1(f1_), f2(f2_) {}

    typedef typename Function::Scalar Scalar;
    static int const components = Function::components;
    typedef Dune::FieldVector<Scalar,components> ValueType;

    template <class Cell>
    int order(Cell const& cell) const
    {
      if (std::max(f1.order(cell),f2.order(cell)) < std::numeric_limits<int>::max())
        return f1.order(cell)+f2.order(cell);
      else
        return std::numeric_limits<int>::max();
    }

    template <class Cell>
    ValueType value(Cell const& cell,
        Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
    {
      return f1.value(cell,localCoordinate) * f2.value(cell,localCoordinate);
    }


  private:
    Function const& f1;
    Function const& f2;
  };
} /* end of namespace Kaskade */

#endif
