/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef NORMS_HH
#define NORMS_HH

/**
 * @file
 * @brief  Some norms for FunctionSpaceElement s or FunctionViews
 * @author Anton Schiela
 */

#include <boost/fusion/include/at_c.hpp>
#include <fem/integration.hh>
#include <fem/functionviews.hh>

namespace Kaskade
{
  /// L_2-norms
  struct L2Norm
  {
  /// Evaluation of square norm
    template<typename Function>
    typename Function::Scalar square(Function const& f) const
    {
      Integral<typename Function::Space> integral;
      FunctionViews::AbsSquare<Function> asf(f);
      return integral(asf);
    }

  /// Evaluation of norm
    template<typename Function>
    typename Function::Scalar operator()(Function const& f) const
    {
      return std::sqrt(square(f));
    }
  };


  /// H1-semi-norms
  struct H1SemiNorm
  {
  /// Evaluation of square norm
    template<typename Function>
    typename Function::Scalar square(Function const& f) const
    {
      Integral<typename Function::Space> integral;
      FunctionViews::GradientAbsSquare<Function> asf(f);
      return integral(asf);
  //    FunctionViews::Gradient<Function> gf(f);
  //    L2Norm l2Norm;
  //    return l2Norm.square(gf);
    }

  /// Evaluation of norm
    template<typename Function>
    typename Function::Scalar operator()(Function const& f) const
    {
      return std::sqrt(square(f));
    }
  };



  /// H1-norms
  struct H1Norm
  {
  /// Evaluation of square norm
    template<typename Function>
    typename Function::Scalar square(Function f)
    {
      H1SemiNorm h1Norm;
      L2Norm l2Norm;
      return h1Norm.square(f)+l2Norm.square(f);
    }

  /// Evaluation of norm
    template<typename Function>
    typename Function::Scalar operator()(Function f)
    {
      return std::sqrt(square(f));
    }
  };

  template<class Space> class LocalIntegral;
  template<class Grid, class T> class CellData;

  /// local (cellwise) H1-semi-norms
  template<class Function>
  [[deprecated("appears to be outdated and unused - candidate for removal")]]
  typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector localH1SemiNorm(Function const& f)
  {
    typedef typename Function::Space Space;
    typedef typename Space::Grid Grid;
    LocalIntegral<Space> localIntegral;
    typename CellData<Grid,typename Function::ValueType>::CellDataVector
      errorIndicator(localIntegral(
                       makeView<FunctionViews::AbsSquare>(makeView<FunctionViews::Gradient>(f))));
    return errorIndicator;
  }

  /// local (cellwise) L2-norms
  template<class Function>
  [[deprecated("appears to be outdated and unused - candidate for removal")]]
  typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector localL2Norm(Function const& f)
  {
    typedef typename Function::Space Space;
    typedef typename Space::Grid Grid;
    LocalIntegral<Space> localIntegral;
    typename CellData<Grid,typename Function::ValueType>::CellDataVector
      errorIndicator(localIntegral(
                       makeView<FunctionViews::AbsSquare>(f)));
    return errorIndicator;
  }


  /**
   * @brief boundaryL2Norm computes the L2-norm of an FE function on the whole boundary of the underlying grid.
   *
   * @param function is the FE function to be integrated.
   * @tparam FEFunction is the type of the integrand.
   * @return value of integral
   *
   */
  template <class FEFunction>
  auto boundaryL2Norm(FEFunction const& function)
  {
    FunctionViews::AbsSquare<FEFunction> functionSquared(function);
    return std::sqrt(integrateOverBoundary(functionSquared));
  }
} // end of namespace Kaskade
#endif
