/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2019-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "fem/diffops/contactConstraints.hh"

namespace Kaskade
{
  template <class Scalar, int dim>
  Mortar<Scalar,dim>::~Mortar()
  {}

  // explicit instantiation
  template class Mortar<double,2>;
  template class Mortar<double,3>;

  // ----------------------------------------------------------------------------------------------

  template <class Scalar, int dim>
  void DiracMortar<Scalar,dim>::updateConstraints(int , Dune::FieldVector<Scalar,dim-1> const& ,
                                                  typename Mortar<Scalar,dim>::Row& vic, Scalar g,
                                                  std::vector<typename Mortar<Scalar,dim>::Row>& rows,
                                                  std::vector<Scalar>& bounds) const
  {
    // That's easy: Every sample is an independent constraint. Just collect them.
    rows.push_back(vic);
    bounds.push_back(g);
  }

  // ----------------------------------------------------------------------------------------------

  namespace ContactConstraintsDetail
  {
    namespace
    {
      DiracMortar<double,2> defaultMortar2;
      DiracMortar<double,3> defaultMortar3;
    }

    template <class Scalar, int dim>
    constexpr DiracMortar<Scalar,dim> const& defaultMortar()
    {
      if constexpr (dim==2)
        return defaultMortar2;
      else
        return defaultMortar3;
    }

    // explicit instantiation
    template DiracMortar<double,2> const& defaultMortar<double,2>();
    template DiracMortar<double,3> const& defaultMortar<double,3>();
  }


  // ----------------------------------------------------------------------------------------------

  template <class Scalar, int dim>
  BezierMortar<Scalar,dim>::BezierMortar(int m_)
  : m(m_)
  {
    assert(m>=0);
  }

  template <class Scalar, int dim>
  int BezierMortar<Scalar,dim>::size(int ) const
  {
    return numberOfMultiindices(dim-1,m);
  }

  template <class Scalar, int dim>
  void BezierMortar<Scalar,dim>::updateConstraints(int n, Dune::FieldVector<Scalar,dim-1> const& xi,
                                                   typename Mortar<Scalar,dim>::Row& Gi,
                                                   Scalar g,
                                                   std::vector<typename Mortar<Scalar,dim>::Row>& rows,
                                                   std::vector<Scalar>& bounds) const
  {
    int const count = size(m);

    // We update all constraints simultaneously if a new sample arrives.
    if (rows.empty())
    {
      rows.insert(begin(rows),count,Row());
      bounds.clear();
      bounds.insert(begin(bounds),count,0.0);
    }
    assert(rows.size()==count);
    assert(bounds.size()==count);
    assert(n>0);

    for (int i=0; i<count; ++i)        // add the integral contribution to all mortar constraints,
    {
      Scalar bw = bezier(m,i,xi) / n; // get the Bezier weight, i.e. the test function value times the integration factor
      for (auto const& p: Gi)         // addition of row by appending duplicates (within the row)
        rows[i].push_back(std::make_pair(p.first,bw*p.second));
      bounds[i] += bw*g;
    }
  }

  // explicit instantiation
  template class BezierMortar<double,2>;
  template class BezierMortar<double,3>;

}
