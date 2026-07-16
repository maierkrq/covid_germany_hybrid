/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2017-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MATRIXPRODUCT_HH
#define MATRIXPRODUCT_HH

#include "threadedMatrix.hh"

namespace Kaskade 
{
  /**
   * \ingroup linalgbasic
   * \relates NumaBCRSMatrix
   * \brief Computes the matrix-matrix product \f$ (A,B) \mapsto AB \f$
   * 
   * \tparam EntryA the type of entries of matrix A. Usually Dune::FieldMatrix<Real,na,ma>.
   * \tparam EntryB the type of entries of matrix B. Usually Dune::FieldMatrix<Real,nb,mb>.
   * 
   * As matrix entries \f$ A_{ij} B_{jk}\f$ have to be multiplied, the entry size must be compatible,
   * i.e. ma==nb. As an extension, if entries of A or B are scalar (i.e. na==ma==1), we treat them
   * implicitly as multiples of the unit matrix of compatible size. 
   */
  template <class EntryA, class EntryB, class IndexA, class IndexB>
  auto operator*(NumaBCRSMatrix<EntryA,IndexA> const& A, NumaBCRSMatrix<EntryB,IndexB> const& B)
  {
    // Compute the resulting entry size, treating the case of scalar entries.
    constexpr bool scalarA = EntryA::rows*EntryA::cols==1;
    constexpr bool scalarB = EntryB::rows*EntryB::cols==1;
    using Entry = std::conditional_t<scalarA,
                                     EntryB,
                                     std::conditional_t<scalarB,
                                                        EntryA,
                                                        Dune::FieldMatrix<typename EntryA::field_type,EntryA::rows,EntryB::cols>>>;
    static_assert(scalarA || scalarB || (int)EntryA::cols==(int)EntryB::rows, // cols and rows are enums, comparison gives warning
                  "Entry size has to match");
    
    assert(A.M()==B.N());
    
    // TODO: parallelize?
    NumaCRSPatternCreator<IndexA> creator(A.N(),B.M(),false,2);
    for (IndexA i=0; i<A.N(); ++i)
    {
      auto rowA = A[i];
      for (auto cai=rowA.begin(); cai!=rowA.end(); ++cai)
      {
        auto rowB = B[cai.index()];
        for (auto cbi=rowB.begin(); cbi!=rowB.end(); ++cbi)
          creator.addElement(i,cbi.index());
      }
    }
    
    // TODO: parallelize
    NumaBCRSMatrix<Entry,IndexA> C(creator);
    for (IndexA i=0; i<A.N(); ++i)
    {
      auto rowA = A[i];
      auto rowC = C[i];
      for (auto cai=rowA.begin(); cai!=rowA.end(); ++cai)
      {
        auto rowB = B[cai.index()];
        for (auto cbi=rowB.begin(); cbi!=rowB.end(); ++cbi)
          rowC[cbi.index()] += normalForm(*cai) * normalForm(*cbi);       // use normalForm() here to interpret 1x1 blocks as scalars
      }
    }
    
    return C;
  }
}

#endif
