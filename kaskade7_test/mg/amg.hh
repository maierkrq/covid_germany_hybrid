/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2018-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef AMG_HH
#define AMG_HH

#include <tuple>

#include "fem/gridmanager.hh"
#include "fem/gridCover.hh"
#include "linalg/threadedMatrix.hh"
#include "mg/prolongation.hh"
#include "fem/cellLocator.hh"
#include "utilities/scalar.hh"

namespace Kaskade
{
  
  using Prolongation = NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  
  /**
   * \ingroup fetransfer
   * \brief Computes an interpolation matrix for transfer between P1 finite element spaces 
   *        on two different grids.
   * 
   * \tparam FromSpace a FEFunctionSpace of continuous piecewise linear elements
   * \tparam ToSpace a FEFunctionSpace of continuous piecewise linear elements
   * 
   * \param from the coarse space
   * \param to the fine space
   * 
   * The computed interpolation matrix transfers coefficient vectors of the coarse space to 
   * coefficient vectors of the fine space such that the function is (approximately) the same. 
   * Thus, it is a representation of (approximate) identity w.r.t. the Lagrangian bases of the 
   * involved spaces.
   * 
   * An example use is in auxiliary space methods, where 
   * 
   * \todo generalize to arbitrary spaces, not only P1?
   */
  template <class FromSpace, class ToSpace>
  Prolongation twoGridProlongation(FromSpace const& from, ToSpace const& to)
  {
    using Real = typename ScalarTraits<typename ToSpace::Scalar>::Real;
    
    // Prolongation matrix has as many rows as the target (fine) space dimension, 
    // and columns as the source (coarse) dimension
    NumaCRSPatternCreator<> creator(to.degreesOfFreedom(),from.degreesOfFreedom());
    
    // interpolate all vertices in the target grid
    ShapeFunctionCache<typename FromSpace::Grid,Real> sfCache;   // create a shape function cache 
    auto eval = from.evaluator(&sfCache,0);                      // to be used by the evaluator
    std::vector<std::tuple<size_t,size_t,Real>> entries;
    
    for (auto const& tv: Dune::vertices(to.gridView()))          // step through all the vertces (i.e. Lagrange 
    {                                                            // basis functions) of the target grid 
      auto x = tv.geometry().center();                           
      auto const& sc = findCell(from.gridView(),x);
      eval.moveTo(sc);
      eval.evaluateAt(sc.geometry().local(x));
      auto sfCoeff = from.linearCombination(eval);               // check the contribution of source dofs to 
                                                                 // that vertex.
      auto row = to.indexSet().index(tv);
      for (auto const& pair: sfCoeff)
      {
        auto col = pair.first;
        creator.addElement(row,col);
        entries.push_back(std::make_tuple(row,col,pair.second[0]));
      }
    }
    
    // Create and fill the prolongator
    NumaBCRSMatrix<Dune::FieldMatrix<Real,1,1>> P(creator);
    for (auto const& entry: entries)
      P[std::get<0>(entry)][std::get<1>(entry)] = std::get<2>(entry);
    
    return P;
  }
  
  
  /**
   * \ingroup fetransfer
   * \brief Computes an interpolation matrix for transfer between P1 finite element spaces on two different (possibly non-nested) grids.
   * 
   * \tparam FromSpace a FEFunctionSpace of continuous piecewise linear elements
   * \tparam ToSpace a FEFunctionSpace of continuous piecewise linear elements
   * 
   * \param from      the coarse space
   * \param to        the fine space
   * \param tolTrunc  tolerance for truncation of P-matrix entries (entries below tol*rowtotal will be deleted)
   * 
   * The computed interpolation matrix transfers coefficient vectors of the coarse space to coefficient vectors of the 
   * fine space such that the function is (approximately) the same. Thus, it is a representation of (approximate)
   * identity w.r.t. the Lagrangian bases of the involved spaces.
   * Since the spaces might be non-nested, there are two undesired cases that can come up:
   * 1. A fine grid vertex is not contained in the coarse cell:
   *    To deal with this we express the fine grid DOF wrt to the DOFs "closest" to the corresponsing vertex.
   * 
   * 2. A coarse grid cell/DOF is not contained in the fine cell: 
   *    This introduces a zero column in the prolongation matrix, which may lead to undefinite Galerkin projections. 
   *    Use the nonNestedProlognationStack() function to remove zero columns from
   *    a resulting prolognation matrix or a prolongation stack.
   * 
   * To retain good convergence of MG methods based on this prolongation matix, a truncation of "small" entries is 
   * required. In each row, we delete all entries, which are smaller than tolTrunc*rowtotal, with rowtotal being the 
   * sum of the absolute values of the row entries. Afterwards we scale each row to retain its original row total.
   * The value of 0.2 seems to be a good guess, which balances the two competing goals of operator complexity and a 
   * good convergence rate. 
   * 
   * 
   * Cf. Dickopf, Krause - "Numerical Study of the Almost Nested Case in a Multlevel Method Based on Non-nested Meshes" (2013)
   * 
   */
  template <class FromSpace, class ToSpace>
  Prolongation nonNestedProlongation(FromSpace const& from, ToSpace const& to, double tolTrunc=0.1)
  {
    using LeafView = typename FromSpace::Grid::LeafGridView;
    using Real = typename ScalarTraits<typename ToSpace::Scalar>::Real;
    assert(("Coarse space should have less DOFs than fine space.", from.degreesOfFreedom()<to.degreesOfFreedom()));

    
    // Prolongation matrix has as many rows as the fine space dimension, 
    // and columns as the coarse dimension    
    NumaCRSPatternCreator<> creator(to.degreesOfFreedom(),from.degreesOfFreedom());
    ShapeFunctionCache<typename FromSpace::Grid,Real> sfCache;              // create a shape function cache 
    auto eval = from.evaluator(&sfCache,0);                                 // to be used by the evaluator
    std::vector<std::tuple<size_t,size_t,Real>> entries;                    // tuple = (row,column,value)
    
    // index set of fine space elements
    const auto& idxCoarse = from.gridView().indexSet();
    
    // cell locator for finding closest cells to point
    CellLocator<LeafView> cellLocator(from.gridView());
    
    // each row might be needed to scale to keep the
    // row total constant after truncation
    std::vector<double> rowscale(to.degreesOfFreedom(),1.0);
    
    // step through all vertices of fine grid
    for( const auto& fVert : Dune::vertices(to.gridView()) )
    {
      // global coordinate of vertex
      auto xglob = fVert.geometry().center();
      
      // find coarse grid cell containing xglob
      // (or closest cell, if not contained in any)
      auto [cell,dist] = cellLocator.closestCell(xglob); 
      
      // check contribution of coarse DOFs to (fine grid) vertex
      eval.moveTo(cell);
      eval.evaluateAt(cell.geometry().local(xglob));
      auto coeff = from.linearCombination(eval);
      
      // get row index
      auto row = to.indexSet().index(fVert);
      
      // get row total, maximal row entry and
      // rowdiff = sum of all truncated values
      double rowtotal(0.), rowMax(0.), rowdiff(0.);
      for(const auto& p : coeff)
      {
        double val = p.second[0];
        rowtotal += fabs(val);
        if(fabs(val) > rowMax)
          rowMax = fabs(val);
      }

      // add entries to pattern creator, truncating small values
      // (those smaller than tolTunc*rowtotal)
      // scaling happens when we insert values into Numa matrix
      for(const auto& p:coeff)
      {
        auto col = p.first;
        double val = p.second[0];
        
        if(fabs(val) >= 1e-10 && fabs(val) >= tolTrunc*rowtotal)        // alternative: if(fabs(val) >= tolTrunc*rowMax)
        {
          creator.addElement(row,col);
          entries.push_back(std::make_tuple(row,col,val));
        } else {
          rowdiff += fabs(val);
        }
      }
      
      // we keep the row total unchanged
      // by scaling the truncated rows
      assert(("Truncation deletes all row entries! Use smaller tolerance or other criterion for truncation.", rowdiff < rowtotal ));
      if(tolTrunc>0 && rowdiff != 0)
        rowscale.at(row) = rowtotal/(rowtotal-rowdiff);
    }
    
    // Create and fill the prolongation matrix
    NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>> P(creator);
    for (auto const& entry: entries)
    {
      P[std::get<0>(entry)][std::get<1>(entry)] = rowscale.at(std::get<0>(entry))*std::get<2>(entry);
    }
    
    return P;
  }
  
  /**
   * \ingroup fetransfer
   * \brief Computes a full column rank prolongation stack from a colum rank deficient prolongation stack.
   * 
   * \param ps vector containing the prolognation stack.
   * 
   * The prolognation stack obtained from non-nested meshes may contain zero columns (in particular if a coarse grid cell is not contained
   * in the fine mesh). 
   * To avoid this, we "erase" the corresponding DOFs on an algebraic level, i.e. by deleting the corresponding column in the prolognation matrix.
   * In the prolognation matrix on the next coarser level in the stack we need to delete the corresponding rows, too.
   * Doing this on each of the levels we obtain a full column rank prolognation stack.
   * 
   *  
   */
  std::vector<Prolongation> 
  nonNestedProlognationStack(std::vector<Prolongation> ps)
  {
    // vector to return, contains the prolongation stack
    // with each P-matrix having full column rank
    std::vector<Prolongation> retPs;
    
    // initalize number of levels and prolongation matrix
    // and vector possibly containing inidices of nonzero columns
    int levels = ps.size();
    Prolongation Pl;
    // initialite nonzero columns vector
    // contaning all indices from 0 to P.M()
    std::vector<size_t> nzCols(ps.at(levels-1).N());
    
    // loop though the stack (not yet freed of zero columns)
    for(int l=levels-1; l>=0; --l)
    {
      // P matrix (with possible zero columns)
      // and vector of columns indices
      Pl = ps.at(l);
      std::vector<size_t> allCols(Pl.M());
      std::iota(allCols.begin(),allCols.end(),0);    
      
      // delete rows corresponding to zero columns of
      // fine space (zeroCols of previous iteration)
      if(nzCols.size() < Pl.N() && l!=levels-1)
        Pl = submatrix<Prolongation>(Pl,nzCols,allCols);
      
      auto& timer = Timings::instance();
      // check for nonzero columns in current matrix
      timer.start("seeking nonzero columns");
      nzCols = nonZeroColumns(Pl);
      timer.stop("seeking nonzero columns");
      std::cout << "Found " << Pl.M() - nzCols.size() << " zero columns in P-matrix on level " << l << ". Removing... ";
      
      // if we have found zero columns "delete" them
      if(nzCols.size() < Pl.M())
      {
        std::vector<size_t> allRows(Pl.N());
        std::iota(allRows.begin(),allRows.end(),0); 
        Pl = submatrix<Prolongation>(Pl,allRows,nzCols);
      }
      std::cout << "done! \n";
      
      // insert in front position of return vector
      retPs.insert(retPs.begin(), Pl);
    }
    return retPs;
  }
  
  
  /**
   * \ingroup fetransfer
   * \brief Computes an interpolation matrix for transfer between P1 finite element spaces on two different (possibly non-nested) grids.
   * 
   * \tparam FromSpace a FEFunctionSpace of continuous piecewise linear elements
   * \tparam ToSpace a FEFunctionSpace of continuous piecewise linear elements
   * 
   * \param from the coarse space
   * \param to the fine space
   * 
   * The computed interpolation matrix transfers coefficient vectors of the coarse space to coefficient vectors of the 
   * fine space such that the function is (approximately) the same. Thus, it is a representation of (approximate)
   * identity w.r.t. the Lagrangian bases of the involved spaces.
   * Since the spaces might be non-nested, there are two undesider case that can come up:
   * 1. A fine grid vertex is not contained in the corase cell.
   *    To deal we this we express the fine grid DOF wrt to the DOFs "closest" to the corresponsing vertex.
   * 
   * 2. In the case a coarse grid cell/DOF is not contained in the fine cell, this introduces a zero column
   *    in the prolongation matrix. Use the nonNestedProlognationStack() function to remove zero columns from
   *    a resulting prolognation stack.
   * 
   * 
   */
  template <class Space>
  std::vector<Prolongation> makeNonNestedPstack(std::initializer_list<Space> spacelist, double tolTrunc=0.1)
  {
    std::cout << "Build prolongation stack based on "  << spacelist.size() << " spaces (levels).\n";
    std::vector<Prolongation> ps;
    for(auto space = spacelist.begin(); space!= spacelist.end()-1; ++space)
    {
      // get prolongation matrix between each pair of spaces
      ps.push_back(nonNestedProlongation(*space, *(space+1), tolTrunc));
    }
    
    // put them together in a prolongation stack
    return nonNestedProlognationStack(ps);
  }
  
  
    /**
   * \ingroup fetransfer
   * \brief Computes a (joint) prolongation stack for a geometry defined by two grids.
   * 
   * \tparam Space1 FEFunctionSpace built on gridManager of first geometry
   * \tparam Space2 FEFunctionSpace built on gridManager of second geometry
   * 
   * \param space1  space of first geometry
   * \param space2  space of second geometry
   * 
   * This function computes a prolongation stack for a hierarchy based on two (hierarchical) spaces. 
   * At least one of the spaces should be hierarchical, i.e. have more than one level. If the number 
   * of levels of the two spaces differ, the hierarchy of the space with fewer levels is "augmeted" 
   * by identity prolongation matrices. For each of the spaces this function calls 
   * makeDeepProlongationStack(...) and concatenates the two prolongation stacks.
   * 
   * 
   */
  template <class Space1,class Space2>
  std::vector<Prolongation> makeJointPstack(Space1 space1, Space2 space2, int coarseLevel=0)
  {
    std::vector<Prolongation> Pstack;
    // concat prolongation for product space (joint geometry)
    std::cout << "grid1:\n";
    std::vector<Prolongation> Pstack1 =  makeDeepProlongationStack(space1,coarseLevel);
    std::cout << "grid2:\n";
    std::vector<Prolongation> Pstack2 =  makeDeepProlongationStack(space2,coarseLevel);
    assert(("At least one grid should be hierarchical.", Pstack1.size()+Pstack2.size()>0));
    // adapt number of levels by appending identity prolongation
    // to grid with fewer levels
    size_t n1 = space1.grid().leafGridView().size(space1.dim), n2 = space2.grid().size(space2.dim);
    if(Pstack1.size()<Pstack2.size())
    {
      while(Pstack1.size()<Pstack2.size())
        Pstack1.push_back(sparseUnitMatrix<double,1>(n1));
    } else if (Pstack1.size() > Pstack2.size())
    {
      while(Pstack1.size()>Pstack2.size())
        Pstack2.push_back(sparseUnitMatrix<double,1>(n2));
    }
    assert(("Number of levels differ on the grids.",Pstack1.size() == Pstack2.size()));
    for(int i=0; i<Pstack1.size(); ++i)
      Pstack.push_back(diagcat(Pstack1[i],Pstack2[i]));
    for(int i=0; i<Pstack.size(); ++i)
      std::cout << "prolongation matrix P(" << i << "," << i+1 << ") has dimensions " << Pstack[i].N() << "x" << Pstack[i].M() << ".\n";
    
    return Pstack;
  }
  
  // ----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup mg
   * \brief creates a semi-geometric multigrid preconditioner based on an auxiliary space
   * \tparam Space the finite element space for which the stiffness matrix A has been computed. 
   *               Currently this is restricted to  a P1 finite element space. 
   *               \todo: implement for general spaces
   * \tparam Entry the type of matrix entries in the stiffness matrix
   * 
   * \return a multigrid stack with prolongations and (Galerkin-projected) stiffness matrices
   */
  template <class Space, class Entry, class Index>
  MultiGridStack<Prolongation,Entry,Index>
  makeAuxiliarySpaceMultigridStack(Space const& space, NumaBCRSMatrix<Entry,Index> const& A, 
                                   double volumeRatio=100)
  {
    assert(A.N() == space.degreesOfFreedom());
    using Grid = typename Space::Grid;
    std::cout << "creating cover grid..."; std::cout.flush();
    GridManager<Grid> coverGridMan(createCoverGrid(space.gridView(),volumeRatio));
    std::cout << " done\n"; std::cout.flush();
    
    // Create a P1 finite element space
    H1Space<Grid> coverSpace(coverGridMan,coverGridMan.grid().leafGridView(),1);
    std::cout << "cover grid vertices: " << coverSpace.degreesOfFreedom() << "\n"; std::cout.flush();
    
    // Create prolongation from cover grid to original grid 
    using Matrix = NumaBCRSMatrix<Entry,Index>;
    std::cout << "creating prolongation..."; std::cout.flush();
    Prolongation P = twoGridProlongation(coverSpace,space);
    std::cout << "done\n"; std::cout.flush();
    std::cout << "creating conjugation..."; std::cout.flush();
    std::vector<Matrix> as{conjugation(P,A,false,true),A};
    std::cout << "done\n"; std::cout.flush();
    std::vector<Prolongation> ps{std::move(P)};

    // It may happen that some vertices of the cover grid do not interact at all with the original grid.
    // Then the Galerkin projection is only positive semidefinite - in particular hard for direct solvers.
    // We add a small positive value to the affected diagonal entries.
    std::cout << "making spd..."; std::cout.flush();
    double maxDiagonal = 0;
    for (Index i=0; i<as[0].N(); ++i)
      maxDiagonal = std::max(maxDiagonal,as[0][i][i].frobenius_norm());
    for (Index i=0; i<as[0].N(); ++i)
      if (as[0][i][i].frobenius_norm() < 1e-8*maxDiagonal)
        as[0][i][i] = 1e-8*maxDiagonal*unitMatrix<typename Entry::field_type,Entry::rows>();
    std::cout << "done\n"; std::cout.flush();
    
    std::cout << "creating multigridstack..."; std::cout.flush();
    
    return MultiGridStack<Prolongation,Entry,Index>(std::move(ps),std::move(as));
  }
}

#endif
