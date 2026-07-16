/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MG_PROLONGATION_HH
#define MG_PROLONGATION_HH

#include <boost/fusion/include/vector.hpp>
#include <boost/timer/timer.hpp>

#include "fem/barycentric.hh"
// #include "fem/deforminggridmanager.hh"
#include "fem/spaces.hh"
#include "linalg/dynamicMatrix.hh"
#include "linalg/localMatrices.hh"
#include "linalg/threadedMatrix.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"

namespace Kaskade
{
  /**
   * \ingroup multigrid
   * \brief Computes an interpolation-based prolongation matrix from a (supposedly) coarser space to a finer space.
   *
   * \tparam CoarseSpace a FEFunctionSpace type
   * \tparam FineSpace a FEFunctionSpace type
   *
   * \param coarseSpace the (supposedly) coarser space (domain)
   * \param fineSpace   the (supposedly) finer space (image)
   *
   * Both function spaces have to be defined on the very same grid view and to be geometrically nested.
   */
  template <class SparseIndex=size_t, class CoarseSpace, class FineSpace>
  NumaBCRSMatrix<Dune::FieldMatrix<typename FineSpace::Scalar,1,1>,SparseIndex> prolongation(CoarseSpace const& coarseSpace, 
                                                                                             FineSpace const& fineSpace)
  {
    namespace bf = boost::fusion;

    Timings& timer = Timings::instance();
    timer.start("space prolongation");

    std::vector<bf::vector<std::vector<size_t>,std::vector<size_t>,
                           DynamicMatrix<Dune::FieldMatrix<typename FineSpace::Scalar,1,1>>>> interpolationData(fineSpace.gridView().size(0));

    // Intermediate shape function values, declared here to prevent frequent reallocations
    typename CoarseSpace::Mapper::ShapeFunctionSet::SfValueArray afValues;

    // Step through all the cells and perform local interpolation on each cell.
    // TODO: do this in parallel
    for (auto const& cell: elements(fineSpace.gridView()))
    {
      auto index = fineSpace.indexSet().index(cell);

      // copy fine and coarse indices on cell - the mapper returns only views
      // TODO: store views in the interpolationData itself
      auto const& fineIndices = fineSpace.mapper().globalIndices(index);
      auto const& coarseIndices = coarseSpace.mapper().globalIndices(index);
      std::copy(fineIndices.begin(),  fineIndices.end(),  std::back_inserter(bf::at_c<0>(interpolationData[index])));
      std::copy(coarseIndices.begin(),coarseIndices.end(),std::back_inserter(bf::at_c<1>(interpolationData[index])));

      auto const& fineSfs = fineSpace.mapper().shapefunctions(cell);
      auto const& coarseSfs = coarseSpace.mapper().shapefunctions(cell);

      // Obtain interpolation nodes of target on this cell
      auto const& iNodes(fineSpace.mapper().shapefunctions(cell).interpolationNodes());

      // Evaluate coarse space global shape functions
      evaluateGlobalShapeFunctions(coarseSpace,cell,iNodes,afValues,coarseSfs);

      // Interpolate fine space global values
      approximateGlobalValues(fineSpace,cell,afValues,bf::at_c<2>(interpolationData[index]),fineSfs);
    }
    timer.stop("space prolongation");

    timer.start("matrix creation");
    // Create a sparse prolongation matrix sparsity pattern. Note that the local prolongation matrices
    // are often sparse (e.g. for Lagrangian elements). We filter out the zero entries in order not to
    // create too densely populated prolongations, which in turn lead to very expensive conjugation
    // computations.
    NumaCRSPatternCreator<SparseIndex> creator(fineSpace.degreesOfFreedom(),coarseSpace.degreesOfFreedom());
    for (auto const& block: interpolationData)
      for (int i=0; i<bf::at_c<0>(block).size(); ++i)
        for (int j=0; j<bf::at_c<1>(block).size(); ++j)
        {
          auto entry = bf::at_c<2>(block)[i][j];
          if (std::abs(entry) > 1e-12)
            creator.addElement(bf::at_c<0>(block)[i],bf::at_c<1>(block)[j]);
        }

    // Now create and fill the matrix itself. Overwrite already existing entries
    // (this means the weighing is just "last one wins". If the coarse space is 
    // a subspace of the fine space (which is usually the case), the prolongation 
    // is nevertheless exact. 
    NumaBCRSMatrix<Dune::FieldMatrix<typename FineSpace::Scalar,1,1>,SparseIndex> p(creator);
    for (auto const& block: interpolationData)
      for (int i=0; i<bf::at_c<0>(block).size(); ++i)
        for (int j=0; j<bf::at_c<1>(block).size(); ++j)
        {
          auto entry = bf::at_c<2>(block)[i][j];
          if (std::abs(entry) > 1e-12)
            p[bf::at_c<0>(block)[i]][bf::at_c<1>(block)[j]] = entry;
        }
    
    timer.stop("matrix creation");
    return p;
  }

  // ---------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Computes a stack of prolongation matrices for higher order finite element spaces.
   *
   * This computes a stack of prolongation matrices for a scale of finite element spaces with
   * increasing polynomial ansatz order. From level to level, the ansatz order is doubled (maybe +1
   * for odd degrees).
   *
   * \tparam Mapper a finite element local to global mapper with scalar shape functions.
   */
  template <class SparseIndex, class Mapper>
  std::vector<NumaBCRSMatrix<Dune::FieldMatrix<typename Mapper::Scalar,1,1>,SparseIndex>> 
  prolongationStack(FEFunctionSpace<Mapper> const& space)
  {
    std::vector<NumaBCRSMatrix<Dune::FieldMatrix<typename Mapper::Scalar,1,1>,SparseIndex>> stack;
    int p = space.mapper().maxOrder();

    if (p > 1) // there is some need for prolongation
    {
      H1Space<typename Mapper::Grid> coarseSpace(space.gridManager(),space.gridView(),p/2);  // a coarser space (lower order)
      stack = prolongationStack(coarseSpace);
      stack.push_back(prolongation<SparseIndex>(coarseSpace,space));
    }

    return stack;
  }

  // ---------------------------------------------------------------------------------------------------------

  /// \cond internals
  namespace ProlongationDetail
  {
    struct Node
    {
      size_t p, q; // global (fine grid) indices of parent vertices
      int level;   // level of the vertex
    };
  }
  /// \endcond

  /**
   * \ingroup multigrid
   * \brief A prolongation operator for P1 finite elements from a coarser grid level to the next finer level.
   *
   * For P1 finite elements on simplicial grids, the prolongation from a coarser level to the next finer is a matrix \f$ P \f$
   * in which each row contains either one entry of value 1 (if the vertex is already contained in the coarser level grid)
   * or two entries that sum up to 1 (if the vertex is created on the finer level by bisecting the edge between two parent
   * nodes). In the latter case, the two entries default to 0.5 each.
   *
   * Storing this as a sparse matrix data structure is quite inefficient. This class is a specialized implementation 
   * of such prolongation matrices.
   */
  class MGProlongation
  {
  public:
    /**
     * \brief Constructor.
     * \param parents a sequence of node information. The parent indices for fine grid nodes shall refer to the fine grid numbering.
     * \param indexInCoarse a vector of length of fine grid vertices. For each coarse grid node, it contains the index of that node
     *                      in the coarse grid.
     * \param nc the number of coarse grid nodes (i.e. number of columns in P)
     * \param fineLevel the level of fine grid nodes
     */
    MGProlongation(std::vector<ProlongationDetail::Node> const& parents, std::vector<size_t> const& indexInCoarse,
                   size_t nc, int fineLevel);

    /**
     * \brief Returns the parent nodes indices.
     *
     * For a fine grid node, this returns the parent node indices in the coarse grid. For a coarse
     * grid node, the "parent indices" coincide and refer to the node index in the coarse grid.
     */
    std::array<size_t,2> const& parents(size_t i) const
    {
      return entries[i];
    }

    /**
     * \brief Matrix vector multiplcation (update mode).
     *
     * This computes \f$ f \leftarrow f + Pc \f$.
     *
     * \param f the fine grid target vector. It has to have the correct size.
     */
    template <class Vector>
    void umv(Vector const& c, Vector& f) const
    {
      assert(f.size()==entries.size());
      assert(c.size()==nc);

      for (size_t row=0; row<entries.size(); ++row)
      {
        auto e = entries[row];
        f[row] += static_cast<typename Vector::field_type>(0.5)*(c[e[0]] + c[e[1]]);
      }
    }

    /**
     * \brief Matrix vector multiplcation.
     *
     * This computes \f$ f \leftarrow Pc \f$, and is faster than but equivalent to
     * \code
     * f = 0;
     * umv(c,f);
     * \endcode
     *
     * \param f the fine grid target vector. It has to have the correct size.
     */
    template <class Vector>
    void mv(Vector const& c, Vector& f) const
    {
      assert(f.size()==entries.size());
      assert(c.size()==nc);

      for (size_t row=0; row<entries.size(); ++row)
      {
        auto e = entries[row];
        f[row] = static_cast<typename Vector::field_type>(0.5)*(c[e[0]] + c[e[1]]);
      }
    }

    /**
     * \brief Transpose matrix vector multiplication
     *
     * This computes \f$ c \leftarrow P^T f \f$.
     *
     * \param c the coarse grid target vector. It is resized to contain the result.
     */
    template <class Vector>
    void mtv(Vector const& f, Vector& c) const
    {
      assert(f.size()==entries.size());
      if (c.size()!=nc)
        c.resize(nc);
      for (size_t i=0; i<c.N(); ++i)
        c[i] = 0;

      for (size_t row=0; row<entries.size(); ++row)
      {
        c[entries[row][0]] += static_cast<typename Vector::field_type>(0.5)*f[row];
        c[entries[row][1]] += static_cast<typename Vector::field_type>(0.5)*f[row];
      }
    }

    /**
     * \brief The number of rows.
     */
    size_t N() const
    {
      return entries.size();
    }

    /**
     * \brief The number of columns.
     */
    size_t M() const
    {
      return nc;
    }

    /**
     * \brief Returns the prolongation in form of a sparse matrix.
     */
    template <class Real>
    NumaBCRSMatrix<Dune::FieldMatrix<Real,1,1>> asMatrix() const 
    {
      NumaCRSPatternCreator<> creator(N(),M(),false,2);
      
      for (size_t i=0; i<N(); ++i)
        creator.addElements(&i,&i+1,begin(entries[i]),end(entries[i]));
      
      NumaBCRSMatrix<Dune::FieldMatrix<Real,1,1>> P(creator);
      for (size_t i=0; i<N(); ++i)
        for (size_t j: entries[i])
          P[i][j] += 0.5;
        
      return P;
    }

    /**
     * \brief Galerkin projection of fine grid matrices.
     * This computes \f$ P^T A P \f$.
     *
     * \tparam Entry the entry type of Galerkin matrix to be projected, in general a quadratic Dune::FieldMatrix type
     * \tparam Index the Index type of the matrix to be projected
     *
     * \param A the quadratic fine level Galerkin matrix of size \f$ N\times N \f$ to be projected
     * \param onlyLowerTriangle if true, only the lower triangular part of symmetric A is touched, and only the lower triangular part of \f$ P^T A P\f$ is created
     */
    template <class Entry, class Index>
    NumaBCRSMatrix<Entry,Index> galerkinProjection(NumaBCRSMatrix<Entry,Index> const& A, bool onlyLowerTriangle = false) const
    {
      assert(A.N()==N() && A.M()==N());
      Timings& timer = Timings::instance();

      // First, create the sparsity pattern of the projected matrix.
      NumaCRSPatternCreator<Index> creator(M(),M(),onlyLowerTriangle);

      // If C = P^T A P, we have that C_{ij} = \sum_{k,l} P_{ki} P_{lj} A_{kl}. Hence, the entry A_{kl} contributes to
      // all C_{ij} for which there are nonzero entries P_{ki} and P_{lj} in the rows k and l of P. Thus we can simply
      // run through all nonzeros A_{kl} of A, look up the column indices i,j of rows k and l of P, and flag C_{ij}
      // as nonzero.

      // Step through all entries of A
      timer.start("h conjugation pattern");
      for (Index k=0; k<A.N(); ++k)
      {
        auto const& is = entries[k];                        // indices i for which Pki != 0
        auto row = A[k];
        for (auto ca=row.begin(); ca!=row.end(); ++ca)
        {
          Index const l = ca.index();
          auto const& js = entries[l];                      // indices j for which Plj != 0
          // add all combinations i,j, but be careful no index must occur twice in any range
          if (is[0]==is[1] && js[0]==js[1])
            creator.addElement(is[0],js[0]);
          else if (is[0]==is[1])
            creator.addElements(std::begin(is),std::begin(is)+1,std::begin(js),std::end(js));
          else if (js[0]==js[1])
            creator.addElements(std::begin(is),std::end(is),std::begin(js),std::begin(js)+1);
          else
            creator.addElements(std::begin(is),std::end(is),std::begin(js),std::end(js));

          if (onlyLowerTriangle && k>l)                     // subdiagonal entry (k,l) of A -> entry (l,k) must be treated implicitly:
            creator.addElements(std::begin(js),std::end(js),std::begin(is),std::end(is)); // add all combinations (j,i)
        }
      }
      timer.stop("h conjugation pattern");

      // An alternative way of computing the sparsity pattern would be to use that the nonzero entries j in column i
      // of C are exactly those for which there is k with (nonzero P_{jk} and there is l with (nonzero P_{il} and A_{lk})).
      // Hence we can obtain the column index set J directly by the following steps:
      // (i) find all l with P_{li} nonzero -> L  [requires to access columns of P - compute the transpose patterns once]
      // (ii) find all k with A_{lk} nonzero for some l in L -> K  [probably sorting K and removing doubled entries would be a good idea here]
      // (iii) find all j with P_{kj} nonzero for some k in K -> J
      // Compared to the above implementation this would have the following (dis)advantages
      // + easy to do in parallel (since write operations are separated)
      // + fewer scattered write accesses to memory
      // - more complex implementation
      // - requires the transpose pattern of P

      // Create the sparse matrix.
      timer.start("matrix creation");
      NumaBCRSMatrix<Entry,Index> pap(creator);
      timer.stop("matrix creation");

      // Fill the sparse matrix PAP. This is done as before by stepping through all Akl entries and scatter
      // Pki*Plj*Akl into PAPij.
      timer.start("matrix conjugation");
      for (Index k=0; k<A.N(); ++k)
      {
        auto const& is = entries[k];
        auto row = A[k];
        for (auto ca=row.begin(); ca!=row.end(); ++ca)
        {
          Index const l = ca.index();
          auto const& js = entries[l];

          if (is[0]==is[1] && js[0]==js[1])               // a coarse grid node - this means four identical contributions
            pap[is[0]][js[0]] += *ca;                     // with factor 1/4. Substitute with one contribution with factor 1.
          else
          {
            auto pap0 = pap[is[0]];
            static typename Entry::field_type const quarter = 0.25;
            pap0[js[0]] += quarter * *ca;
            pap0[js[1]] += quarter * *ca;

            auto pap1 = pap[is[1]];
            pap1[js[0]] += quarter * *ca;
            pap1[js[1]] += quarter * *ca;
          }

          if (onlyLowerTriangle)
            abort();  // not yet implemented
        }
      }
      timer.stop("matrix conjugation");

      return pap;
    }
    

  private:
    // Internally we treat all rows equally: A row with one entry of value one is represented as
    // two entries which sum up to 1 (that happen to reference the same column index...). This allows a
    // very simple and uniform implementation. The vector entries contains the column indices of
    // the two entries in each row.
    std::vector<std::array<size_t,2>> entries;

    // Number of columns (i.e. coarse grid nodes).
    size_t nc;
  };

  std::ostream& operator<<(std::ostream& out, MGProlongation const& p);

  /**
   * \ingroup multigrid
   * \brief Creates a Galerkin projected Matrix \f$ P^T A P \f$ from a prolongation \f$ P \f$ and a
   *        symmetric matrix \f$ A \f$.
   * \param onlyLowerTriangle if true, A contains onl the lower triangular part of the symmetric matrix.
   */
  template <class Entry, class Index>
  auto conjugation(MGProlongation const& p, NumaBCRSMatrix<Entry,Index> const& a, bool onlyLowerTriangle=false)
  {
    return p.template galerkinProjection(a,onlyLowerTriangle);
  }

  // ---------------------------------------------------------------------------------------------------------

  /// \cond internals
  namespace ProlongationDetail
  {
    // For each corner of the given cell, if there are any of them which are not corners of the father cell
    // (i.e. created by bisection of an edge), enter the parent vertices into the parents vector.
    template <class LeafView, class Cell, class Parents>
    void computeParents(LeafView const& leafView, Cell const& cell, Parents& parents)
    {
      // If this cell is a coarse grid cell, none of its vertices have parents, and we're done.
      if (!cell.hasFather())
        return;

      int const dim = LeafView::dimension;

      // For each corner, check its position in the parent. If it has barycentric coordinates with one
      // entry approximately one (all others zero), it is a corner of the father cell and will be treated
      // later (or has already been treated). Otherwise there will be two entries 0.5 (all others zero),
      // and these entries denote father corners which are the parents.
      assert(cell.type().isSimplex());
      auto const& geo = cell.geometryInFather();
      int nCorners = geo.corners();
      for (int i=0; i<nCorners; ++i)
      {
        // TODO: We could first obtain the index of the corner and check whether its parents have
        //       already been determined, and skip the geometric considerations in that case.
        //       This should save some time.
        //       On the other hand, obtaining indices can also be expensive, and in the current
        //       implementation we only have to get them for corners that actually are no corners
        //       in the father cell. This saves some time, too.
        //       We should check which option is faster, and whether it makes a big difference in
        //       the first place.

        // barycentric coordinates of corner in the father cell
        auto b = barycentric(geo.corner(i));

        // find potential parent vertices (those with barycentric coordinates 0.5)
        int pcount = 0;
        int pc[2];
        for (int k=0; k<b.N(); ++k)                         // check all barycentric coordinates
          if (std::abs(b[k]-0.5) < 0.01)                    // if close to 0.5 accept
          {
            pc[pcount] = (k+1) % (dim+1);                   // map barycentric coordinate number to Dune corner number
            ++pcount;                                       // note down that we've found one more parent vertex
          }

        if (pcount == 2)                                    //  corner i is a father edge midpoint
        {
          auto const& is = leafView.indexSet();
          parents.push_back(std::make_pair(is.index(cell.template subEntity<dim>(i)),
                                           Node{is.index(cell.father().template subEntity<dim>(pc[0])),
                                                is.index(cell.father().template subEntity<dim>(pc[1])),
                                                cell.level()}));
        }
      }
    }
    
    // For each corner of the given cell, if there are any of them which are not corners of the father cell
    // (i.e. created by bisection of an edge), enter the parent vertices into the parents vector.
    template <class LeafView, class Cell, class Parents>
    void computeShiftedParents(LeafView const& leafView, Cell const& cell, Parents& parents, 
                               bool interpolateAtShiftedPosition=true)
    {
      // If this cell is a coarse grid cell, none of its vertices have parents, and we're done.
      if (!cell.hasFather())
        return;

      int const dim = LeafView::dimension;

      // For each corner, check its position in the parent. If it has barycentric coordinates with one
      // entry approximately one (all others zero), it is a corner of the father cell and will be treated
      // later (or has already been treated). Otherwise there will be two entries 0.5 (all others zero),
      // and these entries denote father corners which are the direct parents.
      assert(cell.type().isSimplex());
      auto const& geo = cell.geometryInFather();
      int nCorners = geo.corners();
      for (int i=0; i<nCorners; ++i)
      {
        // barycentric coordinates of corner in the father cell
        auto b = barycentric(geo.corner(i));

        // find potential parent vertices (those with barycentric coordinates 0.5)
        int pcount = 0;
        int pc[2];
        for (int k=0; k<b.N(); ++k)                         // check all barycentric coordinates
          if (std::abs(b[k]-0.5) < 0.01)                    // if close to 0.5 accept
          {
            pc[pcount] = (k+1) % (dim+1);                   // map barycentric coordinate number to Dune corner number
            ++pcount;                                       // note down that we've found one more parent vertex
          }

        if (pcount == 2)                                    //  corner i is a father edge midpoint
        {
          auto x = cell.geometry().corner(i);               // global (possibly displaced) vertex position
          auto y0 = cell.father().geometry().corner(pc[0]); // global parent vertex position
          auto y1 = cell.father().geometry().corner(pc[1]); // global parent vertex position
          
          auto const& is = leafView.indexSet();
          std::vector<std::pair<size_t,double>> parentVertices;

          if (!interpolateAtShiftedPosition ||                        // we ignore shifting vertices or
              (0.5*(y0+y1)-x).two_norm() < 1e-5*(y0-y1).two_norm())   // this vertex appears not to be shifted
          {
            parentVertices.push_back(std::make_pair(is.index(cell.father().template subEntity<dim>(pc[0])),0.5));
            parentVertices.push_back(std::make_pair(is.index(cell.father().template subEntity<dim>(pc[1])),0.5));
          }
          else                                                        // this is shifted and we want to respect it
          {
            // compute linear interpolation from father's corners to the current vertex
            // TODO: consider performing the interpolation via the reference element
            Dune::FieldMatrix<double,dim+1,dim+1> A;
            for (int j=0; j<=dim; ++j)
            {
              auto vj = cell.father().geometry().corner(j);
              for (int k=0; k<dim; ++k)
                A[k][j] = vj[k];
              A[dim][j] = 1;
            }
            Dune::FieldVector<double,dim+1> r,w;
            for (int k=0; k<dim; ++k)
              r[k] = x[k];
            r[dim] = 1;
            A.solve(w,r);
            
            for (int j=0; j<=dim; ++j)
              parentVertices.push_back(std::make_pair(is.index(cell.father().template subEntity<dim>(j)),w[j]));
          }
          
          // Now that we have the contributions, insert them in the list.
          parents.push_back(std::make_tuple(static_cast<size_t>(is.index(cell.template subEntity<dim>(i))), 
                                            cell.level(),parentVertices));
        }
      }
    }

    /**
     * \ingroup multigrid
     * \brief Creates a stack of prolongations from parent-child relationships in grids
     * \param level
     */
    std::vector<MGProlongation> makeProlongationStack(std::vector<Node> parents, int maxLevel, int minLevel, size_t minNodes);
    
    std::vector<NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>>
    makeDeformedProlongationStack(std::vector<std::pair<int,std::vector<std::pair<size_t,double>>>>,
                                  int maxLevel, int minLevel, size_t minNodes);
  };
  /// \endcond

  // ---------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Computes a sequence of prolongation matrices for P1 finite elements in hierarchical grids.
   *
   * We assume that each vertex not contained in the coarse grid has been created by bisecting an edge,
   * such that there are exactly two "parent vertices".
   * \tparam Grid the Dune grid type on which the P1 space is defined. The grid has to be a simplicial grid.
   * \param grid the grid itself
   * \param minNodes minimum number of nodes to keep as coarse grid
   */
  template <class GridMan>
  std::vector<MGProlongation> prolongationStack(GridMan const& gridman, int coarseLevel=0, size_t minNodes=0)
  {
    Timings& timer = Timings::instance();

    auto const& grid = gridman.grid();
    auto const& leafView = grid.leafGridView();
    auto const& cellRanges = gridman.cellRanges(grid.levelGridView(0));

    // In parallel compute the parent nodes for edge midpoints cell by cell
    timer.start("parent computation");
    std::vector<std::vector<std::pair<size_t,ProlongationDetail::Node>>> myParents(cellRanges.maxRanges());
    parallelFor([&](int k, int n)
    {
      for (auto const& coarseCell: cellRanges.range(n,k))
        for (auto const& cell: descendantElements(coarseCell,grid.maxLevel()))
          ProlongationDetail::computeParents(leafView,cell,myParents[k]);
    },myParents.size());

    // Now that all parent nodes have been identified, consolidate them in an easily
    // indexable array.
    std::vector<ProlongationDetail::Node> parents(leafView.size(leafView.dimension),ProlongationDetail::Node{0,0,-1});
    for (auto const& ps: myParents)
      for (auto const& p: ps)
        parents[p.first] = p.second;
    timer.stop("parent computation");


    // Now construct the prolongation matrices for all levels.
    timer.start("made stack parents");
    auto ps = ProlongationDetail::makeProlongationStack(std::move(parents),grid.maxLevel(),coarseLevel,minNodes);
    timer.stop("made stack parents");
    return ps;
  }

  // ---------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Computes a sequence of prolongation matrices for P1 finite elements in hierarchical grids.
   *
   * We assume that each vertex not contained in the coarse grid has been created by bisecting an edge,
   * and possibly shifted a moderate distance.
   * 
   * If a vertex is shifted, interpolating their direct parent nodes leads to nonlinear interpolation, i.e. 
   * coarse grid functions which are not piecwise linear on the (deformed) coarse grid cells. This can lead 
   * to high energy of the coarse grid cells. We therefore interpolate the vertex value not only from its 
   * direct parents, but from all corners of a father cell.
   * 
   * \tparam Grid the Dune grid type on which the P1 space is defined. The grid has to be a simplicial grid.
   * \param grid the grid itself
   * \param coarseLevel the grid level to use as coarsest level
   * \param minNodes minimum number of nodes to keep as coarse grid
   * \param interpolateAtShiftedPosition whether or not potential vertex shifts shall be respected or not
   */
  template <class GridMan>
  std::vector<NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>> 
  deformedProlongationStack(GridMan const& gridman, int coarseLevel=0, size_t minNodes=0,
                            bool interpolateAtShiftedPosition=true)
  {
    using std::get;
    using InterpolationNode = std::pair<size_t,double>;           // (index,weight)
    using InterpolationNodes = std::vector<InterpolationNode>;    // complete sequence of interpolation nodes
    
    Timings& timer = Timings::instance();
    
    auto const& grid = gridman.grid();
    auto const& leafView = grid.leafGridView();
    auto const& cellRanges = gridman.cellRanges(grid.levelGridView(0));
    constexpr int dim = GridMan::Grid::dimension;
    
    // Edge midpoints may be shifted on refinement by the deforming grid manager if they are located on the boundary.
    // In that case, linear interpolation by the two parent nodes does not provide a good prolongation. The reason
    // is that low energy modes that can exactly be represented on the (undeformed) coarse grid are distorted and tend to have
    // much higher energy. Consequently, the coarse grid correction cannot take large steps. Here we try to take
    // the node displacement into account and provide a linear interpolation from the the corners of nearby cells.
    // If the shifted edge midpoint is contained in a cell, then we take this cell. Otherwise we average over all 
    // cells that are incident to the edge (in 2D this is exactly one).
    std::vector<std::vector<std::tuple<size_t,int,InterpolationNodes>>> interpolationData(cellRanges.maxRanges());
    parallelFor([&](int k, int n)
    {
      for (auto const& coarseCell: gridman.cellRanges(grid.levelGridView(0)).range(n,k))
        for (auto const& cell: descendantElements(coarseCell,grid.maxLevel()))
          ProlongationDetail::computeShiftedParents(leafView,cell,interpolationData[k],interpolateAtShiftedPosition);
    },interpolationData.size());
      
    // Now that all possible vertex interpolations have been computed, gather them for each vertex.
    std::vector<std::vector<std::pair<int,InterpolationNodes>>> interpolationOptions(leafView.size(dim));
    for (auto const& idk: interpolationData)
      for (auto const& id: idk)
        interpolationOptions[get<0>(id)].push_back(std::make_pair(get<1>(id),get<2>(id)));
    
    // TODO: check whether clearing interpolationData here improves performance. Reasoning: the memory may be hot in
    //       the cache, an reusing it for subsequent allocations might be faster
    
    // Vertices can be placed on the father edge midpoint, then all interpolation options should contain just two 
    // parent cells. Otherwis, there can be at most one father cell containing the vertex. If there is one, we use 
    // this convex combination interpolation. If there is none, we average all extrapolations.
    std::vector<std::pair<int,InterpolationNodes>> interpolation(interpolationOptions.size());
    for (size_t k=0; k<interpolationOptions.size(); ++k)
    {
      auto const& i = interpolationOptions[k];
      
      if (i.empty())  // that's a coarse grid vertex - no interpolation from parents
        continue;

      // Check that all father cells agree whether a vertex is shifted or not.
      bool isEdgeMidpoint = false;
      bool isShifted = false;
      for (auto const& ip: i)
      {
        assert(ip.first == i[0].first); // all levels are equal
        isEdgeMidpoint = isEdgeMidpoint || ip.second.size()==2;
        isShifted = isShifted || ip.second.size()>2;
      }
      assert(isEdgeMidpoint != isShifted);

      if (isEdgeMidpoint)                     // Standard case: vertex is just the midpoint of its 
      {                                       // father edge. Then all interpolation options coincide:
        interpolation[k] = i[0];              // the vertex value is just the mean of the parent node 
        continue;                             // values. We're done.
      }
      
      for (auto const& ip: i)                 // Vertex is shifted. Check whether there is a father cell 
      {                                       // that contains the shifted position. Inside means the 
        bool inner = true;                    // barycentric interpolation weights are all nonnegative.
        for (auto const& p: ip.second)        // If there is a containing father cell, we interpolate 
          inner = inner && p.second>=0;       // from its corners, omitting other father cells which
        if (inner)                            // might need extrapolation.
          interpolation[k] = ip;
      }
      
      if (interpolation[k].second.size()>0)   // we found a father cell containing the shifted vertex
        continue;                             // and inserted the interpolation weights -> done

        
      interpolation[k].first = i[0].first;    // level is the same for all 
      for (auto const& ip: i)                 // and collect all contributing coarse nodes, with plain averaging
        for (auto const& cn: ip.second)       // of all interpolation options
          interpolation[k].second.push_back(std::make_pair(cn.first,cn.second/i.size()));
    }
    
    // TODO: use move semantics for interpolation instead of copy
    timer.start("made stack parents");
    auto ps = ProlongationDetail::makeDeformedProlongationStack(interpolation,grid.maxLevel(),coarseLevel,minNodes);
    timer.stop("made stack parents");
    return ps;
  }
  
  // ---------------------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief Class for multigrid stacks.
   *
   * This provides the storage of and access to prolongations and projected Galerkin matrices,
   * but leaves the construction of these to derived classes.
   */
  template <class Prolongation, class Entry, class Index>
  class MultiGridStack
  {
  public:

    /**
     * \brief Constructor
     *
     * This takes both a stack of prolongations and a matching stack of projected. 
     * \param ps a stack of prolongations 
     * \param as a stack of stiffness matrices. as.size()==ps.size()+1 has to hold.
     * 
     * Usually, it holds that \f$ A_i = P_i^T A_{i+1} P_i \f$, i.e. as[0] contains the coarse grid matrix.
     */
    MultiGridStack(std::vector<Prolongation>&& ps, std::vector<NumaBCRSMatrix<Entry,Index>>&& as)
    : prolongations(std::move(ps)), galerkinMatrices(std::move(as))
    {
      assert(prolongations.size()+1==galerkinMatrices.size());
    }

    /**
     * \brief Constructor
     *
     * This takes a stack of prolongations and creates the stack of projected Galerkin matrices from the given
     * fine grid matrix. Most useful for geometric multigrid, where the prolongations are defined solely in terms of the grid.
     *
     * \see makeMultiGridStack
     */
    MultiGridStack(std::vector<Prolongation>&& ps, NumaBCRSMatrix<Entry,Index>&& A, bool onlyLowerTriangle)
    {
      Timings& timer = Timings::instance();

      prolongations = std::move(ps);
      galerkinMatrices.push_back(std::move(A));

      assert((prolongations.size() == 0) || (prolongations.back().N() == galerkinMatrices[0].N()));
      assert(galerkinMatrices[0].N() == galerkinMatrices[0].M());

      // Step through the prolongations from top level to bottom
      timer.start("matrix projection");
      for (auto pi=prolongations.crbegin(); pi!=prolongations.crend(); ++pi)
        galerkinMatrices.insert(galerkinMatrices.begin(),conjugation(*pi,galerkinMatrices.front(),onlyLowerTriangle));
      timer.stop("matrix projection");
    }

    MultiGridStack(MultiGridStack&& other) = default;

    /**
     * \brief The number of grid levels
     */
    int levels() const
    {
      return galerkinMatrices.size();
    }

    /**
     * \brief Returns the prolongation from given level to next higher one.
     * \param level precondition 0 <= level < levels()-1
     */
    Prolongation const& p(int level) const
    {
      assert(0<=level && level<levels()-1);
      return prolongations[level];
    }

    /**
     * \brief Returns the projected Galerkin matrix on the given level.
     * \param level precondition 0 <= level < levels()
     * 
     * For level==0, the returned matrix can be in an undefined state if it has been modified previously
     * via the non-const coarseGridMatrix() reference.
     * 
     */
    NumaBCRSMatrix<Entry,Index> const& a(int level) const
    {
      assert(0<=level && level<levels());
      return galerkinMatrices[level];
    }

    /**
     * \brief Returns the projected Galerkin matrix on the coarsest level.
     *
     * This is explicitly a mutable reference, such that the coarse grid matrix can be moved
     * from. This is useful, as in the multigrid, the coarsest level matrix is not referenced
     * (the coarse grid preconditioner has its own copy, maybe obtained by moving from here...):
     * \code
     *  auto coarseSolver = makeDirectPreconditioner(std::move(mgStack.coarseGridMatrix()));
     * \endcode
     */
    NumaBCRSMatrix<Entry,Index>& coarseGridMatrix()
    {
      return galerkinMatrices[0];
    }

    void report(std::ostream& out) const
    {
      for (auto const& p: prolongations)
        out << "Prolongation:\n" << p;

      for (auto const& a: galerkinMatrices)
        out << "GalerkinMatrix: \n" << a;
    }

  private:
    std::vector<Prolongation>                 prolongations;    // prolongations grid i -> i+1, i=0,...,n-1
    std::vector<NumaBCRSMatrix<Entry,Index>>  galerkinMatrices; // Galerkin matrices on grid i, i=0,...,n
  };

  template <typename Prolongations, typename Entry, typename Index>
  std::ostream& operator<<(std::ostream& out, MultiGridStack<Prolongations,Entry,Index> const& mgStack) { mgStack.report(out); return out; }

  /**
   * \ingroup multigrid
   * \brief Convenience routine for creating multigrid stacks
   *
   * Given a stack of prolongations and the top level Galerkin matrix, this creates the complete stack
   * including all projected Galerkin matrices.
   *
   * \param ps a vector of prolongation matrices
   * \param A the top level (finest) Galerkin matrix to be projected down
   * \param onlyLowerTriangle if true, only the lower triangular part of A will be referenced
   * 
   * \relates MultiGridStack
   */
  template <class Prolongation, class Entry, class Index>
  MultiGridStack<Prolongation,Entry,Index> makeMultiGridStack(std::vector<Prolongation>&& ps, NumaBCRSMatrix<Entry,Index>&& A,
                                                              bool onlyLowerTriangle)
  {
    return MultiGridStack<Prolongation,Entry,Index>(std::move(ps),std::move(A),onlyLowerTriangle);
  }

  /**
   * \ingroup multigrid
   * \brief convenience routine for creating multigrid stacks based on geometric coarsening for P1 elements
   */
  template <class GridMan, class Entry, class Index>
  MultiGridStack<MGProlongation,Entry,Index> makeGeometricMultiGridStack(GridMan const& gridManager, NumaBCRSMatrix<Entry,Index>&& A,
                                                                         size_t minNodes=10000, bool onlyLowerTriangle=false)
  {
    return MultiGridStack<MGProlongation,Entry,Index>(prolongationStack(gridManager,0,minNodes),std::move(A),onlyLowerTriangle);
  }

  /**
   * \ingroup multigrid
   * \brief Creates stack of prolongations and projected Galerkin matrices.
   *
   * \param A the symmetric sparse matrix
   * \param n stop the coarsening if the number of rows/cols drops below this number
   * \param onlyLowerTriangle if true, only the lower triangular part of symmetric A is accessed
   * 
   * The prolongation/restriction used is particularly simple: Stepping through the nodes, every not yet classified node 
   * is classified as coarse node, and its direct neighbours as fine nodes. Thus every fine node has at least one coarse neighbour,
   * and the coarse nodes are not adjacent. 
   * 
   * For each fine node, the at most two neighbouring coarse nodes with largest connection strength are selected as parents. 
   * From these, the fine node is interpolated.
   * 
   * \todo Currently, the interpolation is just with factor 0.5 from both coarse parents - not really clever for unstructured grids 
   *       or jumping coefficients. Implement a more clever weighting.
   */
  template <class Entry, class Index>
  MultiGridStack<MGProlongation,Entry,Index> makeAlgebraicMultigridStack(NumaBCRSMatrix<Entry,Index>&& A,
                                                                         Index n=0, bool onlyLowerTriangle=false);

  
  /**
   * \ingroup multigrid
   * \brief Convenience routine for creating a deep prolongation stack from coarse to fine grid and on from P1 to higher order.
   * 
   * This relies on a geometric grid hierarchy for the P1 substack.
   * 
   * The multigrid recursion towards coarser levels is terminated as soon as the given coarse level is reached
   * or if the given number of nodes is no longer reached.
   */
  template <typename FineSpace>
  auto makeDeepProlongationStack(FineSpace const& space, int coarseLevel=0, int minNodes=0)
  {
    using Real = typename FineSpace::Scalar;
    using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<Real,1,1>>;
    
    // Create the P1 h-prolongation stack. We use sparse matrices as prolongation format, since this is
    // compatible with the p-prolongation below.
    
    auto prolongations = deformedProlongationStack(space.gridManager(),coarseLevel,minNodes);
    
    // if the space is higher order, add another prolongation layer from P1 to higher order
    if (space.mapper().maxOrder() > 1)
    {
      H1Space<typename FineSpace::Grid,Real> coarseSpace(space.gridManager(),space.gridView(),1);  // a coarser space (lower order)
      prolongations.push_back(prolongation<typename Matrix::size_type>(coarseSpace,space));
    }
    
    return prolongations;
  }
  
  /**
   * \ingroup multigrid
   * \brief convenience routine for creating multigrid stacks based on coarsening by reducing the ansatz order from P to P1
   */
  template <typename FineSpace, typename Matrix>
  auto makePMultiGridStack(FineSpace const& space, Matrix&& A, bool onlyLowerTriangle)
  {
    std::vector<NumaBCRSMatrix<Dune::FieldMatrix<typename FineSpace::Scalar,1,1>, typename Matrix::size_type>> prolongations;

    if (space.mapper().maxOrder() > 1)
    {
      H1Space<typename FineSpace::Grid, typename FineSpace::Scalar> coarseSpace(space.gridManager(),space.gridView(),1);  // a coarser space (lower order)
      prolongations.push_back(prolongation<typename Matrix::size_type>(coarseSpace,space));
    }

    return makeMultiGridStack(std::move(prolongations),std::move(A),onlyLowerTriangle);
  }
  
  /// \cond internals
  namespace ProlongationDetail
  {
    template <typename FineSpace, typename CoarseSpace, typename Matrix>
    auto makePMultiGridStack(FineSpace const& fineSpace, Matrix&& fA, 
                             CoarseSpace const& coarseSpace, std::unique_ptr<Matrix> cA = std::unique_ptr<Matrix>())
    {
      using Prolongation = NumaBCRSMatrix<Dune::FieldMatrix<typename FineSpace::Scalar,1,1>, typename Matrix::size_type>;
      std::vector<Prolongation> prolongations;
      prolongations.push_back(prolongation<typename Matrix::size_type>(coarseSpace,fineSpace));
    
      std::vector<Matrix> matrices;
      if (cA)
        matrices.push_back(std::move(*cA));
      else
        matrices.push_back(conjugation(prolongations[0],fA));
      matrices.push_back(std::move(fA));
    
      return MultiGridStack<Prolongation,typename Matrix::block_type,
                                         typename Matrix::size_type>(std::move(prolongations),std::move(matrices));
    }
  }
  /// \endcond
  
  
  /**
   * \ingroup multigrid 
   * \brief convenience routine for creating multigrid stacks between two spaces.
   */
  template <typename FineSpace, typename CoarseSpace, typename Matrix>
  auto makePMultiGridStack(FineSpace const& fineSpace, Matrix&& fA, CoarseSpace const& coarseSpace)
  {
    return ProlongationDetail::makePMultiGridStack(fineSpace,std::move(fA),coarseSpace);
  }
  
  /**
   * \ingroup multigrid
   * \brief convenience routine for creating multigrid stacks between two spaces.
   *
   * An approximation of the projected Galerkin matrix is provided. Often, this can be assembled much more
   * efficiently than projected.
   */
  template <typename FineSpace, typename CoarseSpace, typename Matrix>
  auto makePMultiGridStack(FineSpace const& fineSpace, Matrix&& fA, CoarseSpace const& coarseSpace, Matrix&& cA)
  {
    return ProlongationDetail::makePMultiGridStack(fineSpace,std::move(fA),coarseSpace,moveUnique(std::move(cA)));
  }

}

#endif
