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

#ifndef LAGRANGESPACE_HH
#define LAGRANGESPACE_HH

#include <cmath>
#include <numeric>
#include <tuple>

#include "boost/multi_array.hpp"

#include "dune/grid/common/capabilities.hh"

#include "fem/converter.hh"
#include "fem/fixdune.hh"
#include "fem/gridcombinatorics.hh"
#include "fem/functionspace.hh"
#include "fem/lagrangeshapefunctions.hh"
#include "fem/scalarspace.hh"

/**
 * \file
 * \brief Lagrange Finite Elements
 * \author Martin Weiser
 *
 */

namespace Kaskade
{
  /**
   * \brief A local to global mapper for scalar Lagrange bases.
   *
   * This is to be used as a policy object for the UniformScalarMapper.
   * 
   * \tparam ScalarType The underlying field type (usually double).
   * \tparam GV The grid view on which the FE space is defined.
   * \tparam restricted 
   */
  template <class ScalarType, class GV, bool restricted=false>
  class LagrangeMapperImplementation
  {
  public:
    typedef ScalarType Scalar;

    typedef GV     GridView;
    typedef typename GridView::Grid     Grid;
    typedef typename GridView::IndexSet IndexSet;
    using Cell = Kaskade::Cell<GridView>;
    
    using ctype = typename Grid::ctype;
    
    static int const dim = Grid::dimension;

    typedef typename LagrangeShapeFunctionSetContainer<ctype,dim,Scalar,restricted>::value_type ShapeFunctionSet;

    /**
     * \brief This is called on grid modifications and can be overwritten if internal data needs to be
     *        recomputed on refinement or coarsening.
     *
     * The default implementation provided here does nothing.
     */
    void update()
    {}
    
    
    /**
     * \brief The index set obtained from gridView().
     */
    IndexSet const& indexSet() const { return gv.indexSet(); }

    /**
     * \brief The number of entities of given geometry type in our support.
     *
     * The default implementation just repoththrts the total number of entities in the grid view,
     * which corresponds to global support.
     */
    size_t size(Dune::GeometryType gt) const
    {
      return indexSet().size(gt);
    }

    /**
     * \brief Tells whether the given cell is contained in the support.
     *
     * The default implementation realizes global support, i.e., all cells are contained.
     */
    bool inSupport(Cell const& cell) const
    {
      return true;
    }

    /**
     * \brief Returns the index of the specified subentity.
     */
    typename IndexSet::IndexType subIndex(Cell const& cell, int subentity, int codim) const
    {
      return indexSet().subIndex(cell,subentity,codim);
    }
    

    GridView const& gridView() const { return gv; }

    /**
     * \param indexSet_ must contain a subset (or all) of the cells in the grid.
     * \param order_ .
     */
    LagrangeMapperImplementation(GridView const& gridView_, int order_)
    : ord(order_)
    , gv(gridView_)
    {
      simplexSFS = &lagrangeShapeFunctionSet<ctype,dim,Scalar,restricted>(Dune::GeometryType(Dune::GeometryType::simplex,dim),ord);
      
      // Cube shape functions are not yet defined for arbitrary order, thus obtaining them may throw. 
      // Ignore this, as the grid may not contain any cubes...
      try
      {
        cubeSFS = &lagrangeShapeFunctionSet<ctype,dim,Scalar,restricted>(Dune::GeometryType(Dune::GeometryType::cube,dim),ord);
      }
      catch(...) 
      {
        cubeSFS = nullptr;
      };
    }

    ShapeFunctionSet const& shapeFunctions(Cell const& cell) const
    {
      if (cell.type().isSimplex())
        return *simplexSFS;
      if (cell.type().isCube())
        if (cubeSFS)
          return *cubeSFS;
        else
          throw LookupException("No Lagrangian shape functions on cubes defined.\n",__FILE__,__LINE__);
      else
        throw LookupException("Not supported geometry type.\n",__FILE__,__LINE__);
    }

    ShapeFunctionSet const& lowerShapeFunctions(Cell const& cell) const
    {
      return lagrangeShapeFunctionSet<ctype,dim,Scalar,restricted>(cell.type(),ord-1);
    }
    
    ShapeFunctionSet const& emptyShapeFunctionSet() const
    {
      static EmptyShapeFunctionSet<ctype,dim,Scalar,ShapeFunctionSet::comps> empty;
      return empty;
    }

    int order() const { return ord; }


    typedef ScalarConverter<Cell,Scalar> Converter;

    /**
     * \brief A class implementing a matrix \f$ K \in \R^{n\times m}\f$ mapping a subset
     * of \f$ m \f$ global degrees of freedom (those given by globalIndices()) to
     * \f$ n \f$ local degrees of freedom (shape functions).
     *
     * For Lagrange elements, this usually realizes just the identity.
     * We also support the case that there are *no* degrees of freedom
     * associated to a cell, i.e. \f$ m=0 \f$, see ContinuousLagrangeMapperSubdomain.
     */
    class Combiner
    {
    public:
      template <class GlobalIndices>
      Combiner(GlobalIndices const& globalIndices, ShapeFunctionSet const& sfs)
      : m(globalIndices.size())
      , n(sfs.size())
      {
        // We support Lagrange spaces where cells either have full support (i.e., shape
        // functions correspond one-to-one to ansatz functions, n==m) or no support (i.e.,
        // there are no dofs at all, m==0, though shape functions may be defined, n>0).
        assert(m==0 || n==m);
      }

      /**
       * @brief In-place computation of \f$ A \leftarrow A K \f$.
       * Since \f$ K \f$ is the identity, this is a no-op.
       */
      template <class Matrix>
      void rightTransform(Matrix& A) const
      {
        assert(A.M()==n);
        if (m==0)
          A.resize(A.N(),0);
      }

      /// In-place computation of row vectors \f$ v \leftarrow v K \f$.
      template <int d, int k>
      void rightTransform(std::vector<VariationalArg<Scalar,d,k>>& v) const
      {
        assert(v.size()==n);
        if (m==0)
          v.clear();
      }

      /// In-place computation of \f$ A \leftarrow K^+ A \f$.
      template <class Matrix>
      void leftPseudoInverse(Matrix& A) const
      {
        assert(A.N()==n);
        if (m==0)
          A.resize(0,A.M());
      }

      /// Implicit conversion to a sparse matrix.
      /// This is just the identity.
      operator Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1>>() const
      {
        Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1>> K(m,m,Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >::random);
        for (int i=0; i<m; ++i)
          K.incrementrowsize(i);
        K.endrowsizes();
        for (int i=0; i<m; ++i)
          K.addindex(i,i);
        K.endindices();
        for (int i=0; i<m; ++i)
          *K[i].begin() = 1;
        return K;
      }

    private:
      int m;  // the number of ansatz functions (dofs) on the current cell
      int n;  // the number of shape functions on the current cell
    };

  private:
    int               ord;
    GridView          gv;
    ShapeFunctionSet const* simplexSFS;
    ShapeFunctionSet const* cubeSFS;
  };


  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


  template <class ScalarType, class GV>
  class DiscontinuousLagrangeMapperImplementation: public LagrangeMapperImplementation<ScalarType,GV>
  {
    typedef LagrangeMapperImplementation<ScalarType,GV> Base;
    typedef typename Base::Cell Cell;

  public:
    DiscontinuousLagrangeMapperImplementation(GV const& gridView, int order):
      Base(gridView,order)
    {}

    // Returns the number of degrees of freedom (global ansatz
    // functions) uniquely associated to the given subentity type.
    int dofOnEntity(Dune::GeometryType gt) const
    {
      int const ord = this->order();
      int const dim = GV::Grid::dimension;


      if (gt.dim() == dim) // dofs live only on codim 0 entities (cells).
      {
        if ( gt.isSimplex() )
        {
          if (dim==1) return ord+1;
          if (dim==2) return (ord+1)*(ord+2)/2;
          if (dim==3) return (ord+1)*(ord+2)*(ord+3)/6;
          assert("Unknown dimension"==0);
        }
        if ( gt.isCube() ) return pow(ord+1,dim);
        if ( gt.isPyramid() || gt.isPrism() || gt.isNone() ) assert( "Not implemented"==0);
      }

      // subentities with codim > 0 do not carry local degrees of freedom
      return 0;
    }

    // Returns the geometry type, subentity number in cell and subentity
    // codimension for the subentity to which the dof is
    // associated. Here we have a discontinuous discretization, where
    // all dofs are associated to the cells.
    template <class ShapeFunction, class Dummy>
    void entityIndex(Cell const& cell, ShapeFunction const& sf, int n,
                     Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Dummy&) const
    {
      gt = cell.type();
      subentity = 0;
      codim = 0;
      indexInSubentity = n;
    }
  };


  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for FEFunctionSpace s of piecewise polynomials of order Order.
   *
   * \tparam ScalarType scalar type
   * \tparam GV grid view
   */
  template <class ScalarType, class GV>
  class DiscontinuousLagrangeMapper:
  public UniformScalarMapper<DiscontinuousLagrangeMapperImplementation<ScalarType,GV> >
  {
    typedef DiscontinuousLagrangeMapperImplementation<ScalarType,GV> Implementation;
    typedef UniformScalarMapper<Implementation> Base;
    typedef DiscontinuousLagrangeMapper<ScalarType,GV> Self;

  public:
    typedef GV GridView;
    typedef LagrangeSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,ScalarType> ShapeFunctionSetImplementation;
    typedef int ConstructorArgument;
    static int const continuity = -1;

    /** \ingroup fem
     *  \brief Type of the FunctionSpaceElement, associated to the FEFunctionSpace
     * 
     * \tparam m number of components
     * 
     */
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };


    DiscontinuousLagrangeMapper(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {}
  };

  template <class ScalarType, class GV>
  class DiscontinuousLagrangeMapperSubdomain:
  public UniformScalarMapper<DiscontinuousLagrangeMapperImplementation<ScalarType,GV> >
  {
    typedef DiscontinuousLagrangeMapperImplementation<ScalarType,GV> Implementation;
    typedef UniformScalarMapper<Implementation> Base;
    typedef DiscontinuousLagrangeMapperSubdomain<ScalarType,GV> Self;

  public:
    typedef GV GridView;
    typedef LagrangeSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,ScalarType> ShapeFunctionSetImplementation;
    typedef int ConstructorArgument;
    static int const continuity = -1;

    /** \ingroup fem
     *  \brief Type of the FunctionSpaceElement, associated to the FEFunctionSpace
     * 
     * \tparam m number of components
     * 
     */
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<Self>,m> type;
    };


    DiscontinuousLagrangeMapperSubdomain(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {}
  };

  //---------------------------------------------------------------------

  template <class ScalarType, class GV, class ShapeFunctionFilter=ScalarSpaceDetail::AllShapeFunctions>
  class ContinuousLagrangeMapperImplementation
  : public LagrangeMapperImplementation<ScalarType,GV,std::is_same<ShapeFunctionFilter,ScalarSpaceDetail::RestrictToBoundary>::value>,
    public ShapeFunctionFilter
  {
    typedef LagrangeMapperImplementation<ScalarType,GV,std::is_same<ShapeFunctionFilter,ScalarSpaceDetail::RestrictToBoundary>::value> Base;

  public:
    typedef typename GV::IndexSet IndexSet;
    typedef ScalarType Scalar;
    typedef GV GridView;
    using Cell = typename Base::Cell;

    ContinuousLagrangeMapperImplementation(GV const& gridView, int order, ShapeFunctionFilter shapeFunctionFilter=ShapeFunctionFilter())
    : Base(gridView,order), ShapeFunctionFilter(shapeFunctionFilter)
    {
    }


    /**
     * \brief Returns the number of degrees of freedom (global ansatz
     *        functions) uniquely associated to the given subentity type.
     */
    int dofOnEntity(Dune::GeometryType gt) const
    {
      int const ord = this->order();

      if(gt.isSimplex())
      {
        if (gt.dim()==0) return 1;
        if (gt.dim()==1) return ord-1;
        if (gt.dim()==2) return std::max(0,(ord-2)*(ord-1)/2);
        if (gt.dim()==3) return std::max(0,(ord-3)*(ord-2)*(ord-1)/6);
        assert("Unknown dimension"==0);
      }
      if(gt.isCube()) return std::max(0.0,pow(ord-1,gt.dim()));
      if(gt.isPyramid() || gt.isPrism()) assert("Not  implemented"==0);

      assert("Unknown geometry type"==0);

      return -1;
    }

    /**
     * \brief Returns the geometry type, subentity number in cell and subentity
     *        codimension for the subentity to which the dof is associated.
     *
     * Here we have a continuous discretization, where all dofs are associated to
     * the entities on which the shape function node is located.
     *
     * \param[out] gt the geometry type of the subentity on which the shape function lies
     * \param[out] subentity the number of the subentity on which the shape function lies
     * \param[out] codim the codimension of the subentity on which the shape function lies
     * \param[out] indexInSubentity
     * \param[out] data
     */
    template <class ShapeFunction, class Data>
    void entityIndex(Cell const& cell, ShapeFunction const& sf, int n,
                     Dune::GeometryType& gt, int& subentity, int& codim, int& indexInSubentity, Data& data) const
    {
      int const dim = Cell::dimension;
      int const ord = this->order();
      IndexSet const& is = this->indexSet();

      int dummy;
      std::tie(dummy,codim,subentity,indexInSubentity) = sf.location();

      // Obtain geometry type of the subentity on which the shape function lives.
      // As we assume that all cells are either simplices or hexahedra, there is only
      // one geometry type per codimension.
      gt = is.types(codim)[0]; // limited to just one geometry type....

      // implementation below assumes simplicial cells
      assert(gt.isSimplex());

      // local index is globally unique in the interior of the cell and
      // on the vertices. Otherwise we need to compute a globally unique
      // numbering based on the global numbering of incident vertices.
      if (codim>0 && codim<dim && ord>2)
      {
        // Obtain the local barycentric index of the shape function.
        std::array<int,dim+1> xi = barycentric(SimplexLagrangeDetail::tupleIndex<dim>(ord,n),ord);

        int nVertices = dim+1-codim;

        // Obtain a globally unique sorting of the barycentric coordinates, and extract
        // the corresponding permutation and selection of the barycentric coordinates
        // in our subentity.
        GlobalBarycentricPermutation<dim> gbp(is,cell);
        int pi[nVertices];
        gbp.barycentricSubsetPermutation(subentity,codim,pi);
        
        // Since here the node is not located on a lower dimensional
        // subentity, we know that all barycentric indices are at least 1. Therefore we
        // may restrict ourselves to the interior nodes by subtracting 1
        // from all barycentric indices. Note that the "order" associated
        // with these interior nodes shrinks by 1+dimension(our subentity),
        // i.e. dim+1-codim !
        int idx[nVertices];
        for (int i=0; i<nVertices; ++i)
          idx[i] = xi[pi[i]]-1;

        // Finally we compute the local index of this permuted node within the 
        // interior of our subentity.
        indexInSubentity = SimplexLagrangeDetail::local(idx,dim-codim,ord-nVertices);
      }

      // consider restrictions to the boundary here
      this->treatBoundary(data,this->gridView(),cell,codim,subentity);
    }

  };

  /**
   * \ingroup fem
   * \brief A degrees of freedom manager for globally continuous FEFunctionSpace s of piecewise polynomials.
   *
   * \tparam ScalarType scalar type
   * \tparam GV grid view
   */
  template <class ScalarType, class GV>
  class ContinuousLagrangeMapper
  : public UniformScalarMapper<ContinuousLagrangeMapperImplementation<ScalarType,GV> >
  {
    typedef ContinuousLagrangeMapperImplementation<ScalarType,GV> Implementation;
    typedef UniformScalarMapper<Implementation> Base;

  public:
    typedef ScalarType Scalar;
    typedef GV GridView;
    typedef typename Base::ShapeFunctionSet ShapeFunctionSet;
    typedef LagrangeSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> ShapeFunctionSetImplementation;
    typedef int ConstructorArgument;
    static int const continuity = 0;

    /** \ingroup fem
     *  \brief Type of the FunctionSpaceElement, associated to the FEFunctionSpace
     * 
     * \tparam m number of components
     * 
     **/
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<ContinuousLagrangeMapper>,m> type;
    };

    /**
     * \brief Constructor.
     * \param gridView the grid view on which to define the space, usually a leaf grid view
     * \param order polynomial ansatz order of shape functions (> 0)
     */
    ContinuousLagrangeMapper(GridView const& gridView, int order):
      Base(Implementation(gridView,order))
    {
      assert(order >= 1);
    }
  };

  // --------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------

  template <class GV, class SupportOracle, class ScalarType=double>
  class ContinuousLagrangeMapperSubdomainImplementation
  : public ContinuousLagrangeMapperImplementation<ScalarType,GV>
  {
    using Base = ContinuousLagrangeMapperImplementation<ScalarType,GV>;
    using Index = typename GV::IndexSet::IndexType;

  public:
    using GridView = GV;

    ContinuousLagrangeMapperSubdomainImplementation(GridView const& gridView, int order,
                                                    SupportOracle&& supportOracle_)
    : Base(gridView,order)
    , supportOracle(std::move(supportOracle_))
    {}

    void update()
    {
      auto const& gv = this->gridView();
      auto const& is = gv.indexSet();
      int const dim = GV::dimension;

      // For each entity in the grid, create a char, which we can use to
      // flag whether this has already been seen or not.
      support.clear();

      for (int codim=0; codim<=dim; ++codim)              // For each entity codimension...
        for(auto const& geoType: is.types(codim))         // ... and each geometry type of that codimension...
        {
          auto& idx = support[geoType];                   // create the flag vector ...
          idx.resize(is.size(geoType),-1);                // ... with as many flags as there are entities of
        }                                                 // this type in the index set, and initialize them to -1.


      // Now step through all entities in the support. This is the closure of the cells within the support.
      // Thus, we step through all the cells and flag them as well as all their subentities.
      Dune::GeometryType cellGeometryType;
      for (auto const& cell : elements(gv))
      {
        cellGeometryType = cell.type();
        size_t const cellIndex = is.index(cell);

        if (supportOracle(cell))
        {
          support[cellGeometryType][cellIndex] = 0;
          auto refElem = referenceElement(cell.geometry());

          // Now step through all subentities (as we have processed the cell itself already, we
          // start with codimension 1). For all the subentities we obtain their geometry type
          // and index and flag them with 1 as within the support.
          for (int codim=1; codim<=dim; ++codim)
            for (int i=0; i<refElem.size(codim); ++i)
              support[refElem.type(i,codim)][is.subIndex(cell,i,codim)] = 0;
        }
      }

      // Now count how many entities of a given geometry type are within our support.
      supportSize.clear();
      for (auto& p: support)
      {
        // Assign global indices to subentities in the support by incrementing a counter.
        long idx = 0;
        for (size_t i=0; i<p.second.size(); ++i)
          if (p.second[i] >= 0)
            p.second[i] = idx++;

        // ... which in the end also yields the total number of entities of a given type.
        supportSize[p.first] = idx;
      }
    }

    /**
     * \brief Returns the number of entities of given geometry type in our support.
     */
    size_t size(Dune::GeometryType gt) const
    {
      return supportSize.find(gt)->second;
    }

    /**
     * \brief Tells whether the given cell is contained in the support.
     *
     * The default implementation realizes global support, i.e., all cells are contained.
     */
    bool inSupport(typename Base::Cell const& cell) const
    {
      auto it = support.find(cell.type());
      if (it==support.end())
        return false;
      return it->second[this->indexSet().index(cell)] >= 0;
    }

    Index subIndex(typename Base::Cell const& cell, int subentity, int codim) const
    {
      auto subentityGeometryType = referenceElement(cell).type(subentity,codim);
      auto it = support.find(subentityGeometryType);

      if (it==support.end())
        throw LookupException("Geometry type of provided cell not contained in the index set.",__FILE__,__LINE__);

      auto const& s = it->second;
      return s[this->indexSet().subIndex(cell,subentity,codim)];
    }

  private:
    SupportOracle supportOracle;
    std::map<Dune::GeometryType,size_t> supportSize;
    std::map<Dune::GeometryType,std::vector<long>> support;
  };


  /**
   * \ingroup fem
   * \brief A local-to-global mapper for continuous finite elements on a subdomain.
   *
   * \tparam GV the grid view on which the space is defined
   * \tparam SupportOracle a callable that for each cell specifies whether it is contained
   *                       in the subdomain (maximal support of the functions) or not
   * \tparam ScalarType
   *
   * Using this mapper is tricky, since the support oracle is usually defined as a lambda
   * function, the type of which has to be provided as a template parameter. Here, class
   * template deduction as of C++ 17 helps:
   * \code
   * FEFunctionSpace space(gridManager,ContinuousLagrangeMapperSubdomain(
   *                                        gridManager.grid().leafGridView(),order,
   *                                        [](auto& cell) { return cell.geometry().center()[0]>0; }));
   * \endcode
   */
  template <class GV, class SupportOracle, class ScalarType=double>
  class ContinuousLagrangeMapperSubdomain
  : public UniformScalarMapper<ContinuousLagrangeMapperSubdomainImplementation<GV,SupportOracle,ScalarType>>
  {
    typedef ContinuousLagrangeMapperSubdomainImplementation<GV,SupportOracle,ScalarType> Implementation;
    typedef UniformScalarMapper<Implementation> Base;

  public:
    typedef ScalarType Scalar;
    typedef GV GridView;
    typedef typename Base::ShapeFunctionSet ShapeFunctionSet;
    typedef LagrangeSimplexShapeFunctionSet<typename GridView::Grid::ctype,GridView::Grid::dimension,Scalar> ShapeFunctionSetImplementation;
    typedef int ConstructorArgument;

    /**
     * \brief Continuity of the functions in this space.
     * If the support is restricted to a proper subdomain, the functions are discontinuous across
     * the inner boundary - hence the low continuity even though inside the subdomain the functions
     * are \f$ C^0 \f$.
     */
    static int const continuity = -1;

    /** \ingroup fem
     *  \brief Type of the FunctionSpaceElement, associated to the FEFunctionSpace
     * 
     * \tparam m number of components
     * 
     **/
    template <int m>
    struct Element
    {
      typedef FunctionSpaceElement<FEFunctionSpace<ContinuousLagrangeMapperSubdomain>,m> type;
    };

    /**
     * \brief Constructor.
     * \param gridView the grid view on which to define the space, usually a leaf grid view
     * \param order polynomial ansatz order of shape functions (> 0)
     * \param supportOracle a callable that for any cell
     */
    ContinuousLagrangeMapperSubdomain(GridView const& gridView, int order, SupportOracle&& supportOracle)
    : Base(Implementation(gridView,order,std::move(supportOracle)))
    {
      assert(order >= 1);
    }
  };
}

#endif
