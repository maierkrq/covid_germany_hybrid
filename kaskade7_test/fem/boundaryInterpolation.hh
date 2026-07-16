/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2017-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef BOUNDARYINTERPOLATION_HH
#define BOUNDARYINTERPOLATION_HH

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include "dune/common/fvector.hh"
#include "dune/grid/common/boundarysegment.hh"
#include "dune/grid/uggrid.hh"

#include "fem/barycentric.hh"
#include "fem/fixdune.hh"
#include "fem/forEach.hh"
#include "fem/gridBasics.hh"
#include "fem/gridcombinatorics.hh"

#include "utilities/detailed_exception.hh"
#include "utilities/power.hh"

namespace Kaskade
{
  /**
   * \ingroup grid
   * \brief A class storing a boundary edge consisting of two vertex indices and associated boundary curve tangents.
   * 
   * The class invariant is that the first vertex index is less than the second one and that the tangent vectors
   * are normalized.
   */
  template <class Index, class Vector>
  class BoundaryEdge
  {
  public:
    using ctype = typename Vector::value_type;
    constexpr static int dim = Vector::dimension;

    /**
     * \brief Constructor taking vertex indices and tangent vectors.
     *
     * Conceptually, the tangent vectors point towards the opposing vertex, i.e. into
     * the "interior" of the edge. Geometrically, if \f$ x_1 \f$ and \f$ x_2 \f$ are
     * the vertices, it must hold that \f$ (x_2-x_1)^T t_1 > 0 \f$ and \f$ (x_1-x_2)^T t_2 > 0 \f$.
     */
    BoundaryEdge(Index v1, Index v2, Vector const& t1, Vector const& t2)
    : vidx{v1,v2}, t{normalize(t1),normalize(t2)}
    {
      if (vidx[0] > vidx[1])
      {
        std::swap(vidx[0],vidx[1]);
        std::swap(t[0],t[1]);
      }
    }
    
    /**
     * \brief Constructor taking only the vertex indices.
     *
     * The tangent vectors are undefined.
     */
    BoundaryEdge(Index v1, Index v2)
    : BoundaryEdge(v1,v2,Vector(1.0),Vector(1.0))
    {
    }
    
    BoundaryEdge()
    : BoundaryEdge(-1,-1)
    {
    }
    
    /**
     * \brief Comparison by lexicographic order of the vertex indices. 
     * 
     * Tangents are not considered for comparison.
     */
    bool operator <(BoundaryEdge const& e) const 
    {
      return vidx < e.vidx;
    }
    
    /**
     * \brief Comparison by vertex indices. 
     * 
     * Tangents are not considered for comparison.
     */
    bool operator ==(BoundaryEdge const& e) const 
    {
      return vidx == e.vidx;
    }
    
    std::array<Index,2> vertices() const 
    {
      return vidx;
    }
    
    Index vertex(int i) const 
    {
      assert(0<=i && i<=1);
      return vidx[i];
    }
    
    Vector tangent(int i) const 
    {
      assert(0<=i && i<=1);
      return t[i];
    }
    
    std::array<Vector,2> const& tangents() const 
    {
      return t;
    }
    
    /**
     * \brief Checks whether the edge is incident to the given vertex.
     */
    bool incident(Index v) const 
    {
      return v==vidx[0] || v==vidx[1];
    }

    /**
     * \brief Applies displacement to the tangent vectors.
     *
     * A deformation \f$ \varphi(x) = x+u(x) \f$ modifies the edge tangent vector \f$ t \f$ to
     * \f$ \varphi_x t = t + u_x t \f$.
     */
    void displace(Dune::FieldMatrix<ctype,dim,dim> const& du0,
                  Dune::FieldMatrix<ctype,dim,dim> const& du1)
    {
      t[0] = normalize(t[0]+du0*t[0]);
      t[1] = normalize(t[1]+du1*t[1]);
    }
    
  private:
    std::array<Index,2> vidx;
    std::array<Vector,2> t;
  };
  
  // -----------------------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------------

  namespace Simplex
  {
    template <class Scalar, int dim>
    using Vector = Dune::FieldVector<Scalar,dim>;

    template <class Scalar, int dim>
    using Barycentric = Dune::FieldVector<Scalar,dim+1>;

    template <class Scalar, int dim>
    using VertexPositions = std::array<Vector<Scalar,dim>,dim+1>;

    template <class Scalar, int dim>
    using EdgeDirections = std::array<BoundaryEdge<int,Vector<Scalar,dim>>,dim*(dim+1)/2>;

    template <int dim>
    using FaceFlags = std::array<bool,dim+1>;


    /**
     * \ingroup fem
     * \brief Computes a boundary displacement based on blended cubic rational Bezier curves.
     *
     * This allows to define a G1 continuous boundary surface. On each vertex \f$ x_i \f$ of
     * a single simplicial surface patch, i.e. a triangle in 3D and a line in 2D, a unit normal
     * vector \f$ n_i \f$ is given.
     * \param vertices the boundary simplex's vertices, as columns of this matrix
     * \param normals the unit normal vectors at the vertices, as columns of this matrix
     * \param b the barycentric coordinates of the point to interpolate, i.e. the weight of the
     *          corresponding vertices/normals. The entries in b shall sum up to 1.
     *
     * \return the displacement vector for the point b
     */
    template <class Scalar, int dim>
    Vector<Scalar,dim> cellInterpolation(VertexPositions<Scalar,dim> const& vertices,
                                         EdgeDirections<Scalar,dim> const& boundaryEdges,
                                         FaceFlags<dim> const& isBoundaryFace,
                                         Barycentric<Scalar,dim> b);

    template <class Scalar, int dim>
    Dune::FieldMatrix<Scalar,dim,dim+1> cellInterpolationDerivative(VertexPositions<Scalar,dim> const& vertices,
                                                                    EdgeDirections<Scalar,dim> const& boundaryEdges,
                                                                    FaceFlags<dim> const& isBoundaryFace,
                                                                    Barycentric<Scalar,dim> b);
  }




  
  /**
   * \ingroup grid
   * \brief A class storing algebraic and geometric information about a vertex-boundaryface incidence in simplicial grids.
   */
  template <class Index, class ctype, int dim>
  struct Corner
  {
    using Vector = Simplex::Vector<ctype,dim>;
    
    /// Indices of the other vertices of the face.
    std::array<Index,dim-1>  vidx;

    /// Positions of the other vertices.
    std::array<Vector,dim-1> vx;

    /// Face normal at the corner.
    Vector                   fn;

    /// Boundary segment index of the face.
    int                      bsi;
  };
  

  /**
   * \ingroup grid
   * \brief A class storing all the corners (vertex-boundaryface incidences) at a boundary vertex.
   */
  template <class Index, class ctype, int dim>
  struct BoundaryStar
  { 
    using Vector = Simplex::Vector<ctype,dim>;
    
    /// Center vertex index.
    Index                                vidx;

    /// Center vertex position.
    Vector                               vx;

    /// Corners of the incident faces.
    std::vector<Corner<Index,ctype,dim>> corners;
  };


  /**
   * \brief Computes the boundary stars (i.e. corners around a vertex) for all boundary vertices of the grid view.
   *
   * \param gv the grid view to use, usually the level 0 coarse grid view.
   *
   * \return a map that maps indices of boundary vertices to their star
   */
  template <class GridView>
  std::map<typename GridView::IndexSet::IndexType,
           BoundaryStar<typename GridView::IndexSet::IndexType,typename GridView::ctype,GridView::dimension>>
  computeBoundaryStars(GridView const& gv)
  {
    constexpr int dim = GridView::dimension;
    using ctype = typename GridView::ctype;
    using Index = typename GridView::IndexSet::IndexType;
    using Vector = Simplex::Vector<ctype,dim>;


    auto const& is = gv.indexSet();

    // Edge tangent and surface normal information is stored for each boundary vertex.
    std::map<Index,BoundaryStar<Index,ctype,dim>> stars;

    // step through all boundary faces, adding their normals up at the vertices
    forEachBoundaryFace(gv,[&](auto const& face)
    {
      if (face.neighbor())    // on periodic boundary
        return;               // - skip

      Vector const n = face.centerUnitOuterNormal();              // get the face normal

      int vSubIdx[dim];                                           // the cell-local indices of vertices
      vertexids(dim,1,face.indexInInside(),vSubIdx);              // incident to the face
      Index vIdx[dim];

      for (int i=0; i<dim; ++i)                                   // for each vertex incident to the face
        vIdx[i] = is.subIndex(face.inside(),vSubIdx[i],dim);      // get the vertex id (later we need them in different order)

      for (int i=0; i<dim; ++i)                                   // for each vertex incident to the face
      {
        Corner<Index,ctype,dim> corner;                           // create the corner (i.e. face-vertex incidence)
        for (int j=0; j<dim-1; ++j)                               // with dim-1 (i.e. codim 1) edge tangent vectors
        {
          corner.vidx[j] = vIdx[(i+j+1)%dim];                     // target vertex of the edge and
          corner.vx[j] = face.inside().geometry().corner(vSubIdx[(i+j+1)%dim]);  // its position
        }
        corner.fn = n;                                            // as well as our face normal
        corner.bsi = face.boundarySegmentIndex();

        auto& star = stars[vIdx[i]];                              // append the generated corner to the list of
        star.vidx = vIdx[i];                                      // corners at this vertex (storing and overwriting
        star.vx = face.inside().geometry().corner(vSubIdx[i]);    // vertex index and position multiple times with the
        star.corners.push_back(corner);                           // same value).
      }
    });

    return stars;
  }

  // ----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup grid
   * \brief The interface for specifying feature edges and feature vertices.
   */
  template <class Index, class ctype, int dim>
  class GeometricFeatureDetector
  {
  public:
    using Vector = Simplex::Vector<ctype,dim>;

    virtual ~GeometricFeatureDetector() {}

    /**
     * \brief Decides whether two boundary edges emanating from a joint vertex of a 2D domain
     * are part of a smooth boundary curve or not.
     *
     * This method is called only when 2D domains are processed.
     *
     * \param vidx the vertex index in the level 0 coarse grid view
     * \param t1 tangent vector of first edge
     * \param t2 tangent vector of second edge
     * \param bsid1 the boundary segment index of the first edge
     * \param bsid2 the boundary segment index of the second edge
     *
     * The tangent vectors are normalized and point away from the vertex. Thus, t1+t2==0 means they
     * form a straight line.
     *
     */
    virtual bool isFeatureVertex(Index vidx, Vector const& t1, Vector const& t2, int bsid1, int bsid2) const = 0;

    /**
     * \brief Decides whether a boundary edge of a 3D domain is part of a smooth boundary patch or not.
     *
     * More precisely, it determines whether an edge-vertex incidence is a feature "edge". Edges can be
     * specified as being smooth at one of its vertices and sharp on the other.
     *
     * This method is called only when 3D domains are processed.
     *
     * \param vidx the vertex index in the level 0 coarse grid view
     * \param n1 the outer normal vector of the face on one side of the edge, evaluated at the vertex
     * \param n2 the outer normal vector of the face on the other side of the edge, evaluated at the vertex
     * \param bsid1 the boundary segment index of the first face
     * \param bsid2 the boundary segment index of the second face
     */
    virtual bool isFeatureEdge(Index vidx, Vector const& n1, Vector const& n2, int bsid1, int bsid2) const = 0;

    /**
     * \brief Decides whether two feature boundary edges of a 3D domain meeting in a joint vertex are part of
     *        a smooth boundary curve or not.
     *
     * This method is called only when 3D domains are processed.
     *
     * \param the vertex index in the level 0 coarse grid view
     * \param t1 tangent vector of first edge
     * \param t2 tangent vector of second edge
     *
     * The tangent vectors are normalized and point away from the vertex. Thus, t1+t2==0 means they
     * form a straight line.
     */
    virtual bool isFeatureCurve(Index vidx, Vector const& t1, Vector const& t2) const = 0;

    /**
     * \brief Provides a surface normal at non-feature vertices (where no feature edge passes through).
     *
     * \param vidx the vertex index
     * \param n a suggested normal, usually an average of incident face normals
     * \return the actual normal \f$ n_a \f$ to use. It holds \f$ n_a^T n > 0 \f$ and \f$ \|n_a\| > 0 \f$.
     *
     * The default implementation just passes the suggested normal on.
     */
    virtual Vector vertexNormal(Index vidx, Vector const& n) const
    {
      return n;
    }
  };


  /**
   * \ingroup grid
   * \brief Specifying geometric features in terms of angles between face normals and edge tangents.
   */
  template <class Index, class ctype, int dim>
  class AngleFeatureDetector: public GeometricFeatureDetector<Index,ctype,dim>
  {
  public:
    using Vector = Simplex::Vector<ctype,dim>;

    /**
     * \brief Constructor.
     * \param featureCurveAngle the angle (in radians) between two edges (feature edges in 3D, boundary edges in 2D)
       *                        that determines whether they form a differentiable
       *                        feature curve. Defaults to roughly 57°.
     * \param featureEdgeAngle the angle (in radians) between adjacent facets' normals that determines whether they share a
       *                       feature edge or a normal (\f$ G^1 \f$ continuous) edge. Defaults to roughly 29°. This is
       *                       only of effect in 3D grids.
     */
    AngleFeatureDetector(ctype featureCurveAngle=1.0, ctype featureEdgeAngle=0.5);

    virtual ~AngleFeatureDetector() {}

    /**
     * \brief Decides whether two boundary edges emanating from a joint vertex of a 2D domain
     * are part of a smooth boundary curve or not.
     *
     * This method is called only when 2D domains are processed.
     *
     * \param vidx the vertex index in the level 0 coarse grid view
     * \param t1 tangent vector of first edge
     * \param t2 tangent vector of second edge
     * \param bsid1 the boundary segment index of the first edge
     * \param bsid2 the boundary segment index of the second edge
     *
     * The tangent vectors are normalized and point away from the vertex. Thus, t1+t2==0 means they
     * form a straight line.
     */
    virtual bool isFeatureVertex(Index vidx, Vector const& t1, Vector const& t2, int bsid1, int bsid2) const;

    /**
     * \brief Decides whether a boundary edge of a 3D domain is part of a smooth boundary patch or not.
     *
     * More precisely, it determines whether an edge-vertex incidence is a feature "edge". Edges can be
     * specified as being smooth at one of its vertices and shart on the other.
     *
     * This method is called only when 3D domains are processed.
     *
     * \param vidx the vertex index in the level 0 coarse grid view
     * \param n1 the outer normal vector of the face on one side of the edge, evaluated at the vertex
     * \param n2 the outer normal vector of the face on the other side of the edge, evaluated at the vertex
     * \param bsid1 the boundary segment index of the first face
     * \param bsid2 the boundary segment index of the second face
     */
    virtual bool isFeatureEdge(Index vidx, Vector const& n1, Vector const& n2, int bsid1, int bsid2) const;

    /**
     * \brief Decides whether two feature boundary edges of a 3D domain meeting in a joint vertex are part of
     *        a smooth boundary curve or not.
     *
     * This method is called only when 3D domains are processed.
     *
     * \param the vertex index in the level 0 coarse grid view
     * \param t1 tangent vector of first edge
     * \param t2 tangent vector of second edge
     *
     * The tangent vectors are normalized and point away from the vertex. Thus, t1+t2==0 means they
     * form a straight line.
     */
    virtual bool isFeatureCurve(Index vidx, Vector const& t1, Vector const& t2) const;

  private:
    ctype featureCurveThreshold;
    ctype featureEdgeThreshold;
  };

  template <class GridView>
  using GridAngleFeatureDetector = AngleFeatureDetector<typename GridView::IndexSet::IndexType,
                                                        typename GridView::ctype,GridView::dimension>;


  /**
   * \brief Detects feature edges and adjusts the edge tangent vectors accordingly.
   *
   * Feature edges are determined according to two criteria:
   * - the normals \f$ n_1, n_2 \f$ of the incident faces deviate from each other significantly, i.e.
   *   \f$ n_1^T n_2 < \alpha \f$, where \f$ \alpha \f$ is the feature edge threshold
   * - the set \f$ \{p_1,p_2\} \f$ of incident faces' patch types occurs in the `featureEdges` list
   *   (the order of \f$ p_1, p_2 \f$ is not relevant).
   *
   * If any of these two criteria is satisfied, the edge is classified as feature edge.
   *
   * If exactly two feature edges meet in one vertex, and their normalized emanating tangent vectors \f$ t_a, t_b \f$ are
   * approximately parallel (i.e. \f$ \|t_a + t_b\| < \beta \f$, where \f$ \beta \f$ is the feature curve threshold),
   * these two tangent vectors are modified such as to be collinear.
   *
   * \param star vertices and faces around a boundary vertex 
   * \param[out] orientedEdges all boundary half-edges, indexed 
   */
  template <class Index, class ctype>
  void adjustEdgeTangents(BoundaryStar<Index,ctype,3>& star, std::map<std::array<Index,2>,Dune::FieldVector<ctype,3>>& orientedEdges,
                          GeometricFeatureDetector<Index,ctype,3> const& features);
  template <class Index, class ctype>
  void adjustEdgeTangents(BoundaryStar<Index,ctype,2>& star, std::map<std::array<Index,2>,Dune::FieldVector<ctype,2>>& orientedEdges,
                          GeometricFeatureDetector<Index,ctype,2> const& features);

  // -----------------------------------------------------------------------------------------------------
  
  /**
   * \ingroup grid
   * \brief A class providing information about domain boundary tangent vectors.
   * 
   * The surface normals are assumed to be continuous within each boundary face, but need not be
   * continuous across non-feature boundary edges. Such a piecewise \f$ G^1 \f$-continuous surface
   * is determined by face-wise interpolation, where for each boundary edge the corresponding
   * tangent vector of the smooth surface is provided.
   *
   * \todo consider to turn this into a function returning the edge set
   */
  template <int dim, class Index, class ctype = double>
  class BoundaryTangents
  {
  public:
    using Vector = Dune::FieldVector<ctype,dim>;
    
    /**
     * \brief Constructor 
     *
     * Feature edges are determined according to two criteria:
     * - the normals \f$ n_1, n_2 \f$ of the incident faces deviate from each other significantly, i.e.
     *   \f$ n_1^T n_2 < \alpha \f$, where \f$ \alpha \f$ is the feature edge threshold
     * - the set \f$ \{p_1,p_2\} \f$ of incident faces' patch types occurs in the `featureEdges` list
     *   (the order of \f$ p_1, p_2 \f$ is not relevant).
     *
     * If any of these two criteria is satisfied, the edge is classified as feature edge.
     *
     * If exactly two feature edges meet in one vertex, and their normalized emanating
     * tangent vectors \f$ t_a, t_b \f$ are approximately parallel (i.e.
     * \f$ \|t_a + t_b\| < \beta \f$, where \f$ \beta \f$ is the feature curve threshold),
     * these two tangent vectors are modified such as to be collinear.
     *
     *
     * \param grid
     * \param features an object telling whether an edge or a vertex are feature edges or vertices. The object
     *                 need only exist during the constructor call.
     */
    template <class Grid>
    BoundaryTangents(Grid const& grid, GeometricFeatureDetector<Index,ctype,dim> const& features)
    {
      static_assert(dim==Grid::dimension, "dimension mismatch.");

      using Vector = Simplex::Vector<typename Grid::ctype,dim>;
      
      // Edge tangent and surface normal information is stored for each boundary vertex.
      auto stars = computeBoundaryStars(grid.levelGridView(0));
      
      // Modify the edge tangents incident to the boundary vertices such that a piecewise G1-continuous surface 
      // is created. If an edge is no feature edge (determined by incident face normal angle), the corner surface 
      // normals coincide afterwards (and the edge tangents are adjusted accordingly). If exactly two feature edges
      // meet, and their tangent vectors are almost parallel, they are made parallel such as to form a G1-continuous
      // feature curve. The modified edge tangent vectors are stored in "orientedEdges", and the stars may be 
      // modified (but are no longer needed anyways).

      std::map<std::array<Index,2>,Vector> orientedEdges;
      for (auto& vertex_star: stars)
        adjustEdgeTangents(vertex_star.second,orientedEdges,features);
      
      // Now combine two oriented edges into one non-oriented edge.
      for (auto const& e: orientedEdges)
      {
        // find oriented edge with opposite direction
        auto const& t2 = orientedEdges[std::array<Index,2>{e.first[1],e.first[0]}];
        boundaryEdges2.insert(BoundaryEdge<Index,Vector>(e.first[0],e.first[1],e.second,t2));
      }
    }
    
  public:
    std::set<BoundaryEdge<Index,Vector>> boundaryEdges2;
  };
  
  /** 
   * \ingroup grid
   * \brief A convenient type alias for BoundaryTangents if the grid view type is known.
   * \relates BoundaryTangents
   */
  template <class GridView>
  using GridBoundaryTangents = BoundaryTangents<GridView::dimension,typename GridView::IndexSet::IndexType,typename GridView::ctype>;

  
  //---------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------
  
  namespace BoundaryInterpolationDetail
  {
    
    // provides an array containing the vertex positions as Dune FieldVectors of the corners of a given cell
    template <class Cell>
    auto getVertexPositions(Cell const& cell)
    {
      constexpr int dim = Cell::dimension;
      std::array<Dune::FieldVector<typename Cell::Geometry::ctype, dim>, dim+1> vertexPositions;
      for (int i=0; i<=dim; ++i)
        vertexPositions[i] = cell.geometry().corner(i);                       // corner position

      return vertexPositions;
    }

    // returns optional (only if given cell has any boundary edge) an array of edges of a given cell;
    // for boundary edges also tangents required for boundary interpolation are provided
    template <class Cell, class GridView, class Vector, class Index>
    auto getCellEdges(Cell const& cell, GridView const& gv,
                      std::set<BoundaryEdge<Index,Vector>> const& boundaryEdges,
                      std::array<Vector, Cell::dimension+1> const& vertexPositions)
    {
      constexpr int dim = Cell::dimension;
      using CellEdges = Simplex::EdgeDirections<typename Cell::Geometry::ctype,dim>;
      using ReturnType = std::optional<CellEdges>;

      // Get global vertex indices.
      std::array<Index, dim+1> vIdx;
      for (int i=0; i<=dim; ++i)
        vIdx[i] = gv.indexSet().subIndex(cell,i,dim);

      std::vector<BoundaryEdge<int,Vector>> cellEdges2;
      bool hasBoundaryEdges = false;
      for (int i=0; i<dim; ++i)
        for (int j=i+1; j<=dim; ++j)
        {
          BoundaryEdge<Index,Vector> edge(vIdx[i],vIdx[j],Vector(0.0),Vector(0.0));
          auto be = boundaryEdges.find(edge);
          if (be!=boundaryEdges.end())                                        // this one is a boundary edge
          {
            std::array<Vector,2> t = be->tangents();
            if (vIdx[i]>vIdx[j])
              std::swap(t[0],t[1]);

            cellEdges2.push_back(BoundaryEdge<int,Vector>(i,j,t[0],t[1]));
            hasBoundaryEdges = true;
          }
          else
          {
            Vector tx = normalize(vertexPositions[j]-vertexPositions[i]);
            cellEdges2.push_back(BoundaryEdge<int,Vector>(i,j,tx,-tx));
          }
        }

      if (!hasBoundaryEdges)                                                    // no boundary edges means the cell deformation
        return ReturnType();                                                    // is the identity -> do nothing

      CellEdges cellEdges;
      std::copy(begin(cellEdges2),end(cellEdges2),begin(cellEdges));
      return ReturnType(cellEdges);
    }

    // returns an array of bools; an entry indicates whether corresponding face lies on the boundary
    // (by local index of the opposite vertex, not by local face index)
    template <class Cell, class GridView>
    auto getIsBoundaryFace(Cell const& cell, GridView const& gv)
    {
      constexpr int dim = Cell::dimension;
      using ReturnType = std::array<bool,dim+1>;

      ReturnType isBoundaryFace;
      std::fill(begin(isBoundaryFace),end(isBoundaryFace),false);

      for (auto const& face: intersections(gv,cell))
        if (face.boundary() && !face.neighbor())
        {
          std::array<int,dim> vSubIdx;                                          // get the cell-local indices
          vertexids(dim,1,face.indexInInside(),&vSubIdx[0]);                    // of the face corners

          std::array<bool,dim+1> isOppositeBoundaryFace;
          std::fill(begin(isOppositeBoundaryFace),end(isOppositeBoundaryFace),true);
          for (int i=0; i<dim; ++i)
            isOppositeBoundaryFace[vSubIdx[i]] = false;

          for (int i=0; i<dim+1; ++i)
            isBoundaryFace[i] |= isOppositeBoundaryFace[i];
        }

      return isBoundaryFace;
    }
    
    
    template <class Cell, class GridView, class Vector, class Index>
    void getChildVertexDisplacements(Cell const& cell, GridView const& gv, 
                                     std::set<BoundaryEdge<Index,Vector>> const& boundaryEdges,
                                     std::vector<Vector>& xs,
                                     bool computeDisplacements)
    {
      constexpr int dim = Cell::dimension;
      
      auto const& gv0 = gv.grid().levelGridView(0);

      // find the vertices of the given coarse grid cell
      std::array<Vector,dim+1> vertexPositions = getVertexPositions(cell);
      
      // find all boundary Edges of the cell
      auto cellEdges = getCellEdges(cell, gv0, boundaryEdges, vertexPositions);
      if(!cellEdges)
        return;
      
      // find whether cell faces lie on boundary, faces are indexed by opposite local vertex index
      std::array<bool,dim+1> isBoundaryFace = getIsBoundaryFace(cell, gv0);
        
      // Now that we have set up the data describing the deformation of the given level 0 cell,
      // step through all child leaf cells and evaluate the vertex position displacement.
      for (auto const& child: descendantElements(cell,std::numeric_limits<int>::max()))    // step through all leaf child cells
        if (child.isLeaf())
          for (int i=0; i<=dim; ++i)                                            // and their corners
          {
            Vector const x = child.geometry().corner(i);                        // and get the global position  
            Index vidx = gv.indexSet().subIndex(child,i,dim);                   // and its global index 
            
            if ( (computeDisplacements && xs[vidx].two_norm2()>0)               // most vertices will be touched from several cells.
                || (!computeDisplacements && (xs[vidx]-x).two_norm2()>0) )      // since computing the cell interpolation is expensive,
              continue;                                                         // we check whether this has already been computed
            
            Vector const xi = cell.geometry().local(x);                         // as well as the local position in coarse grid cell
            auto b = barycentric2(xi);                                          // and the barycentric coordinates
            
            Vector u = Simplex::cellInterpolation(vertexPositions,*cellEdges,isBoundaryFace,b);
            xs[vidx] = computeDisplacements? u: x+u;
          }
    }
  }
  
  /**
   * \ingroup grid
   * \brief Computes a new position for each grid vertex based on boundary interpolation 
   * 
   * The grid has to contain only simplicial cells.
   * If a cell of the grid has more than one boundary face, these need to be separated by a feature
   * edge - otherwise the grid deformation will necessarily create degenerated flat cells on refinement.
   * Cells without boundary face shall have at most one boundary edge (only possible in 3D).
   * 
   * \param gv the grid view for which the vertex positions are computes, usually a leaf grid view
   * \param tangents provides the surface tangent vectors corresponding to boundary edges
   * \param computeDisplacements if true, only the vertex displacements are computed, otherwise their new
   *                             positions will be returned
   * 
   * \return a vector container of the displaced vertex positions (or the displacements,
   *         depending on the `computeDisplacements` parameter)
   */
  template <class GridView>
  std::vector<Dune::FieldVector<typename GridView::ctype,GridView::dimension>>
  interpolateVertexPositions(GridView const& gv, 
                             GridBoundaryTangents<GridView> const& tangents, 
                             bool computeDisplacements=false)
  {
    using namespace BoundaryInterpolationDetail;
    constexpr int dim = GridView::dimension;
    using Vector = Dune::FieldVector<typename GridView::ctype,dim>;
    
    // Before computing displacements of the (few) boundary-affected vertices, fill the positions with 
    // the original positions as default. This will be overwritten by displaced vertices.
    std::vector<Vector> xs(gv.size(dim));         // positions as many as vertices
    if (!computeDisplacements)
      for (auto const& v: vertices(gv))           // as default, enter the original positions for all vertices
        xs[gv.indexSet().index(v)] = v.geometry().center();    
    
    // Find all boundary edges of the coarse grid.
    auto const& gv0 = gv.grid().levelGridView(0);
    
    // Step through all coarse grid cells that are affected by the boundary displacement
    for (auto const& cell: elements(gv0))
      getChildVertexDisplacements(cell,gv,tangents.boundaryEdges2,xs,computeDisplacements);

    return xs;
  }

  
  /**
   * \ingroup grid
   * \brief The BoundaryInterpolationDisplacement class provides a Kaskade WeakFunctionView on the displacement
   *        from grid to actual domain boundary (which is computed by boundary interpolation).
   *
   * It does also implement the interface required for displacements to be used in the BoundaryLocator.
   */
  template <class GridView>
  class BoundaryInterpolationDisplacement
  {
    static constexpr int dim = GridView::dimensionworld;
    using ctype = typename GridView::ctype;
  public:
    using Displacement = GlobalPosition<GridView>;
    using DisplacementDerivative = Dune::FieldMatrix<ctype,dim,GridView::dimension>;

    BoundaryInterpolationDisplacement(BoundaryInterpolationDisplacement&&) = default;

    /**
     * \brief Constructor precomputes for each level-0-cell data that is needed for boundary interpolation
     * \param gv level 0 grid view (coarse level)
     * \param tangents tangents on grid boundary
     */
    BoundaryInterpolationDisplacement(GridView const& gv, GridBoundaryTangents<GridView> const& tangents)
    : gridView(gv)
    {
      using namespace BoundaryInterpolationDetail;

      int nCoarseCells = gridView.size(0);
      vertexPositions.resize(nCoarseCells);
      cellEdges.resize(nCoarseCells);
      isBoundaryFace.resize(nCoarseCells);

      for (auto const& cell: elements(gridView))
      {
        auto cellIdx = gridView.indexSet().index(cell);

        // find the vertices of the given coarse grid cell
        vertexPositions[cellIdx] = getVertexPositions(cell);

        // find all boundary Edges of the cell
        cellEdges[cellIdx] = getCellEdges(cell, gridView, tangents.boundaryEdges2, vertexPositions[cellIdx]);

        // find whether cell faces lie on boundary, faces are indexed by opposite local vertex index
        isBoundaryFace[cellIdx] = getIsBoundaryFace(cell, gridView);
      }
    }

    /**
     * \brief Provides a \f$ G^1 \f$-continuous boundary displacement of deformed grids.
     *
     * Let a (volume) displacement \f$ u \f$ and the corresponding deformation \f$ \varphi \f$ on a
     * simplicial grid with \f$ G^1 \f$-continuous boundary displacement \f$ v \f$ be given. In principle,
     * the composite boundary displacement would be \f[ x \mapsto \varphi(x+v(x))-x = u(x+v(x)) \approx u(x) + u_x(x)v(x). \f]
     * Unfortunately, \f$ \varphi \f$ is in general not \f$ G^1 \f$-continuous, such that the resulting composite boundary
     * deformation is *not* \f$ G^1 \f$-continuous.
     *
     * Therefore we define the composite boundary interpolation by interpolating almost from scratch. We take the
     * \f$ G^1 \f$-continuous boundary edge tangents from the reference grid (which means that feaure edges are
     * independent of the deformation) and rotate them according to  \f$ t \mapsto t + u_x t \f$.
     * Here, \f$ u_x \f$ at a vertex is computed by averaging the surrounding displacement derivatives, which yields a
     * \f$ G^1 \f$-continuous boundary interpolation.
     *
     * \warning As this works on the *coarse* grid and ignores the volume displacement between coarse grid vertices,
     *          this is mostly useful for P1 (if need be, P2) discretizations on a single level.
     */
    template <class VolumeDisplacement>
    void applyVolumeDisplacement(VolumeDisplacement const& u)
    {
      auto const& is = gridView.indexSet();
      std::vector<Dune::FieldMatrix<ctype,dim,dim>> du(is.size(dim),0.0);
      std::vector<int>                              count(du.size(),0);
      std::vector<std::array<size_t,dim+1>>         cellToGlobalVidx(is.size(0));


      // Evaluate the displacement at the coarse grid vertices and store the
      // deformed vertex position. Evaluate the displacement derivative as well,
      // and sum it up for later averaging.
      for (auto const& cell: elements(gridView))
      {
        auto cellIdx = is.index(cell);

        // find the vertices of the given coarse grid cell
        for (int i=0; i<=dim; ++i)
        {
          auto x = cell.geometry().corner(i);
          auto [leafCell,xi] = findCellOnLevel(cell,cell.geometry().local(x));
          vertexPositions[cellIdx][i] = x + u.value(leafCell,xi);    // displaced vertex position phi(x)

          auto vidx = is.subIndex(cell,i,dim);
          assert(0<=vidx && vidx<du.size());
          du[vidx] += u.derivative(leafCell,xi);
          ++count[vidx];

          cellToGlobalVidx[cellIdx][i] = vidx;
        }
      }

      // Compute averages of displacement derivative at coarse grid vertices
      for (size_t i=0; i<du.size(); ++i)
      {
        assert(count[i]>0);
        du[i] /= count[i];
      }

      // Turn al the available tangent vectors according to the averaged displacement derivative.
      // Remember that the vertex numbers in cellEdges are *local to the cell*, not global indices.
      for (size_t cellIdx=0; cellIdx<cellToGlobalVidx.size(); ++cellIdx)
        if (cellEdges[cellIdx])
          for (auto& ed: *cellEdges[cellIdx])
          {
            auto v0idx = cellToGlobalVidx[cellIdx][ed.vertex(0)];
            auto v1idx = cellToGlobalVidx[cellIdx][ed.vertex(1)];
            ed.displace(du[v0idx],du[v1idx]);
          }
    }


    /**
     * \brief Provides the displacement for a point given by cell and local position \f$ \xi \f$.
     *
     * \param cell any cell in the grid, can be leaf or any level below
     */
    Displacement value(Cell<GridView> const& cell, LocalPosition<GridView> const& xiFine) const
    {
      // Find coarse grid cell and local position in coarse cell. There are two ways to do this:
      // (i) going down step by step, looking up the position in the father cells on the way
      // (ii) looking up the global position in the coarse ancestor cell.
      // We take the first route here, because (in deformed grids) the global position need not
      // be contained in the coarse ancestor cell.
      Cell<GridView> coarseCell = cell;
      auto xi = xiFine;
      while (coarseCell.hasFather())
      {
        xi = coarseCell.geometryInFather().global(xi);
        coarseCell = coarseCell.father();
      }
      //      auto const xi = coarseCell.geometry().local(cell.geometry().global(xiFine)); // dangerous if cell is from deformed grid, i.e. global Position can lie outside of coarse Cell

      // We may have a volume displacement applied (implicitly assumed to be given on the vertices only
      // and piecewise linear). Thus, this has to be computed first, before we can add the boundary interpolation
      // correction.
      std::size_t coarseCellIdx = gridView.indexSet().index(coarseCell);
      auto b = barycentric2(xi);
      Displacement u(0.0);
      for (int i=0; i<=dim; ++i)
        u += b[i] * (vertexPositions[coarseCellIdx][i]-coarseCell.geometry().corner(i));
      if (cellEdges[coarseCellIdx])
        u += Simplex::cellInterpolation(vertexPositions[coarseCellIdx],*cellEdges[coarseCellIdx],
                                        isBoundaryFace[coarseCellIdx],b);
      return u;
    }

     /**
     * \brief Provides the global displacement derivative (w.r.t. global variable \f$ x \f$)
     *        for a point given by cell and local position  \f$ \xi \f$.
     *
     * \param cell any cell in the grid, can be leaf or any level below
     */
    DisplacementDerivative derivative(Cell<GridView> const& cell, LocalPosition<GridView> const& xiFine) const
    {
      Cell<GridView> coarseCell = cell;
      auto xi = xiFine;
      while (coarseCell.hasFather())
      {
        xi = coarseCell.geometryInFather().global(xi);
        coarseCell = coarseCell.father();
      }

      // Now we have the *local* position xi. But we need to compute the derivative w.r.t. the *global*
      // position x.
      auto dxi_dx = transpose(coarseCell.geometry().jacobianInverseTransposed(xi));

      std::size_t coarseCellIdx = gridView.indexSet().index(coarseCell);
      auto b = barycentric2(xi);
      auto db_dx = barycentric2Derivative(xi) * dxi_dx;

      Displacement u(0.0);
      DisplacementDerivative du_dx(0.0);

      for (int i=0; i<=dim; ++i)
      {
        u += b[i] * (vertexPositions[coarseCellIdx][i]-coarseCell.geometry().corner(i));
        du_dx += outerProduct(vertexPositions[coarseCellIdx][i]-coarseCell.geometry().corner(i),db_dx[i]);
      }

      if (cellEdges[coarseCellIdx])           // in case of constant (and, incidentally, zero) displacement, the derivative
        du_dx += Simplex::cellInterpolationDerivative(vertexPositions[coarseCellIdx],*cellEdges[coarseCellIdx],
                                                      isBoundaryFace[coarseCellIdx],b) * db_dx;

      return du_dx;
    }

  private:
    GridView const gridView;
    std::vector<Simplex::VertexPositions<ctype,dim>> vertexPositions;
    std::vector<std::optional<Simplex::EdgeDirections<ctype,dim>>> cellEdges; // TODO: optional wastes memory. use unique_ptr instead?
    std::vector<Simplex::FaceFlags<dim>> isBoundaryFace;
  };

  //---------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------

  /**
   * \ingroup grid
   * \brief Creates a composite boundary interpolation displacement from a volume displacement.
   * \relates BoundaryInterpolationDisplacement
   */
  template <class VolumeDisplacement>
  auto makeCompositeBoundaryInterpolation(VolumeDisplacement const& u, 
                                          GeometricFeatureDetector<typename VolumeDisplacement::Space::Grid::LevelGridView::IndexSet::IndexType,
                                                                   typename VolumeDisplacement::Space::Grid::ctype,
                                                                   VolumeDisplacement::Space::Grid::dimension> const& features)
  {
    using Grid = typename VolumeDisplacement::Space::Grid;
    using GridView = typename Grid::LevelGridView;

    GridBoundaryTangents<GridView> boundaryTangents(u.space().grid(),features);
    BoundaryInterpolationDisplacement<GridView> bid(u.space().grid().levelGridView(0),boundaryTangents);
    bid.applyVolumeDisplacement(u);
    return bid;
  }

  //---------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------

  template <class GridView>
  class BoundaryDisplacementByInterpolation
  {
    using Grid = typename GridView::Grid;
    using Index = typename Grid::LevelGridView::IndexSet::IndexType;
    using ctype = typename Grid::ctype;
    static int const dim = Grid::dimension;
    using Features = GeometricFeatureDetector<Index,ctype,dim>;
    
  public:
    using DisplacementDerivative = Dune::FieldMatrix<typename GridView::ctype,GridView::dimensionworld,GridView::dimension>;

    template <class VolumeDisplacement>
    BoundaryDisplacementByInterpolation(VolumeDisplacement const& displacement)
    {
      *this = displacement;
      features = std::make_unique<AngleFeatureDetector<Index,ctype,dim>>(1.0,0.5);
    }

    template <class VolumeDisplacement>
    BoundaryDisplacementByInterpolation(VolumeDisplacement const& displacement, std::unique_ptr<Features> features_)
    : features(std::move(features_))
    {
      *this = displacement;
    }



    template <class VolumeDisplacement>
    BoundaryDisplacementByInterpolation& operator=(VolumeDisplacement const& displacement)
    {
      u = std::make_unique<BoundaryInterpolationDisplacement<GridView>>(makeCompositeBoundaryInterpolation(displacement,*features));
      return *this;
    }

    int order(Cell<GridView> const&) const { return std::numeric_limits<int>::max(); }

    GlobalPosition<GridView> value(Cell<GridView> const& cell, LocalPosition<GridView> const& position) const
    {
      return u->value(cell,position);
    }

    DisplacementDerivative derivative(Cell<GridView> const& cell, LocalPosition<GridView> const& position) const
    {
      return u->derivative(cell,position);
    }

  private:
    std::unique_ptr<BoundaryInterpolationDisplacement<GridView>> u;
    std::unique_ptr<Features> features;
  };

}

#endif
