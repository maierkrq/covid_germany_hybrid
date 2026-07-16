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

#ifndef CONTACTCONSTRAINTS
#define CONTACTCONSTRAINTS

#include <algorithm>
#include <limits>
#include <mutex>
#include <optional>

#include <dune/common/fvector.hh>
#include "dune/istl/bvector.hh"

#include "fem/boundaryLocator.hh"
#include "fem/gridBasics.hh"
#include "fem/quadrature.hh"
#include "linalg/threadedMatrix.hh"
#include "utilities/bezier.hh"
#include "utilities/threading.hh"

namespace Kaskade 
{
  /**
   * \ingroup contact
   * \brief Computes a single multibody contact point constraint.
   * 
   * For a given boundary point \f$ x \f$, the outer normal ray (in direction of the unit outer normal
   * vector \f$ n\f$) is tested for intersection with other boundary faces. More precisely, intersections
   * are considered on the ray \f$ \{ x+tn \mid t\in[-\mathrm{overlap},\infty\mathclose[\} \f$. Positive values
   * of overlap allow to detect contact partners in an infeasible, already overlapping configuration, which
   * is common for penalty methods or may arise due to nonlinearity.
   *
   * If an intersection is found at point \f$ y \f$, the (linearized) non-penetration condition for a displacement
   * \f$ \delta u \f$ is 
   * \f[ n^T (\delta u(x)-\delta u(y)) \le |x-y|. \f]
   * With a basis representation of \f$ \delta u = \sum_j \delta u_j \varphi_j\f$, this leads to a
   * scalar algebraic constraint of the form
   * \f[ \sum_j n^T(\varphi_j(x) -\varphi_j(y)) \delta u_j \le |x-y|. \f]
   * Let \f$ J(z) := \{ j\in \mathbf{N} \mid z \in \mathrm{supp}\varphi_j \} \f$ denote the dofs
   * associated to \f$ z \f$. The return value is the (optional) pair
   * \f[ \left( \{ (j,n^T(\varphi_j(x) -\varphi_j(y))) \mid j \in J(x)\cup J(y) \}, |x-y| \right), \f]
   * which can establish one row of a complete contact constraint of the form \f$ B\delta u \le g \f$.
   *
   * \return an optional pair (vector[(i,c)],g).
   *
   * In order to support nonlinear contact problems (with finite displacement), we consider an already 
   * deformed domain given in terms of a displacement \f$ u \f$, that is provided by the boundary locator.
   * This displacement \f$ u \f$ can, but need not be a finite element function.
   * 
   * \param space the FEFunctionSpace for which the constraints are to be computed
   * \param boundaryLocator a spatial index for looking up boundary intersections with rays
   * \param face the current boundary face
   * \param xi the local coordinate in the face
   * \param overlap states how far into the interior for an opposed face is searched.
   *                If `std::numeric_limits<double>::lowest()` is provided, a default value proportional to the
   *                cell diameter is chosen.
   */
  template <class Space, class DisplacedFace, class Displacement>
  auto getContactConstraint(Space const& space,
                            BoundaryLocator<typename Space::GridView,Displacement> const& boundaryLocator,
                            DisplacedFace const& face,
                            Dune::FieldVector<typename Space::GridView::ctype,Space::GridView::Grid::dimension-1> const& xi,
                            double overlap)
  {
    using Grid = typename Space::Grid;
    using ctype = typename Grid::ctype;
    using RowEntry = Dune::FieldMatrix<ctype,1,Grid::dimension>;
    using EntryType = std::pair<size_t,RowEntry>;
    using ReturnType = std::optional<std::pair<std::vector<EntryType>,ctype>>;

    // Find intersections of faces with the ray starting at the given point and extending in normal direction.
    // In a Newton loop, the current configuration need not be feasible (i.e. without overlap / interpenetration).
    // Hence we pull back a little bit (on the order of a single element diameter) to detect a slight overlap as well.
    auto n = face.unitOuterNormal(xi);


    if (overlap == std::numeric_limits<double>::lowest())
    {
      ctype diameter = std::pow(face.volume(),1.0/DisplacedFace::facedimension);
      overlap = 0.05*diameter;
    }

    auto const x = face.global(xi);
    auto intersection = boundaryLocator.byRay(x-overlap*n,n);

    if (!intersection)                  // no intersection found -> no constraint 
      return ReturnType();
    
    auto [oface,oxi,calpha] = *intersection;  // extract opposing boundary face, local point, and cos(alpha) of normals
    auto const y = oface.global(oxi);

    std::unique_lock lock(GridLocking<Grid>::mutex());
    auto dofs  = space.linearCombination( face.gridFace().inside(), face.gridFace().geometryInInside().global(xi));
    auto odofs = space.linearCombination(oface.gridFace().inside(),oface.gridFace().geometryInInside().global(oxi));
    lock.unlock();

    
    std::vector<EntryType> entries;

    RowEntry nt; nt[0] = n;             // normal as a row vector

    ctype const eps = std::numeric_limits<ctype>::epsilon();
    for (auto const& p: dofs)
      if (p.second.two_norm2() > 1e3*eps)                           // omit negligible contributions. Those are
        entries.push_back(std::make_pair(p.first,p.second[0]*nt));  // typically structurally zero, but have rounding errors
    for (auto const& p: odofs)
      if (p.second.two_norm2() > 1e3*eps)
        entries.push_back(std::make_pair(p.first,-p.second[0]*nt));
    
//  std::cout << "Constraint found:\n";
//  std::cout << "x          = " << x << "\n";
//  std::cout << "x opposite = " << y << "\n";
//  std::cout << "n          = " << n << "\n";
//  std::cout << "n opposite = " << oface.unitOuterNormal(oxi) << "\n";
//  std::cout << "cos(alpha) = " << calpha << "\n";
//  std::cout << "gap        = " << (y-x)*n << "\n";
//  if ((y-x)*n < 0)
//    abort();

    // compute the gap distance (can be negative in case of interpenetration)
    return ReturnType(std::make_pair(entries,(y-x)*n));
  }

  /**
   * \ingroup contact
   * \brief Computes a single multibody contact point constraint.
   *
   * This version is for problems where each body is represented by its own grid.
   *
   * \param space1 the first grid's displacement FE-space with respect to which the contact constraint is to be computed.
   * \param boundaryLocator1 boundary geometry for current point
   * \param space2 the second grid's displacement FE-space with respect to which the contact constraint is to be computed
   * \param boundaryLocator2 a spatial index for looking up boundary intersections with rays
   * \param face the current boundary face
   * \param xi the local coordinate in the face
   * \param overlap states how far into the interior for an opposed face is searched.
   * Negative value means search is only performed in the exterior with corresponding minimal distance.
   */
  template <class Space1, class DisplacedFace1, class Displacement1, int dimw1=Space1::GridView::Grid::dimension, class Space2, class Displacement2, int dimw2=Space2::GridView::Grid::dimension>
  auto getContactConstraint(Space1 const& space1, BoundaryLocator<typename Space1::GridView,Displacement1,Space1::GridView::Grid::dimension> const& boundaryLocator1,
                            Space2 const& space2, BoundaryLocator<typename Space2::GridView,Displacement2,Space2::GridView::Grid::dimension> const& boundaryLocator2,
                            DisplacedFace1 const& face,
                            /*typename BoundaryLocator<typename Space1::GridView,Displacement1,Space1::GridView::Grid::dimension>::Face const& face,*/
                            Dune::FieldVector<typename Space1::GridView::ctype,Space1::GridView::Grid::dimension-1> const& xi,
                            double overlap,
                            int which=0)
  {
    using ctype = typename Space1::GridView::ctype;
    using RowEntry = Dune::FieldMatrix<ctype,1,dimw1>;
    using EntryType = std::pair<size_t,RowEntry>;
    using ReturnType = std::optional<std::tuple<std::vector<EntryType>,std::vector<EntryType>,ctype>>;

    constexpr int dim1 = Space1::GridView::Grid::dimension;
    constexpr int dim2 = Space2::GridView::Grid::dimension;
    
    // face on boundary belonging to boundaryLocator1
    auto n = face.unitOuterNormal(xi);
    if (which)
      n = -n;

    double volume;
    if constexpr(dim1==dimw1) volume = face.gridFace().inside().geometry().volume();        //2D-2D or 3D-3D-contact: face is an intersection
    else                      volume = face.gridFace().geometry().volume();                 //2D grid in 3D world: face is a 2D cell
    ctype diameter = std::pow(volume,1.0/dim1);
    if(overlap == std::numeric_limits<double>::lowest()) 
      overlap = 0.05*diameter;
        
    // Find intersections of faces with the ray starting at the given point and extending in normal direction.
    // In a Newton loop, the current configuration need not be feasible (i.e. without overlap / interpenetration).
    // Hence we pull back a little bit (on the order of a single element diameter) to detect a slight overlap as well.
    auto intersection = boundaryLocator2.byRay(face.global(xi)-overlap*n,n);

    RowEntry nt; nt[0] = n;             // normal as a row vector

    if (!intersection)                  // no intersection found -> no constraint
      return ReturnType();

    // extract opposing boundary face and point
    auto oface = std::get<0>(intersection.value());
    auto oxi = std::get<1>(intersection.value());

    std::vector<EntryType> entries1;
    std::vector<EntryType> entries2;
    ctype const eps = std::numeric_limits<ctype>::epsilon();

    
    if constexpr(dim1==dimw1)
    {
      auto dofs  = space1.linearCombination(face.gridFace().inside(),face.gridFace().geometryInInside().global(xi));
      for (auto const& p: dofs)
        if (p.second.two_norm2() > 1e3*eps)                             // omit negligible contributions. Those are
          entries1.push_back(std::make_pair(p.first,p.second[0]*nt));   // typically structurally zero, but have rounding errors
    }
    else
    {
      auto dofs  = space1.linearCombination(face,xi/*face.geometry().global(xi)*/);
      for (auto const& p: dofs)
        if (p.second.two_norm2() > 1e3*eps)                            // omit negligible contributions. Those are
          entries1.push_back(std::make_pair(p.first,p.second[0]*nt));  // typically structurally zero, but have rounding errors
    }

    if constexpr(dim2==dimw2)
    {
      auto odofs = space2.linearCombination(oface.gridFace().inside(),oface.gridFace().geometryInInside().global(oxi));
      for (auto const& p: odofs)
        if (p.second.two_norm2() > 1e3*eps)
          entries2.push_back(std::make_pair(p.first,-p.second[0]*nt));
    }
    else
    {
      auto odofs = space2.linearCombination(oface,oxi/*oface.geometry().global(oxi)*/);
      for (auto const& p: odofs)
        if (p.second.two_norm2() > 1e3*eps)
          entries2.push_back(std::make_pair(p.first,-p.second[0]*nt));
    }

    auto x = face.global(xi);
    auto y = oface.global(oxi);

    // compute the gap distance (can be negative in case of interpenetration)
    return ReturnType(std::make_tuple(entries1,entries2,(y-x)*n));
  }
  
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------


  /**
   * \brief Abstract base class for computation of constraints from constraint samples.
   *
   * Let \f$ x \f$ denote the vector of displacements. On a boundary face, a set of \f$ n \f$
   * constraint samples is defined at local positions \f$ \xi_i \f$, leading to (pointwise)
   * nonpenetration conditions of the form
   * \f[ \sum_{k\in K_i} G_{ik} x_k \le g_i \f]
   * where \f$ g_i \f$ is the gap between opposing boundary points and \f$ G_i \f$ is the
   * impact of the displacements \f$ x_k \f$ on this gap, i.e. \f$ G_{ik} = \phi_k(\xi_i) \f$,
   * where \f$ \phi_k \f$ is the shape function associated with \f$ x_k \f$.
   * Note that the entries \f$ x_k \f$ may be (column) vector-valued, and the corresponding
   * matrix entries \f$ G_{ik} \f$ are (row) vector-valued.
   *
   * Now the actual constraints of the form \f$ Bx \le b \f$ can be defined by certain
   * linear combinations of the pointwise constraints. Taking them one by one as they are
   * results in pointwise constraints (with Dirac weights),
   * \f[ Gx \le g, \f] see \ref DiracMortar. If the pointwise constraints are multiplied
   * by test functions \f$ \psi_j \f$ and integrated (summed) over the
   * boundary face, different constraints result (see, e.g., \ref BezierMortar):
   * \f[ \forall j: \quad \sum_{i\in I} \psi_j(\xi_i) \sum_{k\in K_i} G_{ik} x_k
   *                      \le \sum_{i\in I} \psi_j(\xi_i) g_i. \f]
   * In any case, linear inequalities of the form \f$ Bx \le b \f$ are defined.
   *
   * This class defines an interface for computing \f$ B \f$ and \f$ b \f$ from \f$ G \f$
   * and \f$ g \f$.
   */
  template <class Scalar, int dim>
  class Mortar
  {
  public:

    /**
     * \brief The type of matrix entries in the constraint matrix \f$ B \f$ in \f$ Bx \le b \f$.
     *
     * The entries are row vectors with the same length as the primal variable \f$ x \f$ entries.
     */
    using Entry = Dune::FieldMatrix<Scalar,1,dim>;
    using Row = std::vector<std::pair<size_t,Entry>>;

    virtual ~Mortar();

    /**
     * \brief The number of resulting constraints if `n` samples are provided.
     */
    virtual int size(int n) const = 0;

    /**
     * \brief Updates the given set of constraints (by \f$ B \f$ and \f$ b \f$) when a
     *        new constraint sample becomes available.
     *
     * This method needs to be implemented in derived classes.
     *
     * \param n The total number of samples that may be provided (actually there may be fewer, because
     *          some points may lack a contact partner).
     * \param xi The face-local coordinate of the current contact sample.
     * \param Gi A sparse representation of the matrix row \f$ G_i \f$.
     * \param gi The gap \f$ g_i \f$
     * \param B The sparse matrix \f$ B \f$ to be updated. The update process can modify values
     *          as well as change the shape of the matrix.
     * \param b The bounds vector \f$ b \f$
     *
     * If B and b are empty on entry, a new face is started. The containers
     * are resized as needed.
     */
    virtual void updateConstraints(int n, Dune::FieldVector<Scalar,dim-1> const& xi,
                                   Row& Gi, Scalar gi,
                                   std::vector<Row>& B, std::vector<Scalar>& b) const = 0;
  };

  // ----------------------------------------------------------------------------------------------

  /**
   * \brief Defines a constraint formulation where each sample is taken as a single constraint of its own.
   *
   * For a detailed description of the concept, see the base class \ref Mortar. This implementation
   * creates constraints that correspond one to one to the provided constraint samples, i.e. the
   * mortar test functions are Dirac functions located at the constraint sample points.
   */
  template <class Scalar, int dim>
  class DiracMortar: public Mortar<Scalar,dim>
  {
  public:
    using typename Mortar<Scalar,dim>::Row;

    /**
     * \brief The number of resulting constraints if `n` samples are provided.
     */
    virtual int size(int n) const { return n; }

    /**
     * \brief The
     */
    virtual void updateConstraints(int n, Dune::FieldVector<Scalar,dim-1> const& xi,
                                   Row& Gi, Scalar gi,
                                   std::vector<Row>& B, std::vector<Scalar>& b) const;
  };

  /// \cond internals
  namespace ContactConstraintsDetail
  {
    template <class Scalar, int dim>
    constexpr DiracMortar<Scalar,dim> const& defaultMortar();
  }
  /// \endcond

  // ----------------------------------------------------------------------------------------------

  /**
   * \brief Defines a constraint formulation where samples are weighted by Bezier test functions.
   *
   * For a detailed description of the concept, see the base class \ref Mortar.
   * Let \f$ g(x_i) \f$ denote the gap evaluated at sample point \f$ x_i \f$, \f$ i=1,\dots,N \f$.
   * With \f$ M \f$ Bezier functions \f$ B_k \f$ (depends on the Bezier order given), we obtain
   * \f$ M \f$ constraints \f[ \sum_{i=1}^N g(x_i) B_k(x_i) \ge 0. \f]
   */
  template <class Scalar, int dim>
  class BezierMortar: public Mortar<Scalar,dim>
  {
  public:
    using typename Mortar<Scalar,dim>::Row;

    /**
     * \brief Creates a Bezier mortar formulation of order `m`.
     *
     * \param m The order of Bezier polynomials to use as constraint test functions (m >= 0).
     */
    BezierMortar(int m);

    /**
     * \brief The number of resulting constraints if `m` samples are provided.
     */
    virtual int size(int n) const;

    /**
     * \brief The
     */
    virtual void updateConstraints(int n, Dune::FieldVector<Scalar,dim-1> const& xi,
                                   Row& vic, Scalar g,
                                   std::vector<Row>& rows, std::vector<Scalar>& bounds) const;
  private:
    int m;
  };

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup contact
   * \brief Defines the sample-based contact constraint formulation.
   */
  struct MortarB
  {
    /**
     * \brief Defines the type of constraint formulation.
     */
    enum class Type { Dirac, Bezier };


    static constexpr MortarB dirac()
    {
      return MortarB{Type::Dirac,0};
    }

    static constexpr MortarB bezier(int order)
    {
      return MortarB{Type::Bezier,order};
    }

    Type type;
    int  mortarOrder;
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------


  /**
   * \ingroup contact
   * \brief Computes contact sample points on faces
   *
   * Equidistant contact sample points are computed.
   *
   * \param gt the type of face geometry, usually line or triangle (only those are supported right now)
   * \param nodesPerEdge the number of sample points per one-dimensional face subentity
   *
   * For lines (2D contact), \f$ n \f$ sample points are created, whereas for triangles (3D contact)
   * \f$ n(n+1)/2 \f$ points are created.
   *
   * \return a vector of sample points with associated (constant) weights
   */
  template <int dim>
  std::vector<Dune::QuadraturePoint<double,dim>> constraintPositions(Dune::GeometryType gt, int nodesPerEdge)
  {
    assert(nodesPerEdge>0);
    
    
    using QP = Dune::QuadraturePoint<double,dim>;
    std::vector<QP> pos;
    
    if (gt.isLine())
    {
      if (nodesPerEdge==1)
      {
        Dune::FieldVector<double,dim> p{0.5};                       // minimal number
        pos.push_back(QP(p,1.0));                                   // -> just center position
      }
      else
      {
        for (int i=0; i<nodesPerEdge; ++i)                          // otherwise equidistant grid on [0,1] including
        {                                                           // the end points
          Dune::FieldVector<double,dim> p{i/(nodesPerEdge-1.0)};    // equidistant positions
          double w = (i==0 || i+1==nodesPerEdge)? 0.5/(nodesPerEdge-1): 1.0/(nodesPerEdge-1); // and trapezoidal rule
          pos.push_back(QP(p,w));                                   // Newton-Cotes does not appear to be really beneficial
        }
      }
    }
    else if (gt.isTriangle())
    {
      if (nodesPerEdge==1)
      {   
        Dune::FieldVector<double,dim> p{1.0/3,1.0/3};               // minimal number of constraints per triangle
        pos.push_back(QP(p,1.0));                                   // -> just center position
      }
      else
      {
        double w = 2.0/(nodesPerEdge*(nodesPerEdge+1));
        for (int ix=0; ix<nodesPerEdge; ++ix)                       // otherwise equidistant cartesian grid of points
        {
          double x = ix/(nodesPerEdge-1.0);
          for (int iy=0; iy<nodesPerEdge-ix; ++iy)                  // with required number of points per edge
          {
            double y = iy/(nodesPerEdge-1.0);
            Dune::FieldVector<double,dim> p{x,y};
            pos.push_back(QP(p,w));                                 // and including the vertices (closed formula)
          }
        }
      }
    }
    else
    {
      std::cerr << "unsupported geometry type of face!\n";
      abort();
    }
    
    return pos;
  }
  
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup contact 
   * \brief Computes a complete set of pointwise multibody contact constraints.
   * 
   * This computes a linear inequality constraint of the form \f$ B \delta u \le b \f$. If satisfied by
   * some small displacement \f$ \delta u \f$, this guarantees (up to linearization and discretization)
   * mutual nonpenetration and self-nonpenetration.
   * 
   * \param boundaryLocator the spatial index of boundary faces
   * \param pointsPerEdge the number of points per edge of the boundary face
   * \param overlap states how far into the interior for an opposed face is searched.
   *                If not given, a default value (small positive value) will be chosen in dependence of the face size.
   * \param mortar the type of constraint agglomeration/averaging
   * \return a pair \f$ (B,b) \f$ of constraint matrix and right hand side
   *
   * \see getContactConstraint
   */
  template <class Space, class Displacement>
  auto getContactConstraints(Space const& space,
                             BoundaryLocator<typename Space::GridView,Displacement> const& boundaryLocator,
                             int pointsPerEdge, double overlap = std::numeric_limits<double>::lowest(),
                             Mortar<typename Space::Grid::ctype,Space::Grid::dimension> const& mortar = ContactConstraintsDetail::defaultMortar<typename Space::Grid::ctype,Space::Grid::dimension>())
  {
    ScopedTimingSection timer("getContactConstraints");
    using Grid = typename Space::Grid;
    using ctype = typename Grid::ctype;
    constexpr int dim = Grid::dimension;
    using Entry = Dune::FieldMatrix<ctype,1,dim>;
    using Row = std::vector<std::pair<size_t,Entry>>;
    
    std::vector<Row> rows;      // for storing
    std::vector<ctype> bounds;
    std::mutex mutex;           // for preventing concurrent access to these buffers
    
    auto faces = boundaryLocator.displacedFaces();

    // In each quadrature point on a boundary face, look in normal direction and find the 
    // corresponding algebraic constraint.
    parallelFor(0,faces.size(),[&](size_t i)
    {
      auto const& face = faces[i];
      auto const qr = constraintPositions<dim-1>(face.type(),pointsPerEdge);
      auto volume = face.volume();


      std::vector<Row> tmpRows;
      std::vector<ctype> tmpBounds;

      for (auto const& qp: qr)
      {
        auto const& xi = qp.position();
        auto c = getContactConstraint(space,boundaryLocator,face,xi,overlap); // find potential contact partner

        if (c)                                                                // there need not be any - then skip
        {
          ctype w = qp.weight()*volume;              // get the quadrature weight
          auto [vic,g] = *c;                         // get
          for (auto& p: vic)                         // and scale the constraint by this weight
            p.second *= w;
          g *= w;

          mortar.updateConstraints(qr.size(),xi,vic,g,tmpRows,tmpBounds);
        }
      }

      // Append the constraints of this face to the complete list of global constraints
      std::lock_guard lock(mutex);

      rows.reserve(size(rows)+size(tmpRows));
      for (auto& tr: tmpRows)
        rows.push_back(std::move(tr));

      bounds.reserve(size(bounds)+size(tmpBounds));
      for (auto& tb: tmpBounds)
        bounds.push_back(std::move(tb));
    });


    // Now convert the data into a CSR matrix
    ScopedTimingSection("convert to NumaBCRSMatrix");
    
    NumaCRSPatternCreator<> creator(rows.size(),space.degreesOfFreedom());
    for (size_t i=0; i<rows.size(); ++i)
      for (auto const& p: rows[i])
        creator.addElement(i,p.first);

      
    NumaBCRSMatrix<Entry> B(creator);
    for (size_t i=0; i<rows.size(); ++i)
      for (auto const& p: rows[i])
        B[i][p.first] += p.second;

      
    Dune::BlockVector<Dune::FieldVector<ctype,1>> b(bounds.size());
    std::copy(begin(bounds),end(bounds),begin(b));

    
    return std::make_pair(B,b);
  }

  /**
   * \ingroup contact
   * \brief Computes a complete pointwise twobody contact constraint.
   *
   * This version is for problems where each body is represented by its own grid.
   *
   * \param boundaryLocator1 the spatial index of boundary faces of first body/grid
   * \param boundaryLocator2 the spatial index of boundary faces of second body/grid
   * \param order the order of the quadrature rule that is used for selecting test points in the faces
   * \param overlap states how far into the interior for an opposed face is searched.
   * \param symmetric if true, this loops over both grids and collects pointwise constraints. If false, only boundary
   *        faces of first grid are considered
   * Negative value means search is only performed in the exterior with corresponding minimal distance.
   * If not given a default value (small positive value) will be chosen.
   * \param space1 if given contact constraints are computed with respect to this space, otherwise the space from boundaryLocator1 is used
   * \param space2 analogously to space1 but for the second grid
   * \return a three tuple \f$ (B1,B2,b) \f$ of constraint matrix blocks (belonging to grids) and right hand side
   */
  template <class Space1, class Displacement1, class Space2, class Displacement2, bool flattened=false>
  auto getContactConstraints(Space1 const& space1, BoundaryLocator<typename Space1::GridView,Displacement1,Space1::GridView::Grid::dimension> const& boundaryLocator1,
                             Space2 const& space2, BoundaryLocator<typename Space2::GridView,Displacement2,Space1::GridView::Grid::dimension> const& boundaryLocator2,
                             int order, double overlap = std::numeric_limits<double>::lowest(), bool symmetric=true)
  {
    constexpr int dim1 = Space1::GridView::Grid::dimension;
    constexpr int dim2 = Space2::GridView::Grid::dimension;
    constexpr int dimw1= Space1::GridView::Grid::dimensionworld;
    constexpr int dimw2= Space2::GridView::Grid::dimensionworld;
    constexpr int dimensionworld = std::max(dimw1,dimw2);
    constexpr int dimension = std::min(dim1, dim2);
    
    using ctype = typename Space1::GridView::ctype;
    constexpr int blockSize = flattened ? 1 : dimensionworld;
    constexpr int blocksPerNode = flattened ? dimensionworld : 1;
    using Entry = Dune::FieldMatrix<ctype,1,blockSize>;
    using Row = std::vector<std::pair<size_t,Dune::FieldMatrix<ctype,1,dimensionworld>>>;
    //  using QR = Dune::QuadratureRule<ctype,dimw1-1>;

    std::vector<Row> rows1;
    std::vector<Row> rows2;
    std::vector<ctype> bounds;

    auto faces = boundaryLocator1.displacedFaces();
    std::vector<typename BoundaryLocator<typename Space1::GridView,Displacement1,dimw1>::DisplacedFace> reducedFaces; //this is either an Intersection or a Cell

    // TODO: parallelFor loop
    for (auto const& face: faces)
    {
      auto volume = face.volume();
      auto const qr = constraintPositions<dimw1-1>(face.type(), order);

      for (auto const& qp: qr)
      {
        auto const& xi = qp.position();
        auto w = qp.weight()*volume;
        auto c = getContactConstraint(space1,boundaryLocator1,space2,boundaryLocator2,face,xi,overlap);
        if (c)
        {
          auto [vic1,vic2,g] = *c;                   // get constraint
          for (auto& p: vic1)                        // and scale it by this weight
            p.second *= w;
          for (auto& p: vic2)                        // and scale it by this weight
            p.second *= w;
          g *= w;
          rows1.push_back(vic1);
          rows2.push_back(vic2);
          bounds.push_back(g);
          reducedFaces.push_back(face);
        }
        
        if constexpr(dimension!=dimensionworld)
        {
          // check second normal direction (0,0,1) and (0,0,-1)
          // for 2D-3D contact or 2D-2D in  a 3D world
          auto cc = getContactConstraint(space1,boundaryLocator1,space2,boundaryLocator2,face,xi,overlap,1);
          if (cc)
          {
            auto [vic1,vic2,g] = *cc;                   // get constraint
            for (auto& p: vic1)                        // and scale it by this weight
              p.second *= w;
            for (auto& p: vic2)                        // and scale it by this weight
              p.second *= w;
            g *= w;
            rows1.push_back(vic1);
            rows2.push_back(vic2);
            bounds.push_back(g);
            reducedFaces.push_back(face);
          }
        }
      }
    }
    
    if(symmetric)
    {
      // loops over second geometry for overconstrained contact formulation
      auto faces2 = boundaryLocator2.displacedFaces();
      std::vector<typename BoundaryLocator<typename Space2::GridView,Displacement2,dimw2>::DisplacedFace> reducedFaces2; //this is either an Intersection or a Cell

      // TODO: parallelFor loop
      for (auto const& face: faces2)
      {
        auto volume = face.volume();
        auto const qr = constraintPositions<dimw1-1>(face.type(), order);

        for (auto const& qp: qr)
        {
          auto const& xi = qp.position();
          auto w = qp.weight()*volume;
          auto c = getContactConstraint(space2,boundaryLocator2,space1,boundaryLocator1,face,xi,overlap);
          if (c)
          {
            auto [vic1,vic2,g] = *c;                   // get constraint
            for (auto& p: vic1)                        // and scale it by this weight
              p.second *= w;
            for (auto& p: vic2)                        // and scale it by this weight
              p.second *= w;
            g *= w;
            rows2.push_back(vic1);
            rows1.push_back(vic2);
            bounds.push_back(g);
            reducedFaces2.push_back(face);
          }
          
          if constexpr(dimension!=dimensionworld)
          {
            // check second normal direction (0,0,1) and (0,0,-1)
            // for 2D-3D contact or 2D-2D in  a 3D world
            auto cc = getContactConstraint(space2,boundaryLocator2,space1,boundaryLocator1,face,xi,overlap,1);
            if (cc)
            {
              auto [vic1,vic2,g] = *cc;                   // get constraint
              for (auto& p: vic1)                        // and scale it by this weight
                p.second *= w;
              for (auto& p: vic2)                        // and scale it by this weight
                p.second *= w;
              g *= w;
              rows2.push_back(vic1);
              rows1.push_back(vic2);
              bounds.push_back(g);
              reducedFaces.push_back(face);
            }
          }
        }
      } 
    }

    NumaCRSPatternCreator<> creator1(rows1.size(),space1.degreesOfFreedom());
    for (size_t i=0; i<rows1.size(); ++i)
      for (auto const& p: rows1[i])
        for (int l=0; l<blocksPerNode; ++l)
          creator1.addElement(i,p.first*blocksPerNode+l);

    NumaBCRSMatrix<Entry> B1(creator1);
    for (size_t i=0; i<rows1.size(); ++i)
      for (auto const& p: rows1[i]) {
        if(!flattened) B1[i][p.first] = p.second;
        else
          for (int l=0; l<blocksPerNode; ++l)
            B1[i][p.first*blocksPerNode+l] = p.second[0][l];
      }

    NumaCRSPatternCreator<> creator2(rows2.size(),space2.degreesOfFreedom());
    for (size_t i=0; i<rows2.size(); ++i)
      for (auto const& p: rows2[i])
        for (int l=0; l<blocksPerNode; ++l)
          creator2.addElement(i,p.first*blocksPerNode+l);

    NumaBCRSMatrix<Entry> B2(creator2);
    for (size_t i=0; i<rows2.size(); ++i)
      for (auto const& p: rows2[i]) {
        if(!flattened) B2[i][p.first] = p.second;
        else
          for (int l=0; l<blocksPerNode; ++l)
            B2[i][p.first*blocksPerNode+l] = p.second[0][l];
      }

    Dune::BlockVector<Dune::FieldVector<ctype,1>> b(bounds.size());
    std::copy(begin(bounds),end(bounds),begin(b));

    return std::make_tuple(B1,B2,b,reducedFaces);
  }

  // ---------------------------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------------------------

  /**
   * \ingroup contact
   * \brief Computes the contact area for the first body.
   *
   *  The threshold value maxGap is not so easy to choose, maybe this is not the best way to compute the contact area.
   *
   * \param boundaryLocator1 the spatial index of boundary faces of first body/grid
   * \param boundaryLocator2 the spatial index of boundary faces of second body/grid
   * \param order the order of the quadrature rule that is used for selecting test points in the faces
   * \param overlap states how far into the interior for an opposed face is searched.
   *                Negative value means search is only performed in the exterior with corresponding minimal distance.
   *                If not given a default value (small positive value) will be chosen.
   * \param maxGap threshold value, distances below this value are treated as "in contact".
   * \return the area
   */
  template <class GridView, class Displacement>
  double getContactArea(BoundaryLocator<GridView,Displacement> const& boundaryLocator1,
                        BoundaryLocator<GridView,Displacement> const& boundaryLocator2,
                        int order, double overlap = std::numeric_limits<double>::lowest(),
                        double maxGap = 1e-4)
  {
    using ctype = typename GridView::ctype;
    using QR = Dune::QuadratureRule<ctype,GridView::Grid::dimension-1>;

    auto faces = boundaryLocator1.displacedFaces();
    QuadratureTraits<QR> quadratureRules;

    double area = 0.0;

    for (auto const& face: faces)
    {
      //  QR const& qr = quadratureRules.rule(face.geometry().type(),order);
      auto const qr = constraintPositions<GridView::Grid::dimension-1>(face.type(), order);


      // usage of standard quadrature rule for integration of discontinuous contact area indicator function is maybe not the best idea, but its convenient
      for (auto const& qp: qr)
      {
        auto const& xi = qp.position();
        auto c = getContactConstraint(boundaryLocator1.displacement()->space(),boundaryLocator1,boundaryLocator1.displacement()->space(),boundaryLocator2,face,xi,overlap);

        if (c && std::get<2>(*c) < maxGap)
        {
          auto dupi = boundaryLocator1.displacement()->derivative(face.gridFace().inside(), face.gridFace().geometryInInside().global(xi));
          dupi += unitMatrix<ctype, GridView::Grid::dimension>();

          auto jt = face.gridFace().geometry().jacobianTransposed(xi).rightmultiplyany(Dune::transpose(dupi));
          ctype intElem = std::sqrt(std::abs(jt.rightmultiplyany(Dune::transpose(jt)).determinant()));
          area += intElem * qp.weight();
        }
      }
    }

    return area;
  }

  // ---------------------------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------------------------

  /**
   * \ingroup contact 
   * \brief Computes a feasible step size for contact constraints.
   * 
   * The contact constraints as, e.g., obtained from \ref getContactConstraints, are linearized, in that they restrict 
   * the displacement of boundary points in \em normal direction. A solution of this linearized contact problem may 
   * be infeasible: other body parts may not have been detected as obstacles in normal direction, but now become relevant.
   * Taking a full step often leads to infeasible deformation states with considerable interpenetration. Reducing the 
   * stepsize appropriately, i.e. taking the largest possible step such that the configuration remains feasible not only 
   * in the linearized, but in the full nonlinear setting, is a promising way.
   * 
   * This function computes a feasible step size that is nearly optimal.
   * 
   * \tparam Displacement a finite element function type for representing volume displacements
   * \tparam BoundaryDisplacement a boundary displacement
   * 
   * \param u               the volume displacement, starting point for line search
   * \param boundaryLocator a spatial index for looking up ray-face intersections. Its internal boundary displacement
   *                        will be modified during the line search
   * \param direction       the step to be taken (linearized contact solution)
   * \param pointsPerEdge   the number of contact sample points per edge
   * \param trialStepSize   how aggressively the method goes ahead. 1=fast but with risk of ending up infeasible, 0=very slow but safe. 
   *                        In the range ]0,1[. Use smaller values for highly nonlinear problems with large displacements.
   * \param overlap         the amount of penetration to be accepted when looking backwards for contact
   */
  template <class Displacement, class BoundaryDisplacement>
  double contactLinesearch(Displacement const& u,
                           BoundaryLocator<typename Displacement::GridView,BoundaryDisplacement>& boundaryLocator,
                           Displacement const& direction,
                           int pointsPerEdge, double trialStepSize=0.2,
                           double overlap = std::numeric_limits<double>::lowest(),
                           double targetOverlap = 0,
                           Mortar<typename Displacement::ctype,Displacement::GridView::dimension> const& mortar =
                             ContactConstraintsDetail::defaultMortar<typename Displacement::ctype,Displacement::GridView::dimension>())
  {
    using GridView = typename Displacement::GridView;
    static_assert(GridView::dimensionworld==Displacement::components);
    assert(pointsPerEdge>=0);
    assert(targetOverlap<std::max(0.0,overlap) || (targetOverlap==0 && overlap==std::numeric_limits<double>::lowest()));
    
    double tmin = 0;        // current position on the full step [0,1].
    double tmax = 1;        // linear model predicts we can take a full step.

    
    assert(boundaryLocator.displacement());


    // Function for computing the minimal remaining gap depending on the
    // step size t we take in direction of the given displacement "direction".
    auto getGap = [&](double t)
    {
      // First get the constraints at the new trial point.
      auto uTrial = u;                                    // uTrial <- u + t*direction
      uTrial.axpy(t,direction);
      boundaryLocator.updateDisplacement(uTrial);         // build spatial index for new deformation
      auto Bb = getContactConstraints(u.space(),boundaryLocator,pointsPerEdge,overlap,mortar);  // get constraints B*dx <= b
      auto Bx(Bb.second); Bb.first.mv(direction.coefficients(),Bx);
      auto const& b = Bb.second;                          // the gap

      // Now step through all constraints and find, for each of them,
      // the maximal feasible step size.
      double gmin = std::numeric_limits<double>::max();
      double dg, taumax = gmin;
      for (int i=0; i<b.N(); ++i)
      {
        // Find the minimal gap
        if (b[i][0] < gmin)
        {
          gmin = b[i][0];
          dg = -Bx[i][0];
        }

        // Find the maximal allowed step size, i.e. max t s.t. (B*direction)[i]*t <= b[i]+targetOverlap.
        // This limits t only if (B*direction)[i] > 0 (otherwise the gap is increased, and we obtain a
        // (hopefully negative) minimal step size.
        if (Bx[i][0] > 0)
          taumax = std::min(taumax,(b[i][0]+targetOverlap)/Bx[i][0]);
      }

      return std::make_tuple(gmin,dg,taumax);
    };

    // step ahead until we find an infeasible step size.
    // This cannot be done in one step up to tmax, because we may overlook constraints.
    double const h = trialStepSize;
    std::cout << "starting linesearch in [" << tmin+h << "," << tmax << "]\n";
    while (tmin < 0.999*tmax)
    {
      double t = std::min(tmax,tmin+h);
      auto [g0,dg0,tau0] = getGap(t);
      //  std::cout << "gap(t=" << t << ") = " << g0 << " [dg = " << dg0 << ", taumax = " << tau0 << "]\n";
      if (g0 < -targetOverlap)
        break;
      else
        tmin = t;
    }

    if (tmin >= 0.999*tmax)
    {
      std::cout << "Found stepsize t=" << tmin << ".\n";
      return tmin;
    } else
    {
      std::cout << "No feasible stepsize found in this interval.\n";
      tmax = std::min(tmax,tmin+h);
    }

    std::cout << "starting linesearch in [" << tmin << "," << tmax << "]\n";
    while (tmax-tmin > std::min(1e-2,tmax/10))
    {
      double t = (tmax+tmin)/2;
      if(t<1e-12)
        break;
      auto [g0,dg0,tau0] = getGap(t);
      //  std::cout << "gap(t=" << t << ") = " << g0 << " [dg = " << dg0 << ", taumax = " << tau0 << "]\n";
     if (g0 < -targetOverlap)
        tmax = t;
      else
        tmin = t;
    }

    std::cout << "Found stepsize t=" << tmin << ".\n";
    return tmin;
  }
  
  /**
   * \ingroup contact
   * \brief Computes a feasible step size for contact constraints.
   *
   * The contact constraints as, e.g., obtained from \ref getContactConstraints, are linearized, in that they restrict
   * the displacement of boundary points in \em normal direction. A solution of this linearized contact problem may
   * be infeasible: other body parts may not have been detected as obstacles in normal direction, but now become relevant.
   * Taking a full step often leads to infeasible deformation states with considerable interpenetration. Reducing the
   * stepsize appropriately, i.e. taking the largest possible step such that the configuration remains feasible not only
   * in the linearized, but in the full nonlinear setting, is a promising way.
   *
   * This function computes a feasible step size that is nearly optimal.
   * 
   * This variant is using two BoundaryLocators to determine contact between objects
   * discretized with two individual grids.
   *
   * \tparam Displacement1            a FunctionSpaceElement  (first object)
   * \tparam BoundaryDisplacement1    a boundary displacement (first object)
   * \tparam Displacement2            a FunctionSpaceElement  (second object)
   * \tparam BoundaryDisplacement2    a boundary displacement (second object)
   *
   * \param u1                        the current displacement of the first object, starting point for the step
   * \param boundaryLocator1          the spatial index for looking up ray-face intersections, fitst object
   * \param direction1                the step to be taken by object 1 (linearized contact solution)
   * \param u2                        the current displacement of the second object
   * \param boundaryLocator2          the spatial index for looking up ray-face intersections, second object.
   * \param direction2                the step to be taken by object 2 (linearized contact solution)
   * \param order                     the order of the quadrature rule that is used for selecting test points in the faces
   * \param trialStepSize             how aggressively the method goes ahead. 1=fast but with risk of ending up infeasible, 0=very slow but safe.
   *                                  In the range ]0,1[. Use smaller values for highly nonlinear problems with large displacements.
   * \param overlap                   the amount of infeasibility that is to be accepted.
   */
  template <class Displacement1, class BoundaryDisplacement1, class Displacement2, class BoundaryDisplacement2>
  double contactLinesearch(Displacement1 const& u1, BoundaryLocator<typename Displacement1::GridView,BoundaryDisplacement1>& boundaryLocator1,
                           Displacement1 const& direction1,
                           Displacement2 const& u2, BoundaryLocator<typename Displacement2::GridView,BoundaryDisplacement2>& boundaryLocator2,
                           Displacement2 const& direction2,
                           int order, double trialStepSize=0.2,
                           double overlap = std::numeric_limits<double>::lowest(),
                           double targetOverlap = 0)
  {
    assert(order>=0);
    assert(0<trialStepSize && trialStepSize<1);
    assert(targetOverlap<overlap || (targetOverlap==0 && overlap==std::numeric_limits<double>::lowest())); 
    assert(Displacement1::GridView::Grid::dimensionworld==Displacement2::GridView::Grid::dimensionworld);
    assert(boundaryLocator1.displacement());
    assert(boundaryLocator2.displacement());

    double tmin = 0;        // current position on the full step [0,1].
    double tmax = 1;        // linear model predicts we can take a full step.

    size_t nCoeff1 = direction1.coefficients().N();
    size_t nCoeff2 = direction2.coefficients().N();

    // For a general documentation of this lambda, compare the single-grid
    // contactLinesearch function above.
    auto getGap = [&](double t)
    {
      auto uTrial1 = u1;
      uTrial1.axpy(t,direction1);
      auto uTrial2 = u2;
      uTrial2.axpy(t,direction2);
      boundaryLocator1.updateDisplacement(uTrial1);                              // build spatial index for new deformation
      boundaryLocator2.updateDisplacement(uTrial2);                            
      auto Bb = getContactConstraints(u1.space(),boundaryLocator1,u2.space(),boundaryLocator2,order,overlap); // get constraints B1*du1 + B2*du2 <= b
      
      auto const& B = horzcat(std::get<0>(Bb),std::get<1>(Bb));
      auto const& g = std::get<2>(Bb);
      typename Displacement1::StorageType combinedDirection(nCoeff1+nCoeff2);
      for(size_t i=0; i< nCoeff1; ++i) 
          combinedDirection[i]=direction1.coefficients()[i];
      for(size_t i=0; i< nCoeff2; ++i) 
          combinedDirection[i+nCoeff1]=direction2.coefficients()[i];
      
      auto Bx(g); B.mv(combinedDirection,Bx); Bx *= -1;

      // Now step through all constraints and find, for each of them,
      // the maximal feasible step size.
      double gmin = std::numeric_limits<double>::max();
      double dg, taumax = gmin;
      for (int i=0; i<g.N(); ++i)
      {
        // find the minimal gap
        if (g[i][0] < gmin)
        {
          gmin = g[i][0];
          dg = Bx[i][0];
        }
        if (Bx[i][0] < 0)
          taumax = std::min(taumax,-(g[i][0]+targetOverlap)/Bx[i][0]);
      }
      
      //  std::cout << "tau max: " << taumax << std::endl;

      return std::make_tuple(gmin,dg,taumax);
    };
    
    // step ahead until we find an infeasible step size.
    // This cannot be done in one step up to tmax, because e may overlook constraints.
    double const h = trialStepSize;
    while (tmin < 0.999*tmax)
    {
      double t = std::min(tmax,tmin+h);
      auto [g0,dg0,tau0] = getGap(t);
std::cout << "gap(t=" << t << ") = " << g0 << " [dg = " << dg0 << ", taumax = " << tau0 << "]\n";
      if (g0 < -targetOverlap)
        break;
      else
        tmin = t;
    }

    if (tmin >= 0.999*tmax)
      return tmin;
    else
      tmax = std::min(tmax,tmin+h);

std::cout << "starting linesearch in [" << tmin << "," << tmax << "]\n";
    while (tmax-tmin > std::min(1e-4,tmax/10))
    {
      double t = (tmax+tmin)/2;
      auto [g0,dg0,tau0] = getGap(t);
//  std::cout << "gap(t=" << t << ") = " << g0 << " [dg = " << dg0 << ", taumax = " << tau0 << "]\n";
      if (g0 < -targetOverlap)
        tmax = t;
      else
        tmin = t;
    }

    return tmin;
  } 
}

#endif
