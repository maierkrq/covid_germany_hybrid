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

#ifndef FEM_LAGRANGESHAPEFUNCTIONS_HH
#define FEM_LAGRANGESHAPEFUNCTIONS_HH

/**
 * @file
 * @brief  define Lagrange type shape functions for simplicial elements of
 *         arbitrary dimension and order
 * @author Martin Weiser
 */


#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <tuple>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include "fem/barycentric.hh"
#include "fem/pshapefunctions.hh"
#include "linalg/dynamicMatrix.hh"
#include "utilities/combinatorics.hh"
#include "utilities/deprecation.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/power.hh"
#include "utilities/threading.hh"

//---------------------------------------------------------------------

namespace Kaskade
{
  // forward declarations
  template <class ctype, int dimension, class Scalar> 
  class LagrangeSimplexShapeFunctionSet;

  /**
   * \cond internals
   */
  namespace SimplexLagrangeDetail {
    

    /**
     * \brief Computes the number of Lagrange shape functions of given order on
     * the unit simplex of dimension dim.  For order<0, the number is 0.
     */
    inline int size(int dim, int order)
    {
      if (order < 0)
        return 0;

      return binomial(dim+order,dim);
    }

    /**
     * \brief Computes sequential index of the shape function / interpolation node in the
     * set of shape functions / interpolation nodes from the given tuple index.
     * \param xi an array of size dim containing the tuple index
     * \param dim the spatial dimension of the simplex
     * \param order the polynomial ansatz order
     */
    inline int local(int const* xi, int dim, int order)
    {
      assert(std::accumulate(xi,xi+dim,0)<=order);

      int local = 0;
      for (int dir=0; dir<dim; ++dir)
        for (int i=0; i<xi[dir]; ++i) {
          local += size(dim-1-dir,order);
          --order;
        }
      return local;
    }

    /**
     * \brief Compute the tuple index of the interpolation point with given sequential index.
     * The tuple index is a dim-tuple of integers with values between 0 and order (including),
     * such that the total sum is at most order.
     * \tparam dim spatial dimension of the simplex
     * \param order polynomial ansatz order
     * \param loc number of the Lagrangian grid point (0<=loc<size(dim,order))
     * \returns an array of size d that contains the tuple index
     */
    template <int dim>
    std::array<int,dim> tupleIndex(int const order, int const loc)
    {
      int m = order;
      int k = loc;
      std::array<int,dim> xi;
      std::fill(begin(xi),end(xi),0);

      for (int dir=0; dir<dim; ++dir) {
        int s = size(dim-1-dir,m);
        while (k>=s && s>0) {
          k -= s;
          --m;
          ++xi[dir];
          s = size(dim-1-dir,m);
        }
      }
      assert(local(&xi[0],dim,order)==loc);
      return xi;
    }

    /**
     * \brief Computes the nodal position in the unit simplex from the tuple index.
     */
    template <int dim>
    Dune::FieldVector<double,dim> nodalPosition(std::array<int,dim> const& xi, int const order)
    {
      // construct equidistant interpolation nodes. max  in denominator is just to prevent GCC from
      // choking on division by 0.
      Dune::FieldVector<double,dim> x;
      for (int i=0; i<dim; ++i)
        x[i] =  order==0? 1.0/(dim+1): static_cast<double>(xi[i])/(order>0? order: 1);
      return x;
    }


    /**
     * \brief Computes the largest codimension of the subentities containing the given barycentric coordinate.
     */
    template <size_t dim_1>
    int codim(std::array<int,dim_1> const& xi)
    {
      // That's simple. If a barycentric coordinate is 0, the node is
      // located on the associated facet (codim=1) of the simplex. The
      // node lies on the intersection of all facets on which is
      // located. Each incident facet increases the codimension by one.
      int count = 0;
      int dim = dim_1 - 1;

      for (int i=0; i<=dim; ++i)
        if (xi[i]==0) ++count;

      if (count > dim) // all entries 0 -> order==0
        return 0;

      return count;
    }

    // Given a barycentric coordinate vector, this function computes the
    // simplex local id of the largest codim subentity on which the node
    // is located. The implementation is based on the reference element
    // documentation of Dune.
    template <size_t dim_1>
    int entity(std::array<int,dim_1> const& xi)
    {
      int dim = dim_1 - 1;

      switch (dim) {
      case 1: // lines
        switch (codim(xi)) {
        case 0: // edge
          return 0;
        case 1: // vertex
          if (xi[0]!=0) return 1;
          if (xi[1]!=0) return 0;
          break;
        default: assert(false);
        }
        break;
        case 2: // triangles
          switch (codim(xi)) {
          case 0: // triangle
            return 0;
          case 1: // edge
            if (xi[0]==0) return 1;
            if (xi[1]==0) return /*2*/0;
            if (xi[2]==0) return /*0*/2; //new Dune 2.0 numbering: 2-oldIndex
            break;
          case 2: // vertex
            if (xi[0]!=0) return 1;
            if (xi[1]!=0) return 2;
            if (xi[2]!=0) return 0;
            break;
          }
          break;
          case 3: // tetrahedra
            switch (codim(xi)) {
            case 0: // tetrahedron
              return 0;
            case 1: // triangle
              if (xi[0]==0) return /*1*/ 2 ; //new Dune 2.0 numbering: 3-oldIndex
              if (xi[1]==0) return /*2*/ 1 ;
              if (xi[2]==0) return /*3*/ 0 ;
              if (xi[3]==0) return /*0*/ 3 ;
              break;
            case 2: // edge
              if (xi[0]==0 && xi[1]==0) return 3;
              if (xi[0]==0 && xi[2]==0) return /*2*/1; //new Dune 2.0 numbering: xchg(1,2)
              if (xi[0]==0 && xi[3]==0) return 5;
              if (xi[1]==0 && xi[2]==0) return 0;
              if (xi[1]==0 && xi[3]==0) return 4;
              if (xi[2]==0 && xi[3]==0) return /*1*/2;
              break;
            case 3: // vertex
              if (xi[0]!=0) return 1;
              if (xi[1]!=0) return 2;
              if (xi[2]!=0) return 3;
              if (xi[3]!=0) return 0;
              break;
            }
            break;
      }
      assert("Unknown dimension/codimension\n"==0);
      return -1;
    }

//     // cycle through all possible tuple indices of sum up to m.
//     template <int dim>
//     void increment(Dune::FieldVector<int,dim>& x, int const m)
//     {
//       for (int i=0; i<dim; ++i) {
//         int sum = std::accumulate(x.begin(),x.end(),0);
//         if (sum==m)  // carry
//           x[i] = 0;
//         else {
//           ++x[i];
//           break;
//         }
//       }
//     }



    template <class ctype, int dimension, class Scalar, bool restricted>
    struct ChooseShapeFunctionSet
    {
      typedef LagrangeSimplexShapeFunctionSet<ctype,dimension,Scalar> type;
    };

    template <class ctype, int dimension, class Scalar>
    struct ChooseShapeFunctionSet<ctype,dimension,Scalar,true>
    {
      typedef RestrictedShapeFunctionSet<LagrangeSimplexShapeFunctionSet<ctype,dimension,Scalar> > type;
    };
  } // End of namespace SimplexLagrangeDetail
  /**
   * \endcond
   */

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  
  template <int dim>
  class LagrangeBase
  {
  public:
    /**
     * \brief Constructor.
     * \param l the barycentric index of the shape function with nonnegative entries. 
     * The order of the shape function is \f$ p = \sum_{i=0}^d l_i \f$.
     */
    LagrangeBase(int order);
    
    /**
     * \brief The number of Lagrange polynomials.
     */
    int size() const 
    {
      return ls.size();
    }
    
    /**
     * \brief The polynomial ansatz order of this Lagrange basis.
     */
    int order() const 
    {
      return myOrder;
    }
    
    std::array<int,dim+1> index(int i) const 
    {
      return ls[i];
    }
    
  protected:
    int myOrder;
    std::vector<std::array<int,dim+1>> ls;    
  };
  
  //---------------------------------------------------------------------
  
  /**
   * \ingroup fem 
   * \brief Lagrange shape functions on an equidistant grid over the simplex.
   * These are polynomial basis functions \f$ p_i \f$ over the unit simplex, together
   * with interpolation nodes \f$ \xi_j \f$, such that \f$ p_i(\xi_j) = \delta_{ij} \f$
   * holds.
   */
  template <int dim, class Real=double>
  class EquidistantLagrange: public LagrangeBase<dim>
  {
  public:
    using Vector = Dune::FieldVector<Real,dim>;
    
    /**
     * \brief Constructor.
     * \param l the barycentric index of the shape function with nonnegative entries. 
     * The order of the shape function is \f$ p = \sum_{i=0}^d l_i \f$.
     */
    EquidistantLagrange(int order);
    
    /**
     * \brief The spatial location of the interpolation node \f$ \xi_i \f$ associated with
     *        the i-th basis function. 
     */
    Vector nodalPosition(int i) const;
    
    Real value(int i, Vector const& xi) const;
    Real derivative(int i, int dir, Vector const& xi) const;
    Real derivative2(int i, int dir1, int dir2, Vector const& xi) const;
    
  private:
    std::vector<Real> normalization;
  };
  
  //---------------------------------------------------------------------
  
  template <int dim, class Real=double>
  class LobattoGridLagrange: public LagrangeBase<dim>
  {
  public:
    using Vector = Dune::FieldVector<Real,dim>;
    
    LobattoGridLagrange(int order);
    Vector nodalPosition(int i) const;
    Real value(int i, Vector const& x) const;
    Real derivative(int i, int dir, Vector const& xi) const;
    Real derivative2(int i, int dir1, int dir2, Vector const& xi) const;
  private:
    std::vector<Real>             lobatto;  // Lobatto nodes on [0,1]
    EquidistantLagrange<dim,Real> basis;    // basis functions used for evaluation
    DynamicMatrix<Real>           vInv;     // inverse of the Vandermonde matrix
  };
  
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  
  /**
   * \brief Scalar Lagrange shape functions on the unit simplex. 
   * 
   * These are polynomial shape functions of fixed given order associated to
   * certain interpolation nodes. A shape function is 1 at its
   * associated node and 0 at all the other nodes.
   *
   * The nodes form a cartesian grid with (order+1) nodes on each
   * axis. The node placement is isotropic (the same node placement on
   * each axis). Each node is referenced by its tuple index: a
   * dim-tuple of integers with values between 0 and order (including),
   * such that the total sum is at most order.
   *
   *
   * \tparam ctype the scalar type for coordinates
   * \tparam dimension the spatial dimension
   * \tparam Scalar the value type (a real number type)
   */
  template <class ctype, int dimension, class Basis, class Scalar=double>
  class LagrangeSimplexShapeFunction: public ShapeFunction<ctype,dimension,Scalar>
  {
  public:
    static int const dim = dimension;
    static int const comps = 1;

    typedef ctype CoordType;
    typedef Scalar ResultType;

    /**
     * \brief Default constructor. 
     * 
     * Constructs the shape function located at the
     * origin of the reference element. Not particularly useful, but
     * it's often convenient to have a default constructible class.
     */
    LagrangeSimplexShapeFunction()
    : local_(-1), codim_(-1), entity_(-1), entityIndex_(-1) {}


    /**
     * \brief Creates a shape function associated to the node with tuple index xi_.
     * The entries in xi_ have to be nonnegative, and their sum must not
     * exceed order.
     */
    explicit LagrangeSimplexShapeFunction(int local, 
                                          std::shared_ptr<Basis const> basis_) 
    : local_(local)
    , basis(basis_)
    {
      int const order = basis->order();
      
      // Compute local subentity indices:
      // Interprete the barycentric index interior to the subentity as
      // complete barycentric index of a lower dimensional simplex
      // with reduced order p=order-(dim-codim_+1).
      if (order==0) 
      {
        codim_ = 0;
        entity_ = 0;
        entityIndex_ = 0;
      } 
      else 
      {
        // construct barycentric indices
        auto xi = basis->index(local_);

        codim_ = SimplexLagrangeDetail::codim(xi);
        entity_ = SimplexLagrangeDetail::entity(xi);

        int yi[dim-codim_+1];
        int j=0;
        for (int i=0; i<=dim; ++i)
          if (xi[i]!=0)
            yi[j++] = xi[i]-1;
        assert(j==dim-codim_+1);
        entityIndex_ = SimplexLagrangeDetail::local(yi,dim-codim_,order-(dim-codim_+1));
        assert(entityIndex_>=0);
      }
    }

    /**
     * \brief Copy constructor.
     */
    LagrangeSimplexShapeFunction(LagrangeSimplexShapeFunction const& other) = default;

    /**
     * \brief Evaluates the shape function at point x. 
     * 
     * The point x is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Dune::FieldVector<Scalar,1> evaluateFunction(Dune::FieldVector<CoordType,dim> const& x) const
    {
      return basis->value(local_,x);
    }

    /**
     * \brief Evaluates the derivative of the shape function. 
     * 
     * The result is r[i][j] = \f$ \partial \phi_i / \partial x_j \f$
     * 
     * The point x is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Dune::FieldMatrix<Scalar,1,dim> evaluateDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      Dune::FieldMatrix<Scalar,1,dim> d;
      for (int dir=0; dir<dim; ++dir)
        d[0][dir] = basis->derivative(local_,dir,x);
      return d;
    }

    /**
     * \brief Evaluates the second derivative of the shape function. 
     * 
     * The result is r[i][j][k] = \f$ \partial^2 \phi_i / (\partial x_j \partial x_k) \f$
     * 
     * The point x is not restricted to be inside the unit simplex, but the meaning of
     * evaluating a shape function outside of the unit simplex is
     * questionable.
     */
    virtual Tensor<Scalar,1,dim,dim> evaluate2ndDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      // exploit symmetry of 2nd derivative
      Tensor<Scalar,1,dim,dim> d;
      for (int dir1=0; dir1<dim; ++dir1)
        for (int dir2=0; dir2<=dir1; ++dir2)
          d[0][dir1][dir2] = basis->derivative2(local_,dir1,dir2,x);

      for (int dir1=0; dir1<dim; ++dir1)
        for (int dir2=dir1+1; dir2<dim; ++dir2)
          d[0][dir1][dir2] = d[0][dir2][dir1];

      return d;
    }

    virtual std::tuple<int,int,int,int> location() const
    {
      return std::make_tuple(basis->order(),codim_,entity_,entityIndex_);
    }

  private:
    int local_;
    int codim_;
    int entity_;
    int entityIndex_;

    std::shared_ptr<Basis const> basis;
  };



  //---------------------------------------------------------------------

  /**
   * \brief A container of Lagrange shape functions of order Ord on the unit
   * simplex of grid dimension. For Ord<0, this set is empty.
   */
  template <class ctype, int dimension, class Scalar>
  class LagrangeSimplexShapeFunctionSet: public ShapeFunctionSet<ctype,dimension,Scalar>
  {
//     using Basis = EquidistantLagrange<dimension,ctype>;
    using Basis = LobattoGridLagrange<dimension,ctype>;
    
  public:
    typedef LagrangeSimplexShapeFunction<ctype,dimension,Basis,Scalar> value_type;
    typedef typename ShapeFunctionSet<ctype,dimension,Scalar>::Matrix Matrix;

    explicit LagrangeSimplexShapeFunctionSet(int order)
    : ShapeFunctionSet<ctype,dimension,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex,dimension))
    , basis(new Basis(order<0? 0: order))
    {
      int const dim = dimension;
      int const n = basis->size();
      this->iNodes.resize(n);

      for (int i=0; i<n; ++i) 
      {
        sf.push_back(std::unique_ptr<ShapeFunction<ctype,dim,Scalar>>(new value_type(i,basis)));
        for (int j=0; j<dim; ++j)
          this->iNodes[i] = basis->nodalPosition(i);
      }
      this->order_ = order;
      this->size_ = n;
    }

    LagrangeSimplexShapeFunctionSet(LagrangeSimplexShapeFunctionSet const& other)
    : ShapeFunctionSet<ctype,dimension,Scalar>(Dune::GeometryType(Dune::GeometryType::simplex,dimension))
    {
      this->iNodes = other.iNodes;
      this->order_ = other.order_;
      this->size_ = other.size_;
      for(size_t i=0; i<other.sf.size(); ++i) 
        sf.push_back(std::make_unique<value_type>(dynamic_cast<value_type const&>(other[i])));
    }

    virtual ShapeFunction<ctype,dimension,Scalar> const& operator[](int i) const
    { 
      return *sf[i]; 
    }
    
    virtual void interpolate(typename ShapeFunctionSet<ctype,dimension,Scalar>::SfValueArray const& A,
                             Matrix& IA) const
    {
      assert(A.N()==sf.size());
      IA = A;
      // For Lagrange shape functions, the interpolation at their nodes is just the identity. Copy!
    }

    void removeShapeFunction(size_t index)
    {
      assert(index>=0 && index < sf.size());
      sf.erase(sf.begin()+index);
      this->iNodes.erase(this->iNodes.begin()+index);
    }

  protected:
    std::vector<std::unique_ptr<ShapeFunction<ctype,dimension,Scalar>>> sf;
    
  private:
    std::shared_ptr<Basis const> basis;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * \brief Scalar Lagrange shape functions on the unit cube. 
   * 
   * These are polynomial shape functions of fixed given order associated to
   * certain interpolation nodes. A shape function is 1 at its
   * associated node and 0 at all the other nodes.
   *
   * The nodes form a cartesian grid with (order+1) nodes on each
   * axis. The node placement is isotropic (the same node placement on
   * each axis). In the current implementation, quasi Chebyshev nodes
   * are used (due to their better interpolation properties compared to
   * equidistant nodes). Each node is referenced by its index: a
   * dim-tuple of integers with values between 0 and order (including),
   * such that the maximum entry is at most order.
   *
   * \tparam ctype the scalar type for coordinates
   * \tparam dim the spatial dimension
   * \tparam Scalar the value type (a real number type)
   * \tparam O the polynomial order
   */
  template <class ctype, int dimension, class Scalar, int O>
  class LagrangeCubeShapeFunction: public ShapeFunction<ctype,dimension,Scalar>
  {
  public:
    static int const dim = dimension;
    static int const comps = 1;
    static int const order = O;

    typedef ctype CoordType;
    typedef Scalar ResultType;

    /**
     * Default constructor. Constructs the shape function located at the
     * origin of the reference element. Not particularly useful, but
     * it's often convenient to have a default constructible class.
     */
    LagrangeCubeShapeFunction() :
      local_(-1)
    {}

    LagrangeCubeShapeFunction(LagrangeCubeShapeFunction const& other)
    : local_(other.local_), codim_(other.codim_), entity_(other.entity_), entityIndex_(other.entityIndex_),
       xi(other.xi), normalization(other.normalization), node(other.node)
    {}

    /**
     * \brief Creates a shape function associated to the node with indices xi_.
     * 
     * The entries in xi_ have to be nonnegative, and their maximum must not
     * exceed order.
     */
    explicit LagrangeCubeShapeFunction(Dune::FieldVector<int,dimension> const& xi_);

    /**
     * \brief Evaluates the shape function at point x. 
     * 
     * The point x is not restricted to be inside the unit cube, but the meaning of
     * evaluating a shape function outside of the unit cube is
     * questionable.
     */
    ResultType evaluateFunction(int /* comp */, Dune::FieldVector<CoordType,dim> const& x) const
    {
      return normalization * evaluateNonNormalized(x);
    }

    virtual Dune::FieldVector<Scalar,1> evaluateFunction(Dune::FieldVector<ctype,dim> const& x) const
    {
      return normalization * evaluateNonNormalized(x);
    }

    /**
     * \brief Evaluates the partial derivative of the shape function in coordinate direction dir.  
     * 
     * 0 <= dir < dim has to hold. The point x is
     * not restricted to be inside the unit cube, but the meaning of
     * evaluating a shape function outside of the unit cube is
     * questionable.
     */
    ResultType evaluateDerivative(int /* comp */, int dir, Dune::FieldVector<CoordType,dim> const& x) const
    {
      assert(0<=dir && dir<dim);
      return normalization * evaluateDerivativeNonNormalized(dir,x);
    }

    virtual Dune::FieldMatrix<ResultType,1,dim> evaluateDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      Dune::FieldMatrix<ResultType,1,dim> d;
      for (int dir=0; dir<dim; ++dir)
        d[0][dir] = normalization * evaluateDerivativeNonNormalized(dir,x);
      return d;
    }

 
    virtual Tensor<ResultType,1,dim,dim> evaluate2ndDerivative(Dune::FieldVector<CoordType,dim> const& x) const
    {
      Tensor<ResultType,1,dim,dim> dd;
      for (int dir1=0; dir1<dim; ++dir1)
        for (int dir2=0; dir2<dim; ++dir2)
          dd[0][dir1][dir2] = normalization * evaluate2ndDerivativeNonNormalized(dir1,dir2,x);
      return dd;
    }

    /**
     * \brief Returns the local index of the shape function.
     */
    int localindex(int /* comp */) const { return local_; }

    /**
     * \brief Returns the codimension of the subentity on which the associated node is located.
     */
    int codim() const { return codim_; }

    /**
     * \brief Returns the subentity number of the smallest codimension
     * subentity on which the associated node is located.
     */
    int entity() const { return entity_; }

    /**
     * \brief Returns the local index on the subentity. 
     * 
     * For higher order shape functions, globally unique numberings on intersections have to be
     * used, therefore this depends on the actual element and an index
     * set that covers this element.
     */
    template <class Cell>
    int entityIndex(Cell const& /* cell */) const 
    {
      // If ordering of dofs is completely local (cell interiors) or
      // globally unique (only one dof on subentity, e.g. vertex or order<=2), return the
      // precomputed value.
      if (entityIndex_ >= 0)
        return entityIndex_;

      // Otherwise we have to compute a globally unique ordering (which
      // will necessarily depend on the actual element).
      // @todo: implement this!
      assert("Not implemented!"==0);
      return -1;
    }

    /**
     * \brief Returns the element-local position of the node associated to the shape function.
     */
    Dune::FieldVector<CoordType,dim> position () const
      {
      Dune::FieldVector<CoordType,dim> pos;
      for (int i=0; i<dim; ++i)
        pos[i] = node[xi[i]];
      return pos;
      }

    virtual std::tuple<int,int,int,int> location() const {
      return std::make_tuple(order,codim_,entity_,entityIndex_);
    }

  private:
    int local_;
    int codim_;
    int entity_;
    int entityIndex_;
    Dune::FieldVector<int,dim> xi;
    ResultType normalization;
    Dune::FieldVector<CoordType,order+1> node;

    ResultType evaluateNonNormalized(Dune::FieldVector<CoordType,dim> const& x) const
    {
      ResultType phi = 1;
      for (int i=0; i<dim; ++i)
        for (int k=0; k<=order; ++k)
          if (k!=xi[i])
            phi *= x[i]-node[k];
      return phi;
    }

    Scalar evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<ctype,dimension> const& x) const;
    
    ResultType evaluate2ndDerivativeNonNormalized(int dir1, int dir2, Dune::FieldVector<CoordType,dim> const& x) const 
    {
      static bool visited = false;
      if (!visited) 
      {
        visited = true;
        std::cerr << "Not implemented: Lagrange cube shape function 2nd derivative at " << __FILE__ << ":" << __LINE__ << "\n";
      }
      return 0;
    }
  };

  // static data member must be defined in order to take it's address (i.e. in std::make_tuple)
  template <class ctype, int dimension, class Scalar, int O>
  int const LagrangeCubeShapeFunction<ctype,dimension,Scalar,O>::order;

  //---------------------------------------------------------------------

  template <class ctype, int dimension, class Scalar, int Ord>
  class LagrangeCubeShapeFunctionSet: public ShapeFunctionSet<ctype,dimension,Scalar>
  {
  public:
    typedef LagrangeCubeShapeFunction<ctype,dimension,Scalar,Ord> value_type;
    typedef typename ShapeFunctionSet<ctype,dimension,Scalar>::Matrix Matrix;    

    LagrangeCubeShapeFunctionSet();

    virtual void interpolate(typename ShapeFunctionSet<ctype,dimension,Scalar>::SfValueArray const& A,
                             Matrix& IA) const
    {
      IA = A;
      // For Lagrange shape functions, the interpolation at their nodes is just the identity. Copy!
    }


    virtual value_type const& operator[](int i) const{ return sf[i]; }
    virtual Dune::GeometryType type() const { return Dune::GeometryType(Dune::GeometryType::cube,dimension); }

  private:
    std::vector<value_type> sf;
  };

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  /**
   * \tparam dim docme \warning what's that separate dimension definition good for?
   */
  template <class ctype, int dimension, class Scalar, int dim = dimension, bool restricted=false>
  class LagrangeShapeFunctionSetContainer: public ShapeFunctionSetContainer<ctype,dimension,Scalar>
  {
  public:
    typedef ShapeFunctionSet<ctype,dimension,Scalar> value_type;

  private:
    typedef std::map<std::pair<Dune::GeometryType,int>,value_type const*> Container;

  public:

    LagrangeShapeFunctionSetContainer()
    {
      // Currently only shape functions for simplices and cubes are defined.
      Dune::GeometryType simplex(Dune::GeometryType::simplex,dim);
      Dune::GeometryType cube(Dune::GeometryType::cube,dim);

      // Create zero order shape functions on the simplex, whithout hierarchic projection (of course).
      // For higher orders, the container is filled on demand.
      auto* sfs = new typename SimplexLagrangeDetail::ChooseShapeFunctionSet<ctype,dimension,Scalar,restricted>::type(0);
      container[std::make_pair(simplex,0)] = sfs;
      sfs->initHierarchicalProjection(0);

      container.insert(typename Container::value_type(std::make_pair(cube,0),new LagrangeCubeShapeFunctionSet<ctype,dimension,Scalar,0>()));
      container.insert(typename Container::value_type(std::make_pair(cube,1),new LagrangeCubeShapeFunctionSet<ctype,dimension,Scalar,1>()));
      container.insert(typename Container::value_type(std::make_pair(cube,2),new LagrangeCubeShapeFunctionSet<ctype,dimension,Scalar,2>()));
      // cube shape functions up to now only to second order implemented!
    }


    virtual ~LagrangeShapeFunctionSetContainer()
    {
      for (typename Container::iterator i=container.begin(); i!=container.end(); ++i)
        delete i->second;
    }

    virtual value_type const& operator()(Dune::GeometryType type, int order) const
    {
      std::lock_guard<std::mutex> lock(mut);
      
      auto i = container.find(std::make_pair(type,order));
      if(i==container.end()) 
        return createShapeFunctionSet(type,order);
      if(i->second->size()<0) 
        throw DetailedException("Invalid shape function set.",__FILE__,__LINE__);
      return *i->second;
    }

  private:
    // A map storing Lagrange shape function sets associated with
    // (reference element type, order). Due to runtime polymorphism,
    // pointers are stored in the map. This class owns the shape
    // function sets pointed to, and has to delete them in the
    // destructor.
    mutable Container container;
    mutable std::mutex mut;
    
    value_type const& createShapeFunctionSet(Dune::GeometryType simplex, int order) const
    {
      if (order<0)
        throw LookupException("Shape function sets of negative order not provided.",__FILE__,__LINE__);
      if (! simplex.isSimplex())
        throw LookupException("Desired FE order not available fon non-simplices.",__FILE__,__LINE__);
      
      // Create the requested shape function set
      auto* sfs = new typename SimplexLagrangeDetail::ChooseShapeFunctionSet<ctype,dimension,Scalar,restricted>::type(order);
      container[std::make_pair(simplex,order)] = sfs;
      
      // For hierarchic projection, we need a lower order shape function set. Get (or create) it.
      auto it = container.find(std::make_pair(simplex,order-1));
      if (it==end(container))
        createShapeFunctionSet(simplex,order-1);
      it = container.find(std::make_pair(simplex,order-1));
      assert(it != end(container));

      // create hierarchic projection
      sfs->initHierarchicalProjection(it->second);
      
      return *sfs;
    }
  };


  /**
   * \brief Returns a Lagrange shape function set for given reference element
   * type and given polynomial order. Singleton function.
   * 
   * \tparam restricted docme!
   * \tparam dim docme! \warning what's that good for?
   */
  template <class ctype, int dimension, class Scalar=double, bool restricted=false, int dim=dimension>
  typename LagrangeShapeFunctionSetContainer<ctype,dimension,Scalar,dim,restricted>::value_type const&
  lagrangeShapeFunctionSet(Dune::GeometryType type, int order)
  {
    // Initialization should be thread-safe (probably C++11 standard, 6.7/4).
    static LagrangeShapeFunctionSetContainer<ctype,dimension,Scalar,dim,restricted> container;
    return container(type,order);
  }
} // end of namespace Kaskade


#endif
