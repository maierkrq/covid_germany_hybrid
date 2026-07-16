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

/**
 * @file
 * @brief  Documentation for Variational Functional Concepts
 * @author Martin Weiser, Anton Schiela
 */
class unspecified;

/**
 * \brief Documentation of the concept of a quadratic variational functional
 * \ingroup concepts
 * The variational functional concept defines the interface that is accessed by the VariationalFunctionalAssembler. 
 * 
 * This is *not* a base class, only a documentation! A convenience base class defining more or less useful default 
 * values for several required member types and functions is \ref FunctionalBase.
 *
 * The mathematical concept represented by this interface is that of a
 * variational functional
 * \f[ J(u) = \int_\Omega f(u_1(x),\nabla u_1(x),\dots,u_n(x),\nabla u_n(x)) \, dx
 *          + \int_{\partial\Omega} g(u_1(x),\dots,u_n(x)) \, ds
 *          + \sum_{F_T\in\mathcal{F}} g_F(u_1(x^+),\dots,u_n(x^+),u_1(x^-),\dots,u_n(x^-)) \, ds
 * \f]
 * with quadratic functions \f$f\f$, \f$g\f$, and \f$g_F\f$ to be evaluated at \f$0\f$.
 * The different variables \f$u_i\f$ may be scalar or vector-valued, and may live in different
 * finite element function spaces. \f$ \mathcal{F} \f$ is a (possibly empty) subset of
 * cell-face incidences in the grid, and \f$x^+\f$ and \f$x^-\f$ refer to the same point on
 * different sides of the face. This way, jumps of discontinuous ansatz functions can be
 * integrated on faces.
 *
 * Some local classes have to be written into VariationalFunctional,
 * where most of the functionality of the functional is defined. These
 * are:
 *
 * - DomainCache   (for evaluation of \f$f\f$ and its derivatives)
 * - BoundaryCache (for evaluation of \f$g\f$ and its derivatives)
 * - InnerBoundaryCache (optional, for evaluation of \f$g_F\f$ and its derivatives)
 * - D1 (information on the block structure of \f$ f' \f$ and \f$g'\f$)
 * - D2 (information on the block structure of \f$f''\f$ and \f$g'\f$)
 *
 * Also read the documentation of LinearizationAt in
 * functional_aux.hh, an adapter class for the usage and
 * the definition of a variational functional, and the NonlinearVariationalFunctional
 * concept. If only quadratic functionals (linear problems) are to be defined,
 * these classes are not necessary.
 *
 * Const methods are required to be thread-safe.
 * Note that the space lists AnsatzVars::Spaces and TestVars::Spaces must coincide!
 *
 * \see LinearizationAt
 * \see NonlinearVariationalFunctional
 */
class VariationalFunctional 
{
public:
  /**
   * \brief The scalar type to be used for this variational problem.
   */
  typedef unspecified Scalar;

  /**
   * \brief A description of the ansatz variables occuring in the variational problem. 
   * 
   * This should be an instance of VariableSetDescription<...>.
   */
  typedef unspecified AnsatzVars;


  /**
   * \brief A description of the test variables occuring in the "variational" problem. 
   * 
   * This should be an instance of VariableSetDescription<...>.
   */
  typedef unspecified TestVars;

  /**
   * \brief A description of the variables defining the linearization point.
   * 
   * In most cases, this is simply AnsatzVars.
   * 
   * In some cases, however, the linear(ized) problem to be solved is posed in a different 
   * ansatz space than the current evaluation point, e.g. hierarchical error 
   * estimators computing corrections in an extension space. Then there is a 
   * difference between ansatz and test variables on one hand and the origin 
   * variables on the other. 
   * 
   * This should be an instance of VariableSetDescription<...>.
   */
  typedef unspecified OriginVars;

  /**
   * \brief The type of problem, either VariationalFunctional or WeakFormulation.
   */
  static ProblemType const type;

  /**
   * The following typedef should be useful, when accessing Entities,
   * it is however not strictly required for a functional to work
   */
  typedef typename Vars::Grid::template Codim<0>::Entity Entity;

  /**
   * \brief This evaluation cache concept defines the interface that is
   * accessed by the assembler. 
   *
   * In the variational functional
   * \f[ J(u) = \int_\Omega f(u_1(x),\nabla u_1(x),\dots,u_n(x),\nabla u_n(x)) \, dx
   *          + \int_{\partial\Omega} g(u_1(x),\dots,u_n(x)) \, dx, \f]
   * the domain cache is responsible for evaluating the integrand \f$ f \f$ and
   * its derivatives. The VariationalFunctionalAssembler performs the integration cell
   * by cell. For each cell it considers, it first calls moveTo with the current cell,
   * allowing the domain cache to conduct lengthy computations that depend only on
   * the cell, but not on the quadrature point inside the cell (e.g., piecewise constant
   * coefficients), and store the results for later multiple use.
   * Next, the assembler runs through all quadrature points inside the current cell, and
   * calls evaluateAt with the local position of the quadrature point. Again, the domain
   * cache has the opportunity to perform expensive calculations that depend on the
   * quadrature point, but not on the directions for which directional derivatives are
   * to be computed (e.g., spatially non-constant coefficients of the PDE). After that,
   * d0, d1, and d2 are called possibly several times with different variational arguments
   * specifying the directions into which directional derivatives are to be computed.
   *
   * Thus, moveTo and evaluateAt need only be implemented if there is something
   * worthwile to compute. Otherwise, one can rely on the do-nothing default
   * implementation of TrivialCacheBase, from which the domain cache can inherit.
   * 
   * This is *not* a base class, only a documentation!
   *
   * Different objects have to be independent w.r.t. parallel
   * execution.
   */
  struct DomainCache 
  {

    /**
     * \brief Construct all data that is constant in an entity. Default
     * implementation in EvalCacheBase: does nothing
     */
    void moveTo(Cell const&)
    {
    }

    /**
     * \brief Construct all data that is constant at one point. Default
     * implementation in EvalCacheBase: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the value \f$f(0)\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$f\f$ at \f$ 0 \f$ with respect to
     * the variable \f$ u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,TestVars::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the second directional derivative of \f$f\f$ at \f$ 0 \f$ with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$u_{\mathrm{col}}\f$ in direction of the test functions given
     * by `argT` and `argA`.
     *
     * For multi-component variables (e.g., displacement in solid mechanics), the given test and ansatz
     * variations `argT` and `argA` will usually be scalar. The returned matrix \f$R\f$ is then defined
     * entry-wise as \f[ R_{ij} = f''(0)[\mathrm{argT}e_i,\mathrm{argA}e_j], \f]
     * i.e. every return entry corresponds to a second directional derivative in directions where just
     * one component is changed in each.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };

  /**
   * \brief Evaluation of boundary conditions
   * 
   * This evaluation cache concept defines the interface that is
   * accessed by the assembler. This is *not* a base class, only a
   * documentation!
   *
   * Different objects have to be independent
   * w.r.t. parallel execution. 
   */
  struct BoundaryCache    
  {
    /**
     * \brief Construct all data that is constant in an entity. 
     * 
     * The intersection iterator type need not be specified, the method can be templated.
     */
    void moveTo(typename Vars::GridView::IntersectionIterator const&);

    /**
     * \brief Construct all data that is constant at one point. 
     * \param xi the local position within the intersection/face
     * \param evaluators a boost::fusion sequence of evaluators corresponding to the sequence of FE spaces
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& xi, Evaluators const& evaluators);

    /**
     * Returns the value \f$g(0)\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$g\f$ at \f$ 0 \f$,
     * with respect to
     * the variable \f$u_{\mathrm{row}}\f$ in direction of the test
     * function given by argT.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,TestVars::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& argT) const;

    /**
     * Returns the second directional derivative of 
     * \f$ g \f$ at \f$ 0 \f$, with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$ u_{\mathrm{col}}\f$ in direction of the test functions given
     * by argT and the ansatz functions, given by argA.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,TestVars::template Components<row>::m,AnsatzVars::template Components<row>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };

  /**
   * \brief Evaluates jump contributions at interior faces, suitable for discontinuous Galerkin methods (optional).
   * 
   * If this structure is defined, the assembler does assemble jump terms
   * \f[ \int_F g(u^{\rm left},u^{\rm right})\,ds \f]
   * at interior faces \f$ F \f$ of the grid.
   * If no InnerBoundaryCache class is defined, the assembler ignores these integrals.
   */
  struct InnerBoundaryCache
  {
    template <class FaceIterator>
    void moveTo(FaceIterator const&);
    
    /**
     * \brief Prepare the cache for evaluation at given quadrature point.
     * 
     * \param x the quadrature point in the local face coordinate system
     * \param evaluators evaluators for quantities on the current cell 
     * \param neighbourEvaluators evaluators for the cell on the opposite side of the face
     */
    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators, Evaluators const& neighbourEvaluators);

    /**
     * \brief Returns the value \f$g(0)\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$g\f$ at \f$ 0 \f$,
     * with respect to the variable \f$u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     * 
     * The test variational argument belongs to the current cell, i.e. the interior side of the face.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,TestVars::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& argT) const;

    /**
     * Returns the second directional derivative of 
     * \f$ g \f$ at \f$ 0 \f$, with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$ u_{\mathrm{col}}\f$ in direction of the test functions given
     * by argT and the ansatz functions, given by argA.
     * 
     * If centerCell is true, both test and ansatz variational arguments are associated on the current cell.
     * Otherwise, the ansatz argument argA belongs to the neighbour cell, i.e. the outside of the face.
     *
     * This method need only be defined for values of row/col for which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,TestVars::::template Components<row>::m,AnsatzVars::::template Components<row>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA, bool centerCell) const;
  };

  /**
   * \brief Provides static information about the right hand side blocks.
   * 
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row block of
   * the first derivative.
   */
  template <int row>
  struct D1
  {
    /**
     * \brief Provides static information about the subvector blocks of the right hand side.
     * 
     * Is true if that block is statically present, i.e. if \f$ f \f$
     * depends on variable \f$ u_{\mathrm{row}} \f$.
     */
    static bool const present;
    
    /**
     * \brief Determines the highest derivatives that need to be available for right hand side evaluation.
     */
    static int const derivatives;
  };

  /**
   * \brief Provides static information about the submatrix blocks of the stiffness block matrix.
   * 
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row-col block of
   * the second derivative.
   */
  template <int row, int col>
  struct D2 
  {
    /**
     * \brief Specifies the presence of a submatrix block.
     * 
     * Is true if that block is statically present, i.e. if variable
     * \f$ u_{\mathrm{row}} \f$ and \f$ u_{\mathrm{col}}\f$ are
     * nonlinearly coupled.
     *
     * For second derivatives (problem type is VariationalFunctional), which are assumed to be always
     * symmetric, only the lower triangular blocks need be present.
     * 
     * Announcing blocks as not present allows the assembler to save both memory and cput time
     * as such blocks are completely ignored.
     */
    static bool const present;

    /**
     * \brief Specifies whether the subblock is symmetric.
     * 
     * Should be true if the block is conceptually symmetric, i.e. if
     * EvaluationCache::d2 with row/col and arg1/arg2 exchanged
     * returns the same value. Note that this does not imply that the
     * Galerkin representation is symmetric, since using different
     * ansatz/test spaces are possible.
     * 
     * Announcing a block to be symmetric allows the assembler to save both memory and cpu time 
     * by accessing only the lower half of the matrix block.
     *
     * Note that if the problem type is VariationalFunctional, all diagonal blocks (row==col)
     * are implicitly assumed to be symmetric.
     */
    static bool const symmetric;

    /**
     * \brief If this flag is true (and symmetric==true), the assembler will enforce positive
     *        semidefiniteness of all local matrices by projecting them onto the cone of psd 
     *        matrices.
     * 
     * If this is true, the variational functional must have a member function
     * \code
     * template <int row, int col>
     * double makePositiveThreshold() const;
     * \endcode
     * that defines the relative limit \f$ \alpha \f$ to which "positivity" is enforced:
     * All eigenvalues \f$ \lambda \f$ are shifted (if necessary) to satisfy
     * \f$ \lambda \ge \alpha \max(|\lambda_{\max}|,|\lambda_{\min}|) \f$.
     * The default implementation in the base class FunctionalBase uses \f$ \alpha = -10\epsilon \f$,
     * with machine precision \f$ \epsilon \f$.
     *
     * While this destroys accuracy in the Taylor approximation, it can be very convenient in 
     * algorithms for nonconvex minimization.
     */
    static const bool makePositive;

    /**
     * \brief Specifies whether only the diagonal of the subblock shall be assembled.
     * Is true if only the diagonal of the second derivative shall be
     * assembled. This is usually false, but can be set to true
     * e.g. for hierarchical error estimation.
     */
    static bool const lumped;
    
    /**
     * \brief Determines the highest derivatives that need to be available for evaluation of this block.
     */
    static int const derivatives;
  };
  
  /**
   * \brief Creates a DomainCache.
   * 
   * \param[in] flags a bit field describing what values will be accessed. See \ref Assembler enums.
   *                  This should be passed on to the domain cache, which then can restrict the 
   *                  computation to values that are actually needed.
   */
  DomainCache createDomainCache(int flags) const {}
  
  /**
   * \brief Creates a BoundaryCache.
   * 
   * \param[in] flags a bit field describing what values will be accessed. See \ref Assembler enums.
   *                  This should be passed on to the boundary cache, which then can restrict the 
   *                  computation to values that are actually needed.
   */
  BoundaryCache createBoundaryCache(int flags) const {}
  
  /**
   * \brief Creates an InnerBoundaryCache.
   *  
   * \param[in] flags a bit field describing what values will be accessed. See \ref Assembler enums.
   *                  This should be passed on to the inner boundary cache, which then can restrict the 
   *                  computation to values that are actually needed.
   * 
   * This is optional and has to be implemented only if the InnerBoundaryCache structure is defined.
   */
  InnerBoundaryCache createInnerBoundaryCache(int flags) const {}


  /**
   * \brief This method defines a suitable quadrature order to be used on the
   * given cell. 
   * 
   * The maximum shape function order occuring in the
   * finite element spaces for this cell is provided for
   * convenience. If boundary is true, this refers to evaluation of
   * \f$g\f$, otherwise to the evaluation of \f$f\f$.
   */
  int integrationOrder(typename Variables::Grid::template Codim<0>::Entity const& cell,
                       int shapeFunctionOrder, bool boundary) const;

  /**
   * This method tells on which inner cell-face incidences \f$ F_T \f$ of the grid view
   * the integration of the
   * variational functional shall be performed, if inner faces shall be treated at all,
   * i.e. an InnerBoundaryCache class is present. The stiffness matrix entries are only
   * allocated for those faces that are reported true, and the assembly is restricted
   * to those cells. Thus, a safe but maybe inefficient version is to always return true,
   * which is the default implementation in FunctionalBase.
   *
   * The method needs to be implemented only if the class InnerBoundaryCache is defined.
   */
  bool considerFace(Face<typename Variables::GridView> const& face) const;

};

// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------

/**
 * \ingroup concepts
 * \brief Documentation of the concept of a nonlinear variational functional
 * The nonlinear variational functional concept defines the interface that is accessed by the LinearizationAt
 * and linearizationAt adapters. It differs from the basic VariationalFunctional concept mainly in how the
 * DomainCache and the BoundaryCache are constructed.
 * 
 * This is *not* a base class, only a documentation!
 *
 * The mathematical concept represented by this interface is that of a
 * variational functional
 * \f[ J(u) = \int_\Omega f(u_1(x),\nabla u_1(x),\dots,u_n(x),\nabla u_n(x)) \, dx 
 *     + \int_{\partial\Omega} g(u_1(x),\dots,u_n(x)) \, dx \f]
 * to be evaluated for a specific fixed value of \f$u\f$.
 * The different variables \f$u_i\f$ may be scalar or vector-valued, and may live in different finite element
 * function spaces.
 *
 * Some local classes have to be written into VariationalFunctional,
 * where most of the functionality of the functional is defined. These
 * are:
 *
 * - DomainCache   (for evaluation of \f$f\f$)
 * - BoundaryCache (for evaluation of \f$g\f$)
 * - D1 (information on the block structure of \f$ f' \f$)
 * - D2 (information on the block structure of \f$f''\f$)
 *
 * Also read the documentation of \ref LinearizationAt in
 * functional_aux.hh, which provides helper classes for the usage and
 * the definition of a variational functional. If only quadratic
 * functionals (linear problems) are to be defined, these classes are
 * not necessary.
 *
 * Const methods are required to be thread-safe.
 * Note that the space lists AnsatzVars::Spaces and TestVars::Spaces must coincide!
 */
class NonlinearVariationalFunctional 
{
public:
  /**
   * \brief The scalar type to be used for this variational problem.
   */
  typedef unspecified Scalar;

  /**
   * \brief A description of the ansatz variables occuring in the
   * variational problem. This should be an instance of
   * VariableSetDescription<...>.
   */
  typedef unspecified AnsatzVars;


  /**
   * A description of the test variables occuring in the "variational"
   * problem. This should be an instance of
   * VariableSetDescription<...>.
   */
  typedef unspecified TestVars;
  
  /**
   * \brief A description of the type of the point of linearization
   */
  typedef unspecified OriginVars;

  /**
   * \brief The type of problem, either VariationalFunctional or
   * WeakFormulation.
   */
  static ProblemType const type;

  /**
   * The following typedef should be useful, when accessing Entities,
   * it is however not strictly required for a functional to work
   */
  typedef typename Vars::Grid::template Codim<0>::Entity Entity;

  /**
   * \brief This evaluation cache concept defines the interface that is accessed by the assembler. 
   * 
   * This is *not* a base class, only a documentation!
   *
   * Different objects have to be independent w.r.t. parallel execution.
   */
  struct DomainCache : public EvalCacheBase    
  {

    /**
     * Constructs an evaluation cache from the associated functional
     * f_, a point of linearization u_, and flags from the assembler (see VariationalFunctionalAssembler)
     * 
     * \param[in] flags A bit field that tells the cache what parts are to be assembled. If flags & Assembler::RHS
     *                  is false, d1() will never be called. If flags & Assembler::VALUE is false, d0() will never
     *                  be called, and d2() will never be called if flags & Assembler::MATRIX is false. This 
     *                  information can be used to restrict the computation to the actually required quantities.
     */
    DomainCache(Functional const& f, typename Vars::VariableSet const& u, int flags)  {};

    /**
     * Construct all data that is constant in an entity. Default
     * implementation in EvalCacheBase: does nothing
     */
    template<class Entity>
    void moveTo(Entity const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default
     * implementation in EvalCacheBase: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the value \f$f(u(x))\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$f\f$  with respect to
     * the variable \f$u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the second directional derivative of \f$f\f$ with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$u_{\mathrm{col}}\f$ in direction of the test functions given
     * by \arg argT and \arg argA.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template Components<row>::m,Variables::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };

  /**
   * \brief Defining boundary conditions.
   * 
   * This evaluation cache concept defines the interface that is
   * accessed by the assembler. This is *not* a base class, only a
   * documentation!
   *
   * Different objects have to be independent
   * w.r.t. parallel execution. 
   */
  struct BoundaryCache : public EvalCacheBase    
  {
    /**
     * \brief Type of the Face Information
     */
    typedef typename AnsatzVars::Grid::LeafIntersectionIterator FaceIterator;

    /**
     * Constructs an evaluation cache from the associated functional f_, a point of 
     * linearization u_, and flags from the assembler 
     * 
     * \param[in] flags A bit field that tells the cache what parts are to be assembled. If flags & Assembler::RHS
     *                  is false, d1() will never be called. If flags & Assembler::VALUE is false, d0() will never
     *                  be called, and d2() will never be called if flags & Assembler::MATRIX is false. This 
     *                  information can be used to restrict the computation to the actually required quantities.
     */
    BoundaryCache(Functional const& f, typename Vars::VariableSet const& u, int const flags=7)  {};

    /**
     * Construct all data that is constant in an entity. Default implementation in EvalBaseCache: does nothing
     */
    template<class FaceIterator>
    void moveTo(FaceIterator const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default implementation in EvalBaseCache: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the value \f$g(u(x))\f$
     */
    Scalar d0() const;

    /**
     * Returns the directional derivative of \f$g\f$,
     * with respect to
     * the variable \f$u_{\mathrm{row}}\f$ in direction of the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the second directional derivative of 
     * \f$g\f$, with
     * respect to the variables \f$u_{\mathrm{row}}\f$ and
     * \f$u_{\mathrm{col}}\f$ in direction of the test functions given
     * by argT and the ansatz functions, given by argA.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template Components<row>::m,Variables::template Components<col>::m>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };


  /**
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row block of
   * the first derivative.
   */
  template <int row>
  struct D1
  {
    /**
     * Is true if that block is statically present, i.e. if \f$ f \f$
     * depends on variable \f$u_{\mathrm{row}}\f$.
     */
    static bool const present;
  };

  /**
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row-col block of
   * the second derivative.
   */
  template <int row, int col>
  struct D2 
  {
    /**
     * Is true if that block is statically present, i.e. if variable
     * \f$u_{\mathrm{row}}\f$ and \f$u_{\mathrm{col}}\f$ are
     * nonlinearly coupled.
     *
     * For second derivatives, which are assumed to be always
     * symmetric, only the lower triangular blocks need be present.
     */
    static bool const present;

    /**
     * Is true if the block is conceptually symmetric, i.e. whether
     * EvaluationCache::d2 with row/col and arg1/arg2 exchanged
     * returns the same value. Note that this does not imply that the
     * Galerkin representation is symmetric, since using different
     * ansatz/test spaces are possible.
     */
    static bool const symmetric;

    /**
     * Is true if only the diagonal of the second derivative shall be
     * assembled. This is usually false, but can be set to true
     * e.g. for hierarchical error estimation.
     */
    static bool const lumped;
  };

  /**
   * This method defines a suitable quadrature order to be used on the
   * given cell. The maximum shape function order occuring in the
   * finite element spaces for this cell is provided for
   * convenience. If boundary is true, this refers to evaluation of
   * \f$g\f$, otherwise to the evaluation of \f$f\f$.
   */
  int integrationOrder(typename Variables::Grid::template Codim<0>::Entity const& cell,
                       int shapeFunctionOrder, bool boundary) const;

  /**
   * This method tells on which inner cell-face incidences \f$ F_T \f$ of the grid view
   * the integration of the
   * variational functional shall be performed, if inner faces shall be treated at all,
   * i.e. an InnerBoundaryCache class is present. The stiffness matrix entries are only
   * allocated for those faces that are reported true, and the assembly is restricted
   * to those cells. Thus, a safe but maybe inefficient version is to always return true,
   * which is the default implementation in FunctionalBase.
   *
   * The method needs to be implemented only if the class InnerBoundaryCache is defined.
   */
  bool considerFace(Face<typename Variables::GridView> const& face) const;
};

// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------

/**
 * \ingroup concepts
 * \brief Documentation of the concept of a parabolic parabolic equation
 * 
 * The ParabolicEquation concept defines the interface that is accessed by the 
 * SemiLinearizationAt and semiLinearization adapters. It differs from the basic 
 * VariationalFunctional concept mainly in how the DomainCache and the BoundaryCache 
 * are constructed.
 * 
 * This is *not* a base class, only a documentation!
 *
 * The mathematical concept represented by this interface is that of a
 * weak formulation for a parabolic equation
 * \f[ \int_\Omega v^T B \dot u  \,dx = \int_\Omega (v^T f(u,\nabla u)+ \nabla v : \tilde f(u,\nabla u) \, dx 
 *     + \int_{\partial\Omega} v^Tg(u) \, dx \quad \forall v\f]
 * to be evaluated for a specific fixed value of \f$u\f$ (which may consist of multiple variables).
 * The different variables \f$u_i\f$ and \f$ v_i \f$ may be scalar or vector-valued, and may 
 * live in different finite element function spaces.
 *
 * Some local classes have to be written into ParabolicEquation,
 * where most of the functionality of the equation is defined. These
 * are:
 *
 * - DomainCache   (for evaluation of \f$f\f$ and \f$\tilde f\f$)
 * - BoundaryCache (for evaluation of \f$g\f$)
 * - D1 (information on the block structure of \f$ f' \f$)
 * - D2 (information on the block structure of \f$f''\f$)
 *
 * Also read the documentation of \ref SemiLinearizationAt in
 * functional_aux.hh, which provides helper classes for the usage and
 * the definition of a variational functional. 
 *
 * Const methods are required to be thread-safe.
 * Note that the space lists AnsatzVars::Spaces and TestVars::Spaces must coincide!
 */
class ParabolicEquation 
{
public:
  /**
   * \brief The scalar type to be used for this variational problem.
   */
  typedef unspecified Scalar;

  /**
   * \brief A description of the ansatz variables occuring in the
   * weak formulation. This should be an instance of
   * VariableSetDescription<...>.
   */
  typedef unspecified AnsatzVars;


  /**
   * A description of the test variables occuring in the weak formulation. 
   * This should be an instance of VariableSetDescription<...>.
   */
  typedef unspecified TestVars;
  
  /**
   * \brief A description of the type of current evaluation state.
   */
  typedef unspecified OriginVars;

  /**
   * \brief The type of problem, either VariationalFunctional (for symmetric problems) 
   * or WeakFormulation.
   */
  static ProblemType const type;

  /**
   * The following typedef should be useful, when accessing Entities,
   * it is however not strictly required for a ParabolicEquation to work
   */
  using Cell = Kaskade::Cell<typename AnsatzVars::GridView>;

  /**
   * \brief This evaluation cache concept defines the interface that is accessed by the assembler. 
   * 
   * This is *not* a base class, only a documentation!
   *
   * Different objects have to be independent w.r.t. parallel execution.
   */
  struct DomainCache : public EvalCacheBase    
  {

    /**
     * Constructs an evaluation cache from the associated functional
     * f_, a point of linearization u_, and flags from the assembler (see VariationalFunctionalAssembler)
     * 
     * \param[in] u   evaluation point (used for b1, b2, and d1)
     * \param[in] flags A bit field that tells the cache what parts are to be assembled. 
     *                  If flags & Assembler::RHS is false, d1() will never be called. If 
     *                  flags & Assembler::VALUE is false, d0() will never
     *                  be called, and d2() will never be called if flags & Assembler::MATRIX 
     *                  is false. This information can be used to restrict the computation to 
     *                  the actually required quantities.
     */
    DomainCache(Functional const& f, typename OriginVars::VariableSet const& u, int flags)  {};

    /**
     * Construct all data that is constant on a cell. Default
     * implementation in EvalCacheBase: does nothing
     */
    void moveTo(Cell const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default
     * implementation in EvalCacheBase: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     * Returns the component \f$f_{\mathrm{row}}\f$ multiplied by the test
     * function given by arg.
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template components<row>> 
    d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     * Returns the directional derivative of \f$f_{\mathrm{row}}\f$ with
     * respect to the variables \f$u_{\mathrm{col}}\f$ in direction of the test 
     * functions given by \arg argA, and multiplied by the test function given
     * by \arg argT.
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template components<row>,
                             Variables::template components<col>>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;

    /**
     * Returns the product \f$v_T B_{\rm row,col} v_A\f$ with
     * the ansatz function \f$ v_A\f$ given by \arg argA, and the test function 
     * \f$ v_T \f$ given by \arg argT.
     *
     * This method need only be defined for values of row/col for
     * which B2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template components<row>,
                             Variables::template components<col>>
    b2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };

  /**
   * \brief Defining boundary conditions.
   * 
   * This evaluation cache concept defines the interface that is
   * accessed by the assembler. This is *not* a base class, only a
   * documentation!
   *
   * Different objects have to be independent
   * w.r.t. parallel execution. 
   */
  struct BoundaryCache : public EvalCacheBase    
  {
    /**
     * \brief Type of the Face Information
     */
    typedef typename AnsatzVars::Grid::LeafIntersectionIterator FaceIterator;

    /**
     * Constructs an evaluation cache from the associated functional f_, a point of 
     * linearization u_, and flags from the assembler 
     * 
     * \param[in] flags A bit field that tells the cache what parts are to be assembled. If flags & Assembler::RHS
     *                  is false, d1() will never be called. If flags & Assembler::VALUE is false, d0() will never
     *                  be called, and d2() will never be called if flags & Assembler::MATRIX is false. This 
     *                  information can be used to restrict the computation to the actually required quantities.
     */
    BoundaryCache(Functional const& f, typename OriginVars::VariableSet const& u, int const flags=7)  {};

    /**
     * Construct all data that is constant in an entity. Default implementation in EvalBaseCache: does nothing
     */
    template<class FaceIterator>
    void moveTo(FaceIterator const&) 
    {
    }

    /**
     * Construct all data that is constant at one point. Default implementation in EvalBaseCache: assert(false)
     */
    template<class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators) 
    {
    }

    /**
     */
    template <int row, int dim>
    Dune::FieldVector<Scalar,Variables::template Components<row>::m> d1(VariationalArg<Scalar,dim> const& arg) const;

    /**
     *
     * This method need only be defined for values of row/col for
     * which D2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template components<row>,
                             Variables::template Components<col>>
    d2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;

    /**
     *
     * This method need only be defined for values of row/col for
     * which B2<row,col>::present is true.
     */
    template <int row, int col, int dim>
    Dune::FieldMatrix<Scalar,Variables::template components<row>,
                             Variables::template Components<col>>
    b2(VariationalArg<Scalar,dim> const& argT, VariationalArg<Scalar,dim> const& argA) const;
  };


  /**
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row block of
   * the first derivative.
   */
  template <int row>
  struct D1
  {
    /**
     * Is true if that block is statically present, i.e. if \f$ f \f$
     * depends on variable \f$u_{\mathrm{row}}\f$.
     */
    static bool const present;
  };

  /**
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row-col block of
   * the second derivative.
   */
  template <int row, int col>
  struct D2 
  {
    /**
     * Is true if that block is statically present, i.e. if variable
     * \f$u_{\mathrm{row}}\f$ and \f$u_{\mathrm{col}}\f$ are
     * nonlinearly coupled.
     *
     * For second derivatives, which are assumed to be always
     * symmetric, only the lower triangular blocks need be present.
     */
    static bool const present;

    /**
     * Is true if the block is conceptually symmetric, i.e. whether
     * EvaluationCache::d2 with row/col and arg1/arg2 exchanged
     * returns the same value. Note that this does not imply that the
     * Galerkin representation is symmetric, since using different
     * ansatz/test spaces are possible.
     */
    static bool const symmetric;

    /**
     * Is true if only the diagonal of the second derivative shall be
     * assembled. This is usually false, but can be set to true
     * e.g. for hierarchical error estimation.
     */
    static bool const lumped;
  };

  /**
   * This block info concept defines the interface that is accessed by
   * the assembler. This is *not* a base class, only a documentation!
   *
   * The class provides static information about the row-col block of
   * the matrix factor \f$ B \f$ in front of the time derivative.
   */
  template <int row, int col>
  struct B2 
  {
    /**
     * \brief Presence flag.
     * Is true if that block is statically present, i.e. \f$ B_{rc} \ne 0 \f$.
     */
    static bool const present;

    /**
     * Is true if the block is conceptually symmetric, i.e. whether
     * EvaluationCache::b2 with row/col and arg1/arg2 exchanged
     * returns the same value. Note that this does not imply that the
     * Galerkin representation is symmetric, since using different
     * ansatz/test spaces are possible.
     */
    static bool const symmetric;

    /**
     * Is true if only the diagonal of the second derivative shall be
     * assembled. This is usually false, but can be set to true
     * e.g. for hierarchical error estimation.
     */
    static bool const lumped;
    
    /**
     * Often, the matrix factor \f$ B \f$ in front of the time derivative
     * is just a constant matrix, but it may depend on the solution
     * \f$ u \f$. In the former case, computations can be simplified,
     * which is enabled by setting this flag to true.
     */
    static bool const constant;
  };

  /**
   * \brief Reports the current simulated time for which the evaluations are done.
   */
  double time() const;

  /**
   * \brief Sets a new time for evaluating PDE coefficients.
   *
   * \return a mutable reference to the parabolic equation object
   */
  ParabolicEquation& setTime(double t);

  /**
   * \brief Notifies the PDE in which time interval evaluations will be done.
   *
   * This can be larger than the time step \f$ \tau \f$, e.g., in extrapolation
   * or deferred correction methods, where several shorter Euler steps are performed
   * and combined to a higher order time step.
   *
   * The PDE formulation may (or may not) make use of this information, e.g., by
   * bounding some terms within that time range.
   */
  void temporalEvaluationRange(double t0, double t1);

  /**
   * This method defines a suitable quadrature order to be used on the
   * given cell. The maximum shape function order occuring in the
   * finite element spaces for this cell is provided for
   * convenience. If boundary is true, this refers to evaluation of
   * \f$g\f$, otherwise to the evaluation of \f$f\f$.
   */
  int integrationOrder(typename Variables::Grid::template Codim<0>::Entity const& cell,
                       int shapeFunctionOrder, bool boundary) const;
};

// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------

/**
 * \ingroup functional
 *
 * \brief Concept for providing block information to hierarchical error estimator.
 * 
 * Concept that the fourth template parameter of the
 * HierarchicErrorEstimator must model. Note that this is no base
 * class, just an interface documentation.
 */
struct HierarchicErrorEstimatorD2Info 
{
  /**
   * Member template that models the VariationalFunctional::D2 concept.
   */
  template <int row, int col> class D2;
};
