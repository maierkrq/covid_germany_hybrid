/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * \defgroup concepts Concepts
 * \brief STL-like concepts for types to be used as template parameters
 *
 * The concepts formulate and document the requirements for certain nontrivial types that can be used as template
 * parameters in several parts of Kaskade.
 */

/// \cond internals
/// \defgroup implementation Implementation Details
/// \endcond

/** \defgroup problemdefinition Problem Definition
 *  \brief Classes, needed to define a problem
 *
 *  This group of classes helps to define minimization problems or
 *  operator equations in finite element spaces.  for this purpose one
 *  has to define a variational functional, or an operator in weak
 *  form, which maps certain variables to certain residual spaces
 *  (often dual spaces).  Often one has to deal with systems of
 *  equations.
 *
 *  In order to set up a PDE problem, the following steps are usually required:
 *
 *  1. Obtain a grid and let a \ref GridManager handle subsequent mesh refinements and prolongation.
 *
 *  2. Define suitable finite element spaces required for the discretization of the problem (see \ref fem).
 *
 *  3. Define the variables ocurring in the problem. For equations (either scalar or vectorial), this will
 *     usually be only one. For systems of equations, there will be several variables (see \ref variables).
 *
 *  4. Define the variational functional (or the operator) and its derivatives (see \ref functional
 *     and the \ref VariationalFunctional concept).
 *
 *  5. Instantiate an assembler that can perform assembly of right hand sides and stiffness matrices
 *     (see \ref assembly).
 *
 * The following tasks are supported by the following classes
 *
 *  - Standard FEM Spaces: \ref DiscontinuousLagrangeMapper, \ref ContinuousLagrangeMapper, \ref NedelecMapper
 *  - Definition of Spaces: FEFunctionSpace
 *  - Definition of elementary properties of varibles: VariableDescription
 *  - Grouping several variables and spaces together into a whole: VariableSetDescription
 *    via \ref fusion "boost::fusion::vector"
 *  - Representation of such a group: VariableSet
 *  - Representation of a single component: FunctionSpaceElement
 *  - Definition of a functional/equation: VariationalFunctional
 *
 */

/**
 * \ingroup problemdefinition
 * \defgroup diffops Differential operators
 * \brief Predefined standard differential operators to be used as building blocks.
 *
 * This group of convenience classes helps to define variational and weak problems. It contains
 * a number of common differential operators (and their derivatives), e.g. for diffusion problems
 * or for elasticity.
 */

/**
 * \ingroup diffops
 * \defgroup diffopsElasto Elastomechanics
 * \brief Predefined standard differential operators and material laws to be used as building blocks for elastomechanics problems.
 */

/**
 * \ingroup diffopsElasto 
 * \defgroup contact Contact Problems
 * \brief Components for defining elastomechanical contact problems
 */

/**
 * \ingroup diffops
 * \defgroup boundaryconditions Boundary conditions
 * \brief Predefined standard boundary conditions to be used in boundary caches.
 */

/**
 * \ingroup problemdefinition
 * \defgroup variables Variables
 *
 * \brief Classes to define groups of variables
 *
 * Defining variables and sets of variables
 * ----------------------------------------
 *
 * In a PDE system, several variables occur. In the Stokes equation, e.g., velocity and pressure. In optimal
 * control, it is the state \f$ y \f$, the control \f$ u \f$ , and one or more multipliers, e.g., \f$ \lambda \f$.
 * These variables can live in different spaces, or in the same. In optimal control of elliptic equations, it is
 * usually \f$ y,\lambda \in H^1(\Omega) \f$ and \f$ u \in L^2(\Omega)\f$.
 *
 * For setting up the PDE problem, we need to define which variables occur, and in which space each lives. First
 * we define the spaces that are involved:
 * \code
 * auto spaces = makeSpaceList(&h1Space,&l2Space);
 * \endcode
 * Next we define the variables that occur, and tell which space they live in by indexing into the space list.
 * Additionally, the number of components are required, and a name can be given:
 * \code
 * auto variables = boost::fusion::make_vector( Variable<SpaceIndex<0>,Components<3>>("y"),
 *                                              Variable<SpaceIndex<1>,Components<1>>("u"),
 *                                              Variable<SpaceIndex<0>,Components<3>>("lambda"));
 * \endcode
 * As a treat, you don't have to remember whether SpaceIndex or Components comes first - they can be switched.
 *
 * And finally, both are tied together:
 * \code
 * auto variableSetDescription = makeVariableSetDescription(spaceList,variables);
 * \endcode
 * This variable set description tells the assembler how to treat the problem, and usually tells the
 * problem definition the finite element spaces in which the variables live.
 *
 * Accessing variables and associated spaces
 * -----------------------------------------
 * At some point usually one needs the index or type of the space to which a variable in a variable set
 * belongs. This is provided by \ref spaceIndex and \ref Space, respectively:
 * \code
 * int spaceIdx = spaceIndex<VariableSetDescription,varIdx>;
 * using SpaceType = Space<VariableSetDescription,varIdx>;
 * \endcode
 * In domain and boundary caches, all or some of the variables need to be evaluated from a sequence of evaluators,
 * usually in the \ref evaluateA method. This is provided by the \a value and \a derivative
 * methods of the variable set:
 * \code
 * auto u = vars.template value<uIdx>(evaluators);
 * auto du = vars.template derivative<uIdx>(evaluators);
 * \endcode
 *
 *
 */

/** \ingroup problemdefinition
 * \defgroup functional Functionals
 * \brief Classes for defining PDE problems in terms of variational functionals or weak formulations
 *
 * PDE problems in Kaskade7 are formulated in terms of variational functionals or weak formulations. In order 
 * to assemble right hand sides and matrices, an appropriate class modeling the  \ref VariationalFunctional concept
 * must be provided. Adaptors and helper classes supporting the definition of such classes are gathered here.
 * 
 * In particular, the definition of nonlinear or time-dependent PDEs is simplified using
 * \ref Kaskade::LinearizationAt or \ref Kaskade::SemiLinearizationAt in combination with the
 * \ref NonlinearVariationalFunctional concept.
 */

/**
 * \ingroup problemdefinition
 * \defgroup fem Finite Elements
 *
 * \brief Classes to define finite element spaces and functions on them
 *
 * In this group the problem is addressed of creating and maintaining finite element spaces and functions that live there.
 * The following tasks are fulfilled:
 *
 * - Finite Element Space: FEFunctionSpace
 * - Functions on a FEMSpace: FunctionSpaceElement
 * - Grid consistency issues: GridManager
 *
 * Objects of these classes are connected by a boost::signal in case that the grid is modified (adapted)
 *
 * Standard FEM spaces are implemented:
 *
 * - DiscontinuousLagrangeMapper
 * - DiscontinuousHierarchicMapper
 * - ContinuousLagrangeMapper
 * - ContinuousHierarchicMapper
 * - NedelecMapper
 * - ConstantMapper
 *
 * Access to shape functions
 *
 * - ShapeFunction
 * - ShapeFunctionSet
 *
 */

/**
 * \ingroup fem
 * \defgroup fefunctions Working with finite element functions
 * \brief Tools for representing and working with finite element functions.
 */

/**
 * \ingroup fefunctions
 * \defgroup fetransfer Function interpolation and transfer
 * \brief Tools for  converting finite element functions.
 */

/**
 * \ingroup fetransfer
 * \defgroup functionviews
 * \brief Predefined function views.
 */

/**
 * \ingroup fem
 * \defgroup grid Grid Management
 * \brief Tools for managing grids and their refinement, and for iterating over grids.
 */


/** \defgroup utilities Utilities
 *  \brief Miscellaneous support routines and classes
 */

/** \ingroup utilities
 * \defgroup combinatorics Combinatorics
 * \brief Functions for counting and enumerating combinatorial structures.
 */

/** \ingroup utilities
 * \defgroup threading Multithreading
 * \brief Support routines and data structures for multithreaded execution.
 */

/** \ingroup utilities
 * \defgroup exceptions Exceptions
 * \brief A hierarchy of exception classes providing more or less detailed information.
 */

/** \defgroup solution Solution Tools
 */

/**
 * \ingroup solution
 * \defgroup assembly Assembly
 * \brief Classes and functions for constructing matrix representations of differential operators.
 */

/** \ingroup solution
 *  \defgroup linalg Linear Algebra
 *  \brief Classes and functions for linear algebra.
 */

/**
 * \ingroup linalg
 * \defgroup linalgbasic Basic linear algebra
 * \brief Classes and functions for basic vector and matrix arithmetics
 */

/**
 * \ingroup linalg
 * \defgroup linalgsolution Solvers
 * \brief Classes and functions for solving linear equation systems
 */

/**
 * \ingroup linalgsolution
 * \defgroup direct Direct solvers
 * \brief Classes and functions for elimination methods
 */

/**
 * \ingroup linalgsolution
 * \defgroup iterative Iterative solvers
 * \brief Classes and functions for iterative solution of linear equation systems
 *
 * This includes on one hand actual \em solvers (Krylov methods such as \ref Kaskade::Pcg) and 
 * related functions (such as termination criteria),
 * and on the other hand \em preconditioners, both simple stationary iterations (such as Kaskade::JacobiPreconditioner) 
 * and multilevel preconditioners (such as Kaskade::AdditiveMultiGrid or Kaskade::MultiplicativeMultiGrid).
 */

/**
 * \ingroup iterative
 * \defgroup multigrid
 * \brief Classes and functions for multigrid preconditioners and solvers
 *
 * This includes general multiplicative V-cycles and additive multigrid methods, as well as several convenience functions for defining
 * concrete preconditioners by specifying particular ingredients (smoothers and coarse grid preconditioners, in particular).
 * Multigrid methods for arbitrary finite element order on simplicial grids are available.
 *
 * If unsure which multigrid method to create, try the default construction makeMultigrid, which tries to come up with
 * reasonable ingredients.
 */

 /**
  * \ingroup iterative
  * \defgroup preconditioners
  * \brief Preconditioners for iterative solvers, including smoothers for multigrid methods
  */



/** \ingroup solution
 *  \defgroup adapt Adaptivity
 *  \brief Tools for error estimation and local mesh refinement.
 *
 *  The issue of adaptivity is currently in a rather experimental state with several ideas coexisting.
 *
 */

/** \ingroup adapt
 * \defgroup refcrit Refinement Criteria
 * \brief Tools for selecting cells for refinement.
 *
 * Given an error indicator for each cell, the question remains as to which cells actually to refine. It is
 * obvious that cells with larger error indicator get refined whereas cells with smaller one are not marked,
 * so the question boils down to selecting a \em threshold for error indicators. Several choices have been
 * proposed in the literature, and some of them are implemented here.
 */

/**
 * \ingroup solution
 * \defgroup timestepping Timestepping
 *
 * This interface allows to define and solve time-dependent problems,
 * in particular parabolic PDEs.
 */

/** \ingroup solution
 *  \defgroup alg Algorithms
 */

/** \ingroup alg
 *  \defgroup abstract Abstract Algorithmic Objects
 *
 *  Here the basic elements are represented as abstract base classes,
 *  on which algorithms with dynamic polymorphism can be built.
 */

/**
 * \defgroup IO IO
 * \brief Input and Output functionality
 */

/**
 * \ingroup IO
 * \defgroup gridInput Grid input
 * \brief Grid creation or input
 */

/**
 * \ingroup IO
 * \defgroup entropycoding Entropy coding
 * \brief Tools for entropy coding, i.e. storing sequences of symbols with different frequencies in a (small) bitstream.
 */

/** \mainpage
 *
 *  Based on the <A HREF="http://www.dune-project.org"> DUNE Library
 *  </A>, Kaskade 7 is a toolbox for the solution of linear and
 *  nonlinear function space problems. Further, it makes heavy use of the <A
 *  HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/index.html">
 *  Boost Fusion Library </A>.
 *
 *  Some important groups of classes of Kaskade
 *
 *  - \ref problemdefinition "Definition of problems"
 *  - \ref solution "Tools for the solution"
 *
 */

/** \page fusion Some Remarks on Boost Fusion
 *
 * The <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/index.html"> boost fusion library </A> supports the use of
 *  <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/container.html"> heterogenous containers </A> or
 *  <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/view.html"> views </A>
 * implementing <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/sequence.html"> sequence concepts </A> and
 * <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/algorithm.html"> algorithms </A>,
 * working on them.
 * In Kaskade these containers,
 * in particular boost::fusion::vector is used to define a VariableSet, which consists of variables of varying types.
 *
 * There are two types of "functions" avaliable in boost::fusion:
 *
 * - proper functions: they work on data, placed inside the containers
 * - metafunctions: they work on the types of the data
 *
 * Usually each function is accompanied by a metafunction, which
 * returns the return type of the function.  To be able to store the
 * output of functions, one has to use metafunctions to determine the
 * storage type.
 *
 * Example: a function f that operates on a container v of type V and generates a container t of type T
 * \code
 * typedef T result_of::f<Arg,V>::type
 * T t = f<Arg>(v)
 * \endcode
 *
 * \b Some \b often \b used \b functions
 *
 *
 * - <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/container/vector.html"> vector </A> a heterogenous vector
 * - <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/sequence/intrinsic/functions/at_c.html"> at_c<N> </A>
 * returns a reference to the N'th element of a vector.
 *    Metafunctions: <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/doc/html/fusion/sequence/intrinsic/metafunctions/at_c.html"> at_c </A>
 * (type of the reference),
 *    <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/sequence/intrinsic/metafunctions/value_at_c.html">
 * value_at_c </A> (type, referenced on)
 * - <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/doc/html/fusion/algorithm/iteration/functions/for_each.html"> for_each </A> apply a
functor to every element of a container.
 * Metafunction: <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/doc/html/fusion/algorithm/iteration/functions/for_each.html"> for_each </A>
 * -  <A HREF="http://www.boost.org/doc/libs/1_59_0/libs/fusion/doc/html/fusion/algorithm/transformation/functions/transform.html"> transform </A>
 *
 */




 /**
  * \defgroup problems Specific Application Problems
  * \brief Classes that define and solve concrete problems.
  */

