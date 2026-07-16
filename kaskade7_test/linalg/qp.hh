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

#ifndef QP_HH
#define QP_HH

#include <memory>
#include <tuple>

#include "dune/istl/bvector.hh"

#include "linalg/dynamicMatrix.hh"
#include "linalg/threadedMatrix.hh"

namespace Kaskade
{
  /**
   * \ingroup iterative
   * \defgroup qp QP Solvers
   * \brief Classes and methods for solving linear-quadratic programs
   *
   * A general linear-quadratic program is of the (standard) form
   * \f[ \min_x \frac{1}{2} x^T A x + c^T x \quad \textsf{s.t. $Bx\le b$, $Cx=0$}. \f]
   *
   * The classes provided here are intended for subsets of QPs with certain structure, in particular
   * without equality constraints. There are classes intended to be used directly for solving QPs
   * and supporting classes.
   *
   * ## General inequalities
   *
   * For general inequality constrained problems of the type
   * \f[ \min_x \frac{1}{2} x^T A x + c^T x \quad \textsf{s.t. $Bx\le b$} \f]
   * with positive (semi-) definite \f$ A \f$, there are two solvers:
   * - \ref QPSolver, which transforms the QP into its Lagrange dual problem and solves this via a projected gradient method
   * - \ref QPPenalizedSolver, which solves a quadratically penalized problem
   * - \ref QPALSolver, which uses quadratic penalization and a Rockafellar augmentation
   * - \ref QPMultiGrid,
   *
   * ## Bound constraints
   *
   * For bound constrained problems of the form
   * \f[ \min_x \frac{1}{2} x^T A x + c^T x \quad \textsf{s.t. $x\le b$} \f]
   * with positive semi-definite \f$ A \f$, there are
   * - \ref QPBoundSolver, realizing a projected gradient method with dense linear algebra
   * - \ref QPDirectSparse, a projected gradient method with sparse linear algebra
   */



  enum class QPConvexificationStrategy { DONOTHING, INCREMENTALADDITION, EIGENVALUEADDITION };

  /**
   * \brief A flag that determines whether the QP solver shall work with dense or sparse linear algebra.
   */
  enum class QPStructure { DENSE, SPARSE };

  //-----------------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

  /**
   * \ingroup qp
   * \brief Computes the KKT residual for approximate solutions of certain QPs.
   *
   * For QPs of the form \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } Bx \le b, \f]
   * this function computes the KKT residual and returns the values
   * \f[ [ \|Ax + c + B^T \lambda \|, \|Bx-b\|_+ , \|\lambda\|_- , \lambda^T (b-Bx) ]. \f]
   */
  template <class MatrixA, class MatrixB, class VectorX, class VectorB>
  std::array<typename ScalarTraits<typename MatrixA::field_type>::Real,4>
  checkKKT(MatrixA const& A, MatrixB const& B,VectorX const& c, VectorB const& b,
           VectorX const& x, VectorB const& lambda);


  //-----------------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup qp
   * \brief An iterative solver for small instances of bound constrained quadratic programs.
   * 
   * The solver addresses problems of the form
   * \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } x \le b, \f]
   * where \f$ A \f$ is positive definite. The implementation  is intended to solve 
   * _small_ problems of just a couple of variables, typically less than 30.
   * 
   * The actual matrix operations have to be provided by a derived class.
   * 
   * \tparam Real the real field type of matrix/vector entries.
   * \tparam Implementation the type of the derived class (Barton-Nackman trick)
   * The implementation class needs to overwrite the solveSubmatrix method.
   * 
   * ### Solution algorithm ###
   * 
   * A simple projected gradient method is used.
   */
  template <class Real, class Implementation>
  class QPBoundSolverBase
  {
  public:
    using Vector = Dune::BlockVector<Dune::FieldVector<Real,1>>;
    using Matrix = DynamicMatrix<Dune::FieldMatrix<Real,1,1>>;

    virtual ~QPBoundSolverBase()
    {}
    
    /**
     * \brief Solves the QP for the given right hand side vectors.
     *
     * Solves the QP with a projected gradient method
     * starting at the feasible point \f$ x = \min(0,b)\f$.
     *
     * \param c the linear term of the objective
     * \param b the constraints right hand side
     * \param tol the tolerance for termination criterion
     *
     * The iteration is terminated as soon as the projected gradient \f$ s \f$ is small,
     * i.e. \f$ \|s\|_2 \le \mathsf{tol} \f$.
     *
     * \return a tuple \f$ (x,iter) \f$ of solution vector, and iteration count
     *
     * If the QP is infeasible, an InfeasibleProblemException is thrown.
     */
    std::tuple<Vector,int> solve(Vector const& c, Vector const& b, double tol) const;

    /**
     * \brief Solves the QP for the given right hand side vectors.
     *
     * \param x an initial guess for the solution. If an unfeasible starting point is provided, it will be
     *          projected to the feasible set.
     * \param c the linear term of the objective
     * \param b the constraints right hand side
     * \param tol the tolerance for termination criterion
     *
     * The iteration is terminated as soon as the projected gradient \f$ s \f$ is small,
     * i.e. \f$ \|s\|_2 \le \mathsf{tol} \f$.
     *
     * \return a tuple \f$ (x,iter) \f$ of solution vector, and iteration count. If the iteration count
     *         equals the maximum iteration count (\see setMaxSteps), the returned solution vector need
     *         not satisfy the tolerance criterion.
     *
     * If the QP is infeasible, an InfeasibleProblemException is thrown.
     */
    std::tuple<Vector,int> solve(Vector const& x, Vector const& c, Vector const& b, double tol) const;

    /**
     * \brief Performs a single projected gradient step.
     * 
     * This is a generalization of a Richardson smoother with linesearch, and can
     * be used when the QP is just a subproblem, e.g. in multigrid methods, and 
     * rather well conditioned.
     */
    void gradientStep(Vector const& c, Vector const& b, Vector& x) const;
    
    /**
     * \brief defines the minimal number of steps to perform in solving.
     *
     * The default value is 0. If the minimum iteration count is set to less than the maximum
     * iteration count, the number of steps performed is undefined.
     */
    Implementation& setMinSteps(int count);

    /**
     * \brief defines the maximum allowed number of steps to perform in solving.
     *
     * The default value is 200. If the maximum iteration count is below the minimum
     * iteration count, the number of steps performed in solving is undefined.
     */
    Implementation& setMaxSteps(int count);

  protected:
    /**
     * \brief Constructor (to be used by derived classes only)
     * \param dim the dimension of the primal variables, i.e. the length of vector \f$ x \f$
     */
    QPBoundSolverBase(int dim)
    : dimension(dim) 
    , minSteps(0)
    , maxSteps(200)
    {}
    
    /** 
     * \brief solves a subsystem \f$ \tilde A x = b \f$.
     * 
     * This needs to be overwritten by derived classes.
     * 
     * For a set \a idx of indices, the submatrix \f$ \tilde A \f$ in Matlab notation is A(idx,idx).
     * \param idx the indices defining the submatrix 
     * \param b the right hand side (of the same length as the index vector)
     * \return the solution \f$ x \f$
     */
    Vector solveSubmatrix(std::vector<int> const& idx, Vector const& b) const;

    /**
     * \brief Performs the update \f$ s \gets s + aAe_i \f$, where \f$ e_i \f$ is the i-th unit vector.
     *
     * The default implementation performs a complete matrix-vector product, which is
     * performance-wise suboptimal. Override this for better performance.
     */
    virtual void updateWithColumn(Real a, int i, Vector& s) const;
    
  private:
    int dimension;      // size of primal variables (i.e. length of x and size of A)
    int minSteps;       // minimal number of steps to perform when solving
    int maxSteps;       // upper bound on iteration count
    
    bool getProjectedGradient(Vector const& x, Vector const& c, Vector const& b, Vector& g, Vector& s) const;
    void getProjectedNewton(Vector const& x, Vector const& c, Vector const& b, Vector& g, Vector& s) const;
    
    /**
     * \brief Updates \f$ x \f$ by a line search along the projected search direction \f$ s \f$.
     * \param x starting point 
     * \param b bounds
     * \param g gradient at \f$ x \f$
     * \param s search direction
     * \return true if the active set changed
     */
    bool projectedLineSearch(Vector& x, Vector const& b, Vector& g, Vector& s,
      Vector const* c=nullptr
    ) const;
    void checkNaN(Vector const& v) const;
    
    virtual void printMatrix(std::ostream& out) const = 0;

  };
  
  //-----------------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

  /**
   * \ingroup qp
   * \brief An iterative solver for small instances of bound constrained quadratic programs.
   * 
   * The solver addresses problems of the form
   * \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } x \le b, \f]
   * where \f$ A \f$ is a positive definite dense matrix. The implementation  is intended to solve 
   * _small_ problems of just a couple of variables, typically less than 30.
   * 
   * \tparam Real the real field type of matrix/vector entries. Explicit instantiations are
   * provided for float and double.
   */
  template <class R=double>
  class QPBoundSolver: public QPBoundSolverBase<R,QPBoundSolver<R>>
  {
    using Base = QPBoundSolverBase<R,QPBoundSolver<R>>;
    
  public:
    using Real   = R;
    using Vector = typename Base::Vector;
    using Matrix = DynamicMatrix<Dune::FieldMatrix<Real,1,1>>;
    
    /**
     * \brief Constructs the solver, providing the matrix \f$ A \f$.
     * 
     * The argument matrix is copied internally, and is not referenced later.
     * 
     * Throws a NonpositiveMatrixException if \f$ A \f$ is not positive definite 
     * and \arg convex is QPConvexificationStrategy::DONOTHING.
     * 
     * \todo write a move constructor
     */
    QPBoundSolver(Matrix const& A, QPConvexificationStrategy convex=QPConvexificationStrategy::DONOTHING);
        
  private:
    Matrix A;

    friend Base;
    void mv(Vector const& x, Vector& Ax) const { A.mv(x,Ax); }
    Real frobenius_norm() const { return A.frobenius_norm(); }
    Vector solveSubmatrix(std::vector<int> const& idx, Vector const& b) const;
    
    virtual void printMatrix(std::ostream& out) const { out << "A = \n" << A; }
    virtual void updateWithColumn(Real a, int i, Vector& s) const override;
  };
  
  //-----------------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

  /**
   * \ingroup qp
   * \brief A solver for sparse, medium-sized instances of bound constrained quadratic programs.
   * 
   * The solver addresses problems of the form
   * \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } x \le b, \f]
   * where \f$ A \f$ is a positive definite sparse matrix. A direct solver is used to 
   * solve projected subproblems.
   * 
   * 
   * \tparam Real the real field type of matrix/vector entries. Explicit instantiations are
   * provided for float and double.
   */
  template <class R=double>
  class QPDirectSparse: public QPBoundSolverBase<R,QPDirectSparse<R>>
  {
    using Base = QPBoundSolverBase<R,QPDirectSparse<R>>;
     
    public:
      using Real = R;
      using Vector = typename Base::Vector;       // Dune::BlockVector<Dune::FieldVector<Real,1>>
      using SparseMatrix =  NumaBCRSMatrix<Dune::FieldMatrix<Real,1,1>>;
      
    /**
     * \brief Constructs the solver, providing the matrix \f$ A \f$.
     * 
     * The argument matrix is copied internally, and is not referenced later.
     * 
     * Throws a NonpositiveMatrixException if \f$ A \f$ is not positive definite 
     * and \arg convex is QPConvexificationStrategy::DONOTHING.
     * 
     * \todo write a move constructor
     */
      QPDirectSparse(SparseMatrix const& A, QPConvexificationStrategy convex=QPConvexificationStrategy::DONOTHING);
      
    private:
      SparseMatrix A;
      void mv(Vector const&x, Vector& Ax) const { A.mv(x,Ax); }
      Real frobenius_norm() const;
      friend Base;      
      Vector solveSubmatrix(std::vector<int> const& idx, Vector const& b) const;
      
      virtual void printMatrix(std::ostream& out) const { out << "sparse matrix A in triplet format: \n" << A; }
  };
  
  //--------------------------------------------------------------------------------------------------------------

  /**
   * \ingroup qp
   * \brief An iterative solver for particular instances of bound constrained quadratic programs.
   * 
   * The solver addresses problems of the form
   * \f[ \min_x \frac{1}{2} x^T Qx + c^T x \quad \text{s.t. } x \le b, \f]
   * where \f[ Q = \begin{bmatrix} A+\gamma B^T B & -\gamma B^T \\ -\gamma B & \gamma I \end{bmatrix}\f]
   * is a positive definite block
   * structured matrix. Problems of this type rarely arise on their own, mostly from a penalty or augmented
   * Lagrangian relaxation of general QPs (see \ref QPPenalizedSolver).
   *
   * The implementation  is intended to solve problems where \f$ A \f$ is _small_ and dense , with 
   * just a couple of variables, typically less than 30.
   * 
   * \tparam structure specifies which type of linear algebra to use:
   *                   - DENSE uses dense linear algebra and is intended to solve problems where \f$ A \f$
   *                     is small, with usually not more than 30 variables
   *                   - SPARSE uses sparse solvers, and is intended for problems of moderate size
   * \tparam Real the real field type of matrix/vector entries. Explicit instantiations are
   * provided for float and double.
   */
  template <QPStructure sparsity=QPStructure::DENSE, class R=double>
  class QPSlackSolver: public QPBoundSolverBase<R,QPSlackSolver<sparsity,R>>
  {
    using Base = QPBoundSolverBase<R,QPSlackSolver<sparsity,R>>;
    using Entry = Dune::FieldMatrix<R,1,1>;

  public:
    using Real = R;
    using Vector = typename Base::Vector;             // Dune::BlockVector<Dune::FieldVector<Real,1>>
    using Matrix = std::conditional_t<sparsity==QPStructure::DENSE,
                                      DynamicMatrix<Entry>,
                                      NumaBCRSMatrix<Entry>>;
    
    /**
     * \brief Constructs the solver, providing the matrix \f$ A \f$ and \f$ B \f$.
     * 
     * The argument matrices are copied internally, and are not referenced later.
     * 
     * Throws a NonpositiveMatrixException if \f$ A \f$ is not positive definite
     * and \arg convex is QPConvexificationStrategy::DONOTHING.
     *
     * \todo write a move constructor
     */
    QPSlackSolver(Matrix const& A, Matrix const& B, Real gamma, QPConvexificationStrategy convex=QPConvexificationStrategy::DONOTHING);
    
  private:
    Matrix A;
    Matrix B;
    Real gamma;

    friend Base;
    void mv(Vector const& x, Vector& Ax) const;
    Real frobenius_norm() const;
    Vector solveSubmatrix(std::vector<int> const& idx, Vector const& b) const;

    virtual void printMatrix(std::ostream& out) const { out << "A = \n" << A << "\nB = \n" << B << "\ngamma = " << gamma << "\n"; }
  };

  
  //--------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------

   /**
    * \ingroup qp
    * \brief An iterative solver for small instances of a particular class of quadratic programs.
    *
    * The solver addresses problems of the form
    * \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } Bx \le b, \f]
    * where \f$ A \f$ is positive definite. The implementation works with dense matrices
    * and is intended to solve _small_ problems of just a couple of variables, typically
    * less than 30, and moreover _well conditioned_ problems.
    *
    * \tparam Real the real field type of matrix/vector entries.
    * Explicit instantiations are provided for float and double.
    * \tparam d the dimension of (vector valued) optimization variables
    *
    * ### Solution algorithm ###
    *
    * The solver forms the Lagrange dual problem
    * \f[ \max_\lambda - \frac{1}{2} \lambda^TBA^{-1}B^T\lambda + (BA^{-1}c+b)^T \lambda \quad\text{s.t. } \lambda \le 0, \f]
    * which is bound constrained and can be solved by a projected gradient method.
    * As the primal problem is convex, strong duality holds, and the unique primal
    * minimizer \f$ x^* \f$ can be computed from a dual maximizer \f$ \lambda^* \f$ by
    * \f[ x^* = A^{-1}(B^T\lambda^*-c). \f]
    *
    *
    */
   template <int d, class Real=double>
   class QPSolver
   {
     using MatrixM = DynamicMatrix<Dune::FieldMatrix<Real,1,1>>;

   public:
     using VectorX = Dune::BlockVector<Dune::FieldVector<Real,d>>;
     using VectorB = Dune::BlockVector<Dune::FieldVector<Real,1>>;
     using MatrixA = DynamicMatrix<Dune::FieldMatrix<Real,d,d>>;
     using MatrixB = DynamicMatrix<Dune::FieldMatrix<Real,1,d>>;
     using Self = QPSolver<d,Real>;

     /**
      * \brief Constructs the solver, providing the matrices \f$ A \f$ and \f$ B \f$.
      *
      * The argument matrices are copied internally, and are not referenced later.
      *
      * Throws a NonpositiveMatrixException if \f$ A \f$ is not positive definite
      * and `convex` is QPConvexificationStrategy::DONOTHING.
      *
      * \todo write a move constructor
      */
     QPSolver(MatrixA const& A, MatrixB const& B, QPConvexificationStrategy convex=QPConvexificationStrategy::DONOTHING);

     QPSolver(Self const& qp);

     Self& operator=(Self const& qp);

     /**
      * \brief Solves the QP for the given right hand side vectors.
      *
      * Solves the QP by solving the dual QP with a projected gradient method
      * starting at the dually feasible point \f$ \lambda = 0\f$.
      *
      * \todo add option to provide an initial guess for \f$ \lambda \f$.
      *
      * \param c the linear term of the objective
      * \param b the constraints right hand side
      * \param tol the tolerance for termination
      *
      * The iteration is terminated as soon as the projected gradient \f$ s \f$ of the dual problem
      * is small, i.e. \f$ \|s\|_2 \le \mathsf{tol}\f$.
      *
      * \return a tuple \f$ (x,\lambda,iter) \f$ of solution vector, multiplier, and iteration count
      *
      * If the QP is infeasible, an InfeasibleProblemException is thrown.
      */
     std::tuple<VectorX,VectorB,int> solve(VectorX const& c, VectorB const& b, double tol) const;

     /**
      * \brief Provides access to the internal copy of \f$ B\f$.
      */
     MatrixB const& matrixB() const
     {
       return B;
     }
     
     /**
      * \brief Given info on whether dense or sparse arithmetic is used.
      * 
      */ 
      static constexpr QPStructure sparse = QPStructure::DENSE;

   private:
     MatrixA Ainv;
     MatrixB B;
     std::unique_ptr<QPBoundSolver<Real>> dual;
   };

  //--------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------

   /**
    * \ingroup qp
    * \brief An iterative solver for small instances of a particular class of quadratic programs.
    *
    * The solver addresses problems of the form
    * \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } Bx \le b, \f]
    * where \f$ A \f$ is positive semidefinite, and approximates the solution by minimizing
    * the Powell-Hestenes-Rockafellar augmented Lagrangian
    * \f[ \min_{x} \frac{1}{2} x^T Ax + c^T x
    *     + \frac{\gamma}{2}\left\|Bx-\Big(b-\frac{\lambda}{\gamma}\Big)\right\|_+^2 \f]
    * using a semismooth Newton method with exact linesearch.
    *
    * The penalty parameter \f$ \gamma \f$ affects both, solution quality and solution effort. For
    * \f$ \gamma\to\infty \f$, the original problem is better approximated, but the penalized problem's condition
    * gets worse.
    *
    * A good multiplier estimate \f$\lambda\f$ improves the solution accuracy significantly. Compared to the
    * standard augmented Lagrangian where a linear term \f$ \lambda^T (Bx-b) \f$ is added, the PHR version with
    * shiftet constraints has the advantage that in the inner bound constrained QPs that are to be solved the
    * strongly active constraints remain strongly active. In the standard formulation, they become weakly active
    * in the course of the AL iteration, which can induce frequent switching of constraints between active and
    * inactive.
    *
    * \tparam d the dimension of (vector valued) optimization variables
    * \tparam sparsity the type of linear algebra data structures to use, either DENSE or SPARSE
    * \tparam Real the real field type of matrix/vector entries.
    *              Explicit instantiations are provided for float and double.
    *
    */
   template <int d, QPStructure sparsity = QPStructure::DENSE, class Real=double>
   class QPPenalizedSolver
   {
     using AEntry = Dune::FieldMatrix<Real,d,d>;
     using BEntry = Dune::FieldMatrix<Real,1,d>;

   public:
     using VectorX = Dune::BlockVector<Dune::FieldVector<Real,d>>;
     using VectorB = Dune::BlockVector<Dune::FieldVector<Real,1>>;
     using MatrixA = std::conditional_t<sparsity==QPStructure::DENSE,
                                        DynamicMatrix<AEntry>,
                                        NumaBCRSMatrix<AEntry>>;
     using MatrixB = std::conditional_t<sparsity==QPStructure::DENSE,
                                        DynamicMatrix<BEntry>,
                                        NumaBCRSMatrix<BEntry>>;

     /**
      * \brief Constructs the solver, providing the matrices \f$ A \f$ and \f$ B \f$.
      *
      * The argument matrices are copied internally, and are not referenced later.
      *
      * \todo write a move constructor
      */
     QPPenalizedSolver(MatrixA const& A, MatrixB const& B);

     /**
      * \brief Constructs the solver, providing the matrices \f$ A \f$ and \f$ B \f$.
      *
      * The argument matrices are copied internally, and are not referenced later.
      *
      * \todo write a move constructor
      */
     QPPenalizedSolver(MatrixA const& A, MatrixB const& B, Real gamma);

     /**
      * \brief Sets the minimum number of projected gradient steps to perform.
      */
     QPPenalizedSolver<d,sparsity,Real>& setMinSteps(int i);

     /**
      * \brief Sets the maximum number of projected gradient steps to perform.
      */
     QPPenalizedSolver<d,sparsity,Real>& setMaxSteps(int i);

     /**
      * \brief Sets the penalty parameter \f$ \gamma \f$.
      *
      * The default on construction is \f$ \gamma = 10^3 \f$.
      */
     QPPenalizedSolver<d,sparsity,Real>& setPenalty(Real gamma);
     
     /**
      * \brief Returns the used penalty parameter \f$ \gamma \f$.
      */
     Real getPenalty() const { return gamma; }

     /**
      * \brief Solves the QP approximately for the given right hand side vectors.
      *
      * The starting point is \f$ x = 0 \f$.
      *
      * \param c the linear term of the objective
      * \param b the constraints right hand side
      * \param lambda the multiplier estimate. Due to the problem structure, \f$ \lambda \le 0 \f$ holds.
      * \param tol the tolerance for termination
      *
      * The iteration is terminated as soon as the projected gradient \f$ s \f$
      * is small, i.e. \f$ \|s\|_2 \le \mathsf{tol}\f$. Together with the solution vector,
      * a multiplier estimate based on the first order update  \f$ \lambda \gets \lambda - (Bx-s) \f$ is computed.
      *
      * \return a tuple \f$ (x,\delta\lambda,\mathrm{iter}) \f$ of solution vector, multiplier update, and iteration count.
      *
      * If the QP is infeasible, an InfeasibleProblemException is thrown.
      */
     std::tuple<VectorX,VectorB,int> solve(VectorX const& c, VectorB const& b, VectorB const& lambda, double tol) const;


     /**
      * \brief Solves the QP approximately for the given right hand side vectors.
      *
      * \param x the start iterate
      * \param c the linear term of the objective
      * \param b the constraints right hand side
      * \param lambda the multiplier estimate. Due to the problem structure, \f$ \lambda \ge 0 \f$ holds.
      * \param tol the tolerance for termination
      *
      * The iteration is terminated as soon as the projected gradient \f$ s \f$
      * is small, i.e. \f$ \|s\|_2 \le \mathsf{tol}\f$. Together with the solution vector,
      * a multiplier estimate based on the first order update  \f$ \lambda \gets \lambda - (Bx-s) \f$ is computed.
      *
      * \return a tuple \f$ (x,\delta\lambda,\mathrm{iter}) \f$ of solution vector, multiplier update, and iteration count
      *
      * If the QP is infeasible, an InfeasibleProblemException is thrown.
      */
     std::tuple<VectorX,VectorB,int> solve(VectorX const& c, VectorB const& b,
                                           VectorB lambda, double tol, VectorX const& x) const;

     std::tuple<VectorX,int> solve(VectorX const& c, VectorB const& b, double tol, VectorX const& x) const;
     std::tuple<VectorX,int> solve(VectorX const& c, VectorB const& b, double tol) const;

     /**
      * \brief Computes a least squares multiplier estimate.
      *
      * From the adjoint equation \f[ Ax+c - B^T \lambda = 0 \f] and a given approximate solution \f$ x \f$,
      * the multiplier can be estimated by minimizing the adjoint equation residual. This results
      * in \f[ BB^T \lambda = B(Ax+c). \f]
      * This estimate, projected back to \f$ \lambda \le 0 \f$, is computed here.
      *
      * The least squares multiplier estimate is well suited for obtaining an initial multiplier for
      * the augmented Lagrangian method, if only a good guess of the solution is available. Examples of
      * this are discretized inequality constraints, where the multiplier structure can change, rendering
      * previously computed multiplier useless (or requires at least a suitable prolongation).
      *
      * In contrast, the first order multiplier *update* works only if \f$ x \f$ is a minimizer of the
      * penalized problem with *known* multiplier.
      *
      * \todo Check if \f$ BA^{-1}B^T \lambda = BA + A^{-1}c \f$ is better or more appropriate, as it
      *       minimizes the residual in the problem's energy norm.
      *
      */
     VectorB multiplierEstimate(VectorX const& x, VectorX c) const;
     
    /**
      * \brief Given info on whether dense or sparse arithmetic is used
      *       (i.e. the choice of template parameter used in this class).
      * 
      */ 
      static constexpr QPStructure sparse = sparsity;

   private:
     Real gamma = 1e-3;
     MatrixA A;
     MatrixB B;

     int minSteps = 0;
     int maxSteps = 10000;

     std::vector<int> allCols;
   };

   // ---------------------------------------------------------------------------------------------
   // ---------------------------------------------------------------------------------------------
   
   /**
    * \ingroup qp
    * \brief An augmented Lagrangian solver for small instances of a particular class of quadratic programs.
    *
    * The solver addresses problems of the form
    * \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } Bx \le b, \f]
    * where \f$ A \f$ is positive semidefinite, and approximates the solution by minimizing
    * the Powell-Hestenes-Rockafellar augmented Lagrangian
    * \f[ \min_{x} \frac{1}{2} x^T Ax + c^T x
    *     + \frac{\gamma}{2}\left\|Bx-\Big(b-\frac{\lambda}{\gamma}\Big)\right\|_+^2 \f]
    * using a semismooth Newton method with exact linesearch and first-order multiplier updates.
    *
    * The penalty parameter \f$ \gamma \f$ affects both, solution quality and solution effort. For
    * \f$ \gamma\to\infty \f$, the augmented Lagrangian iteration converges faster, but the
    * penalized problem's condition gets worse. A reasonable default value is chosen based on
    * theoretical estimates of the convergence rate, aiming at a contraction factor of 0.01
    * (a very rough estimate, but should work as a default).
    *
    * A good multiplier estimate \f$\lambda\f$ improves the solution accuracy significantly. Compared to the
    * standard augmented Lagrangian where a linear term \f$ \lambda^T (Bx-b) \f$ is added, the PHR version with
    * shiftet constraints has the advantage that in the inner penalized QPs that are to be solved the
    * strongly active constraints remain strongly active. In the standard formulation, they become weakly active
    * in the course of the AL iteration, which can induce frequent switching of constraints between active and
    * inactive.
    *
    * \tparam d the dimension of (vector valued) optimization variables
    * \tparam sparsity the type of linear algebra data structures to use, either DENSE or SPARSE
    * \tparam Real the real field type of matrix/vector entries.
    *              Explicit instantiations are provided for float and double.
    *
    */
   template <int d, QPStructure sparsity = QPStructure::DENSE, class Real=double>
   class QPALSolver
   {
     using AEntry = Dune::FieldMatrix<Real,d,d>;
     using BEntry = Dune::FieldMatrix<Real,1,d>;

   public:
     using VectorX = Dune::BlockVector<Dune::FieldVector<Real,d>>;
     using VectorB = Dune::BlockVector<Dune::FieldVector<Real,1>>;
     using MatrixA = std::conditional_t<sparsity==QPStructure::DENSE,
                                        DynamicMatrix<AEntry>,
                                        NumaBCRSMatrix<AEntry>>;
     using MatrixB = std::conditional_t<sparsity==QPStructure::DENSE,
                                        DynamicMatrix<BEntry>,
                                        NumaBCRSMatrix<BEntry>>;

     /**
      * \brief Constructs the solver, providing the matrices \f$ A \f$ and \f$ B \f$.
      *
      * The argument matrices are copied internally, and are not referenced later.
      *
      * \todo write a move constructor
      */
     QPALSolver(MatrixA const& A, MatrixB const& B);

     /**
      * \brief Sets the minimum number of augmented Lagrangian steps to perform.
      * 
      * \param i the minimum number of steps, needs to be nonnegative
      * 
      * If i is larger than the maximum number of steps, the latter ist increased
      * to the same value.
      */
     QPALSolver<d,sparsity,Real>& setMinSteps(int i);

     /**
      * \brief Sets the maximum number of projected gradient steps to perform.
      * 
      * \param i the maximum number of steps, needs to be nonnegative
      * 
      * If i is less than the minimum number of steps, the latter ist decreased
      * to the same value.
      */
     QPALSolver<d,sparsity,Real>& setMaxSteps(int i);

     /**
      * \brief Sets the penalty parameter \f$ \gamma \f$.
      *
      * The default on construction is \f$ \gamma =  \f$.
      */
     QPALSolver<d,sparsity,Real>& setPenalty(Real gamma);

     Real getPenalty() const { return qp.getPenalty(); }

     /**
      * \brief Solves the QP approximately for the given right hand side vectors.
      *
      * The starting point is \f$ x = 0 \f$.
      *
      * \param c the linear term of the objective
      * \param b the constraints right hand side
      * \param tol the tolerance for termination
      *
      * The iteration is terminated as soon as the projected gradient \f$ s \f$
      * is small, i.e. \f$ \|s\|_2 \le \mathsf{tol}\f$. Together with the solution vector,
      * a multiplier estimate based on the first order update  \f$ \lambda \gets \lambda - (Bx-s) \f$ is computed.
      *
      * \return a tuple \f$ (x,\lambda,\mathrm{iter}) \f$ of solution vector, multiplier estimate, and iteration count.
      *
      * If the QP is infeasible, an InfeasibleProblemException is thrown.
      */
     std::tuple<VectorX,VectorB,int> solve(VectorX const& c, VectorB const& b, double tol) const;


     /**
      * \brief Solves the QP approximately for the given right hand side vectors.
      *
      * \param x the start iterate
      * \param c the linear term of the objective
      * \param b the constraints right hand side
      * \param tol the tolerance for termination
      *
      * The iteration is terminated as soon as the projected gradient \f$ s \f$
      * is small, i.e. \f$ \|s\|_2 \le \mathsf{tol}\f$. Together with the solution vector,
      * a multiplier estimate based on the first order update  \f$ \lambda \gets \lambda - (Bx-s) \f$ is computed.
      *
      * \return a tuple \f$ (x,\delta\lambda,\mathrm{iter}) \f$ of solution vector, multiplier update, and iteration count
      *
      * If the QP is infeasible, an InfeasibleProblemException is thrown.
      */
     std::tuple<VectorX,VectorB,int> solve(VectorX const& c, VectorB const& b, double tol, VectorX x) const;
     
     /**
      * \brief Given info on whether dense or sparse arithmetic is used
      *       (i.e. the choice of template parameter used in this class).
      * 
      */ 
      static constexpr QPStructure sparse = sparsity;

   private:
     int minSteps = 0;
     int maxSteps = 10000;
     
     QPPenalizedSolver<d,sparsity,Real> qp;
   };


}

#endif
