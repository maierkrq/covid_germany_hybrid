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

#ifndef QPMG_HH
#define QPMG_HH

#include <memory>
#include <tuple>

#include "dune/istl/bvector.hh"

#include "linalg/dynamicMatrix.hh"
#include "linalg/threadedMatrix.hh"
#include "mg/prolongation.hh"
#include "utilities/enums.hh"

namespace Kaskade
{

  /**
   * \ingroup qp
   * \brief An interface for gathering multigrid solver statistics.
   *
   * This interface provides also a trivial, i.e. do-nothing, implementation of
   * all virtual functions for entering values.
   */
  template <int d, class Real=double>
  class MGSolverStatistics
  {
  public:
    using VectorX = Dune::BlockVector<Dune::FieldVector<Real,d>>;
    using VectorB = Dune::BlockVector<Dune::FieldVector<Real,1>>;

    ~MGSolverStatistics() {}

    virtual void enterPreSmoothingCorrection(int level, VectorX const& dx) {}
    virtual void enterCoarseGridCorrection(int level, VectorX const& dx) {}
    virtual void enterPostSmoothingCorrection(int level, VectorX const& dx) {}

    /**
     * \brief At the end of a V-cycle, specify the iteration number, the correction, and the resulting energy
     *        at the new iterate.
     *        
     * \param activeSetChanges the number of inequality constraints that changed their activity state in this iteration
     */
    virtual void enterVCycleCorrection(int iter, VectorX const& dx, Real energy, 
                                      VectorX const& res1, VectorX const& res2,
                                      VectorB const& cons, int activeSetChanges, 
                                      std::vector<VectorX> &correctionStack,
                                      std::vector<std::tuple<double,double,double>> &coarseRes) {}
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------



   /**
    * \ingroup qp
    * \brief A base class for multigrid solver for QPs. 
    *
    * This base class implements a multigrid solver for problems of the form
    * \f[ \min_x \frac{1}{2} x^T Ax + c^T x \quad \text{s.t. } Bx \le b, \f]
    * where \f$ A \f$ is positive semidefinite defined on a hierarchy of grids,
    * i.e. for a given prolongation stack defining the (geometric) relation between
    * a sequence of grids.
    * 
    * To the QP on fine grids we apply a block Gauss-Seidel method, which treats 
    * DOFs coupled by the constraints \f[ Bx \le b \f] jointly. These local QPs are 
    * then solved by the smoother on the corresponsing level.
    * On the corse grid the system is solver either by a direct solver (dense or sparse)
    * or by a few smoothing steps.
    * 
    *
    * \tparam d               the dimension of (vector valued) optimization variables
    * \tparam Prolongation    defines the data type of the prolongation matrices defining the hierarchy
    * \tparam Smoother        the smoother for the local QPs on all grids but the coarse one
    * \tparam CoarseSmoother  the smoother for the coarse grid
    *
    */
  template <int d, class Prolongation, class Smoother, class CoarseSmoother=Smoother, class Real=double>
  class QPMultiGridBase
  {
    using EntryA = Dune::FieldMatrix<Real,d,d>;
  public:
    using VectorX = Dune::BlockVector<Dune::FieldVector<Real,d>>;
    using VectorB = Dune::BlockVector<Dune::FieldVector<Real,1>>;
    using MatrixA = NumaBCRSMatrix<EntryA>;
    using MatrixB = NumaBCRSMatrix<Dune::FieldMatrix<Real,1,d>>;
    using Self = QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>;
    
    // check if coarse smoother is of equal type as rest of smoothers
    // if so, SmootherType = Smoother, else use std::variant<Smoother,CoarseSmoother>
    static constexpr bool smoothersEqual = std::is_same<Smoother,CoarseSmoother>::value;
    using SmootherType = std::conditional_t<smoothersEqual,
                                            Smoother,
                                            std::variant<Smoother,CoarseSmoother>>;

    /**
     * \brief Constructs the solver, providing the matrices \f$ A \f$ and \f$ B \f$.
     *
     * \param A the Hessian. The matrix is referenced internally and has to exist during the lifetime of the smoother.
     * \param B the constraints. The matrix is referenced internally and has to exist during the lifetime of the smoother.
     *
     * \param blocks
     *    - block=true: treat blocks of dof coupled by constraints jointly (default)
     *    - block=false: treat all dofs individually (only for testing and comparison purposes)
     *    \warning This is intended for testing purposes only and may be removed in future.
     * 
     * \param directOnCoarse
     *    - directOnCoarse=true   use direct solver for problem on coarse grid
     *    - directOnCoarse=false  use a few iterative (smoothing) steps on coarse grid (block Jacobi/Gauss Seidel)
     */
    QPMultiGridBase(MatrixA A, MatrixB const& B, std::vector<Prolongation>&& prolongations,
                    Real smootherRegularization=0, bool blocks=true, bool directOnCoarse=true);

    /**
     * \brief Destructor.
     */
    virtual ~QPMultiGridBase();

    /**
     * \brief Computes a single V-cycle correction for approximately solving the QP for the given right hand side vectors.
     *
     * \param x initial guess for the solution
     * \param c linear objective term
     * \param b constant constraint term
     * \return the correction dx
     */
    VectorX step(VectorX const& x, VectorX c, VectorB b) const;

    /**
     * \brief Approximately solves the QP for the given right hand side vectors up to a given tolerance.
     *
     * \param x initial guess for the solution
     * \param c linear objective term
     * \param b constant constraint term
     * \param tol tolerance
     * \return a tuple \f$ (x,\mathrm{iter}) \f$ of solution vector, multiplier update, and iteration count
     */
    std::tuple<VectorX,int> solve(VectorX x, VectorX const& c, VectorB const& b, double tol, int vcycles) const;


    /**
     * \brief Defines the grid level to be used as coarse grid.
     *
     * \param coarselevel the coarse grid level, between 0 and the max grid level
     *
     * If coarselevel is outside its allowed range, it is projected onto that range.
     */
    Self& setCoarseLevel(int coarselevel);

    Self& setSmoothings(int pre, int post);

    /**
     * \brief Sets the number of coarse corrections to compute and apply.
     *
     * \param n number of coarse corrections
     *
     * - n=0: just the smoother on the fine grid is performed, no multigrid at all (for testing purposes only).
     * - n=1: V-cycle (default)
     * - n=2: W-cycle
     */
    Self& setCoarseCorrections(int n);

    /**
     * \brief Defines the smoothing strategy to use for unconstrained degrees of freedom.
     *
     * In contact problems, those are usually the interior nodes, and the ones without
     * potential contact partner. The smoother can operate in parallel (Jacobi) or sequentially
     * (Gauss-Seidel). In the parallel way, an exact linesearch is performed after the smoothing
     * in order to guarantee descent.
     */
    Self& setBulkMode(ParallelMode m);

    /**
     * \brief Sets the logging facility for reporting solver statistics.
     *
     * \param logger a pointer to the logger object. This can be a nullpointer, in which case logging is stopped.
     *
     * Ownership of the logger object remains with the caller. The logger object
     * needs to exist for the lifetime of the QPMultiGrid class or until it is replaced
     * by a different logger (or none).
     */
    Self& setLogger(MGSolverStatistics<d,Real>* logger);

  protected:
    /**
     * \brief Compute the problem energy of current iterate for logging purposes.
     */
    virtual double computeEnergy(MatrixA const& A, MatrixB const& B,
                                 VectorX const& c, VectorB const& b, VectorX const& x) const = 0;
     /**
     * \brief Compute the problem gradient of current iterate for adaptive termination criterion.
     */                 
    virtual std::tuple<VectorX,VectorB,int> gradient(MatrixA const& A, MatrixB const& B,
                                                     VectorX const& c, VectorB const& b, VectorX const& x) const = 0;
                                                     
    /**
     * \brief Compute the problem residual of current iterate for adaptive termination criterion.
     */                 
    virtual std::tuple<VectorX,VectorX,VectorB,int> residual(MatrixA const& A, MatrixB const& B,
                                                             VectorX const& c, VectorB const& b, VectorX const& x) const = 0;
    
    /**
     * \brief Compute the optimal step size in direction of dx.
     */
    virtual double qpLinesearch(MatrixA const& A, MatrixB const& B, VectorX const& c,
                                VectorB const& b, VectorX const& dx) const = 0;

    
    // vector of smoothers, one for each level
    // note: SmootherType might be either of Type
    // Smoother or std::variant<Smoother,CoarseSmoother>
   std::vector<SmootherType> smoothers;

    /**
     * \brief Provides a reference to the fine grid constraints.
     */
    MatrixB const& B() const;

  private:
    /**
     * \brief Solves a (possibly penalized) QP problem by a multigrid V-cycle
     *
     * \param b upper bound for \f$ s \f$, has to be nonnegative (\f$ b\ge 0 \f$)
     *
     * \return the solution \f$ (x,s) \f$
     */
    VectorX mgRecursion(int level, int lpre, int lpost, VectorX c, VectorB b) const;

    MultiGridStack<Prolongation,EntryA,size_t> mgStackA;
    std::vector<MatrixB> mgStackB;
    mutable std::vector<VectorX> correctionStack;                                 // vector with coarse grid corrections on each level (interpolated to fine grid)
    mutable std::vector<std::tuple<double,double,double>> coarseResidual;         // vector with norm of residual on each level (after computing correction)

    mutable Real gamma;
    int coarseLevel = 0;
    int pre = 3;
    int post = 3;
    int mid = 1;

    std::unique_ptr<MGSolverStatistics<d,Real>> dummyLogger; // this dummy will be used if no other logger is given
    MGSolverStatistics<d,Real>* logger;
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  // forward declaration
  template <int d, class Real, QPStructure sparsity=QPStructure::DENSE>
  class QPPenaltySmoother;



  /**
   * \ingroup qp
   * \brief A multigrid solver for convex QPs.
   *
   * This solves QPs of the form
   * \f[ \min \frac{1}{2} x^T A x + c^T x \quad \text{s.t.} \quad Bx\le b. \f]
   *
   * On the finest level, this works with the augmented Lagrangian formulation
   * \f[ \min \frac{1}{2} x^T A x + c^T x + \frac{\gamma}{2} \|(Bx-(b-\lambda/\gamma))_+\|^2  \f]
   * for some multiplier \f$ \lambda \ge 0 \f$, and an appropriate multiplier update.
   *
   * The multilevel minimization works with level-dependent penalty factors. On each level, the
   * following computations (\f$\text{MGM}(A,B,x,c,b,\gamma)\f$ are performed:
   * - pre-smoothing: apply a block Jacobi smoother to the augmented Lagrangian.
   * - multiplier update (optional: is interleaving multigrid with multiplier updates a good idea?)
   * - coarse grid correction: compute residuals  \f$ c_r \f$ and \f$ b_r \f$ and compute a correction
   *   \f[ (\delta x, \delta \lambda) = \text{MGM}(P^TAP,BP,0,P^Tc_r,b_r,\gamma/4) \f]
   * - post-smoothing: apply correction to primal variable and multiplier and perform
   *   block Jacobi
   * - multiplier update
   *
   * Rationale: The level-dependence of the penalty factor is based on the following idea.
   * Enforcing the constraint by the penalty should allow sufficient slack such as not to
   * impede convergence of Gauß-Seidel. On coarser levels, GS is supposed to take larger steps
   * (the whole reason for multigrid) approximately by a factor of two for each grid level,
   * and so the penalty should increase the permitted steps by a factor of two as well.
   * This means decreasing the penalty factor by four.
   *
   * The restriction of multipliers is not necessary for asymptotically optimal complexity,
   * as the number of constraints is \f$ \mathcal{O}(N^{1-1/d}) \f$ and therefore the
   * computational work of a single V-cycle is \f$ \mathcal{O}(N + N^{1-1/d}\log N) = \mathcal{O}(N) \f$.
   * Nevertheless, restricting the constraints in an appropriate way (open how to do this) might
   * be beneficial for performance.
   * 
   * \tparam d               the dimension of (vector valued) optimization variables
   * \tparam Prolongation    defines the data type of the prolongation matrices defining the hierarchy
   * \tparam Real            floating point format
   * \tparam Smoother        the smoother for the local QPs on all grids but the coarse one
   * \tparam CoarseSmoother  the smoother for the coarse grid
   * 
   */
  template <int d, class Prolongation, class Real=double, class Smoother=QPPenaltySmoother<d,Real>, class CoarseSmoother=QPPenaltySmoother<d,Real>>
  class QPMultiGrid
  : public QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>
  {
    using EntryA = Dune::FieldMatrix<Real,d,d>;
    using Base = QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>;

  public:
    using VectorX = Dune::BlockVector<Dune::FieldVector<Real,d>>;
    using VectorB = Dune::BlockVector<Dune::FieldVector<Real,1>>;
    using MatrixA = NumaBCRSMatrix<EntryA>;
    using MatrixB = NumaBCRSMatrix<Dune::FieldMatrix<Real,1,d>>;
    using Self = QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>;

    /**
     * \brief Constructs the solver, providing the matrices \f$ A \f$ and \f$ B \f$.
     *
     * \param A the Hessian. The matrix is copied internally and need not exist after solver construction.
     * \param B the constraints. The matrix is referenced internally and has to exist during the lifetime of the solver.
     *
     * \param blocks
     *    - block=true: treat blocks of dof coupled by constraints jointly (default)
     *    - block=false: treat all dofs individually (only for testing and comparison purposes)
     *    \warning This is intended for testing purposes only and may be removed in future.
     * 
     * \param directOnCoarse
     *    - directOnCoarse=true   use direct solver for problem on coarse grid
     *    - directOnCoarse=false  use a few iterative (smoothing) steps on coarse grid (block Jacobi/Gauss Seidel)
     */
    QPMultiGrid(MatrixA const& A, MatrixB const& B, std::vector<Prolongation>&& prolongations,
                 Real smootherRegularization=0, bool blocks=true, bool directOnCoarse=true);

    /**
     * \brief Destructor.
     */
    ~QPMultiGrid();

    /**
     * \brief Approximately solves the QP for the given right hand side vectors up to a given tolerance.
     *
     * \param x initial guess for the solution
     * \param c linear objective term
     * \param b constant constraint term
     * \param lambda approximate multiplier
     * \param tol tolerance
     * \return a tuple \f$ (x,\delta\lambda,\mathrm{iter}) \f$ of solution vector, multiplier update, and iteration count
     */
    std::tuple<VectorX,VectorB,int> solve(VectorX x, VectorX const& c, VectorB const& b, VectorB const& lambda, double tol, int vcycles) const;

//     /**
//      * \brief A least-squares multiplier update.
//      *
//      * This can be used if a good initial guess for \f$ x \f$ is available, but not for \f$ \lambda \f$.
//      * A common situation where this is typically the case is finite strain contact mechanics, where the
//      * constraints and hence the multiplier changes with deformation, after having solved a QP linearization.
//      */
//     void leastSquaresMultiplierUpdate(VectorX const& c, VectorB const& b, VectorX const& x, VectorB& lambda) const;
//

    /**
     * \brief A first order update for the Lagrange multiplier
     * 
     * Compute first order multiplier update \lambda <- (\lambda + \gamma*(B*x-b))_+,
     * or equivalently \lambda <- \lambda + max(-\lambda, \gamma*(B*x-b)).
     * This makes only sense, however, if we have actually minimized - otherwise
     * there is no new information in the constraint violation. Thus the accuracy
     * of minimization should be sufficient.
     *
     * \param x primal variable 
     * \param b constant constraint term
     * \param lambda approximate multiplier
     * \return multiplier update \f$ \delta\lambda \f$ 
     */
    VectorB firstOrderUpdate(VectorX const& x, VectorB const& b, VectorB const& lambda) const;


    /**
     * \brief Sets the penalty factor \f$ \gamma \f$.
     *
     * The default on construction is \f$ \gamma = 10^3 \f$.
     *
     * \warning set methods in the base class return a reference to the base class,
     *          which impairs chaining of setters. Call the derived class setters first, e.g.
     *          \code
     *          mg.setPenalty(1e3).setCoarseCorrections(2);
     *          \endcode
     * 
     * \param gamma the penalty parameter
     * \param levelDependent if set true, we use level dependent penaltys,
     *                       i.e., gamma on finest level and gamma/4 recursively
     *                       on coarser levels
     */
    Self& setPenalty(Real gamma, bool levelDependent=false);

    Self& setConstraintsMode(ParallelMode m);
    Self& setGlobalMode(ParallelMode m);
    
    /**
     * \brief Returns the current penalty factor \f$ \gamma \f$.
     * 
     */
    double getPenalty();

  protected:
    /**
     * \brief Computes the energy, given by the augmented Lagrangian function
     *    \f[ \frac{1}{2} x^T A x + c^T x + \frac{\gamma}{2} \|(Bx-b)_+\|^2  \f]
     * 
     */
    virtual double computeEnergy(MatrixA const& A, MatrixB const& B,
                                 VectorX const& c, VectorB const& b, VectorX const& x) const;
    /**
     * \brief Computes the gradient of the augmented Lagrangian function
     *    \f[ g = Ax + c + \gamma B^T(Bx-b)_+  \f]
     *  
     * \return a tuple  (\f$ g \f$,con,nActive)  with the gradient, the constraint
     *         violation  con = \f$ (Bx-b)_+ \f$ and the number of active constraints nActive
     */                      
    virtual std::tuple<VectorX,VectorB,int> gradient(MatrixA const& A, MatrixB const& B,
                                                     VectorX const& c, VectorB const& b, VectorX const& x) const;
                                                     
   /**
     * \brief Computes the residual of the augmented Lagrangian function
     *    \f[ res = Ax + c + \gamma B^T(Bx-b)_+  \f]
     *    as a sum of both addends
     *  
     * \return a tuple  (\f$ res1, res2 \f$,con,nActive)  with the summands
     *         \f res1 = Ax+c, res2 = \gamma* B^T(Bx-b)_+ \f
     *         and the constraint violation  con = \f$ (Bx-b)_+ \f$ and the number of active constraints nActive
     */                      
    virtual std::tuple<VectorX,VectorX,VectorB,int> residual(MatrixA const& A, MatrixB const& B,
                                                             VectorX const& c, VectorB const& b, VectorX const& x) const;
    /**
     * 
     * 
     * 
     */
    virtual double qpLinesearch(MatrixA const& A, MatrixB const& B, VectorX const& c,
                                VectorB const& b, VectorX const& dx) const;


  private:
    mutable Real gamma;
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  // forward declaration
  template <int d, class Real, QPStructure sparsity=QPStructure::DENSE>
  class QPSmoother;

  template <int d, class Prolongation, class Real=double, class Smoother=QPSmoother<d,Real>, class CoarseSmoother=QPSmoother<d,Real>>
  class QPMultiGridStrict
  : public QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>
  {
    using Base = QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>;

  public:
    using Self = QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>;
    using typename Base::MatrixA;
    using typename Base::MatrixB;
    using typename Base::VectorX;
    using typename Base::VectorB;

    /**
     * \brief Constructs the solver, providing the matrices \f$ A \f$ and \f$ B \f$.
     *
     * \param A the Hessian. The matrix is copied internally and need not exist after solver construction.
     * \param B the constraints. The matrix is referenced internally and has to exist during the lifetime of the solver.
     *
     * \param blocks
     *    - block=true: treat blocks of dof coupled by constraints jointly (default)
     *    - block=false: treat all dofs individually (only for testing and comparison purposes)
     *    \warning This is intended for testing purposes only and may be removed in future.
     * 
     * \param directOnCoarse
     *    - directOnCoarse=true   use direct solver for problem on coarse grid
     *    - directOnCoarse=false  use a few iterative (smoothing) steps on coarse grid (block Jacobi/Gauss Seidel)
     */
    QPMultiGridStrict(MatrixA const& A, MatrixB const& B, std::vector<Prolongation>&& prolongations,
                      Real smootherRegularization=0, bool blocks=true, bool directOnCoarse=true);

    /**
     * \brief Destructor.
     */
    ~QPMultiGridStrict();
    
    
    /**
     * \brief Approximately solves the QP for the given right hand side vectors up to a given tolerance.
     *
     * \param x initial guess for the solution
     * \param c linear objective term
     * \param b constant constraint term
     * \param lambda approximate multiplier
     * \param tol tolerance
     * \return a tuple \f$ (x,\delta\lambda,\mathrm{iter}) \f$ of solution vector, multiplier update, and iteration count
     */
    std::tuple<VectorX,int> solve(VectorX x, VectorX const& c, VectorB const& b, double tol, int vcycles) const;


  protected:
    virtual double computeEnergy(MatrixA const& A, MatrixB const& B,
                                 VectorX const& c, VectorB const& b, VectorX const& x) const;
                                 
   /**
     * \brief Computes the gradient of the (non-augmented) Lagrangian
     *    \f[ g = Ax + c \f]
     *  
     * \return a tuple  (\f$ g \f$,con,nActive)  with the gradient, the constraint
     *         violation  con = \f$ (Bx-b)_+ \f$ and the number of active constraints nActive
     */                      
    virtual std::tuple<VectorX,VectorB,int> gradient(MatrixA const& A, MatrixB const& B,
                                                     VectorX const& c, VectorB const& b, VectorX const& x) const;
                                                     
                                                     
                                                        /**
   * \brief Computes the residual of the augmented Lagrangian function
     *    \f[ res = Ax + c + \gamma B^T(Bx-b)_+  \f]
     *    as a sum of both addends
     *  
     * \return a tuple  (\f$ res1, res2 \f$,con,nActive)  with the summands
     *         \f res1 = Ax+c, res2 = \gamma* B^T(Bx-b)_+ \f
     *         and the constraint violation  con = \f$ (Bx-b)_+ \f$ and the number of active constraints nActive
     */                      
    virtual std::tuple<VectorX,VectorX,VectorB,int> residual(MatrixA const& A, MatrixB const& B,
                                                             VectorX const& c, VectorB const& b, VectorX const& x) const;

    virtual double qpLinesearch(MatrixA const& A, MatrixB const& B, VectorX const& c,
                                VectorB const& b, VectorX const& dx) const;
  };
}

#endif
