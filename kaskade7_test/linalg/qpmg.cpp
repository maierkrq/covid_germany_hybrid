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

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <iterator>
#include <mutex>
#include <math.h>

#include "dune/grid/config.h"

#include "fem/firstless.hh"
#include "fem/fixdune.hh"
#include "io/matlab.hh"
#include "linalg/apcg.hh"
#include "linalg/conjugation.hh"
#include "linalg/dynamicMatrixOps.hh"
#include "linalg/matrixProduct.hh"
#include "linalg/qp.hh"
#include "linalg/qpLinesearch.hh"
#include "linalg/qpmg.hh"
#include "linalg/symmetricOperators.hh"
#include "linalg/direct.hh"
#include "linalg/threadedMatrix.hh"
#include "linalg/triplet.hh"
#include "linalg/crsutil.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/power.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"

namespace Kaskade
{
  // ----------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------

  namespace
  {

    // From the constraints matrix, extract groups of primal variables and associated constraints. It returns a triple
    // (p,c,u) such that
    // p[i] is a vector of primal variables all of which share at least one constraint (sorted ascendingly)
    // c[i] is a vector of all the constraints that are affected by primal variables in p[i] (sorted ascendingly)
    // u    is a vector of all structurally unconstrained primal variables  (sorted ascendingly)
    template <class Matrix>
    std::tuple<std::vector<std::vector<size_t>>,
               std::vector<std::vector<size_t>>,
               std::vector<size_t>> getDofSets(Matrix const& B, bool blocks)
    {
      using DofSet = std::vector<size_t>;
      std::vector<DofSet> dofSets(B.N());

      if (blocks)
      {
        // We shall couple dofs sharing a joint constraint. For each constraint, get the dofs that affect it
        // (they are encountered in increasing order in the row of B).
        for (size_t j=0; j<B.N(); ++j)
        {
          auto const& row = B[j];
          dofSets[j].reserve(row.size());
          for (auto ci=row.begin(); ci!=row.end(); ++ci)
            dofSets[j].push_back(ci.index());
        }
      }
      else
      {
        // We shall not couple dofs by constraints. Thus, any constrained dof is treated as
        // a dof set of its own.
        for (size_t j=0; j<B.N(); ++j)
        {
          auto const& row = B[j];
          for (auto ci=row.begin(); ci!=row.end(); ++ci)
            dofSets.push_back(DofSet{ci.index()});
        }
      }

      // These dof sets may be redundant (if rows of B have the same sparsity structure), thus we
      // remove duplicates. For that we first sort the dof sets lexicographically decreasing. This
      // implies that a precedes b if a is a superset of b.
      std::sort(begin(dofSets),end(dofSets),std::greater<DofSet>());

      // Now remove duplicates of dof sets. For duplicate removal, sets are considered equal if one is
      // a superset of the other. As supersets have been sorted to the beginning, and unique keeps the
      // first of equal elements, the largest set is kept.
      auto last = std::unique(begin(dofSets),end(dofSets),[](DofSet const& a, DofSet const& b)
      {
        auto n = std::min(a.size(),b.size());
        return equal(begin(a),begin(a)+n,begin(b));
      });
      dofSets.erase(last,end(dofSets));

      // For each dof set, consider all the affected constraints. For that, we need the transpose
      // pattern of B. TODO: rewrite NumaCRSPattern to provide useful access to the sparsity pattern
      // Note that Bt[i] is a sorted vector of constraint indices (ascending)
      std::vector<DofSet> Bt(B.M());
      for (size_t j=0; j<B.N(); ++j)
        for (auto ci=B[j].begin(); ci!=B[j].end(); ++ci)
          Bt[ci.index()].push_back(j);

      // For each dof set, compute the union of all affected constraints of all the dofs.
      // The resulting set of constraint indices is sorted ascendingly.
      std::vector<DofSet> conSets(dofSets.size());
      for (size_t k=0; k<dofSets.size(); ++k)
      {
        // compute the union by concatenating the constraints, sorting, and removing duplicates
        auto& con = conSets[k];
        for (size_t i: dofSets[k])
          con.insert(end(con),begin(Bt[i]),end(Bt[i]));
        std::sort(begin(con),end(con));
        con.erase(std::unique(begin(con),end(con)),end(con));
      }

      // Now gather all structurally unconstrained degrees of freedom.
      DofSet uncon;
      for (size_t i=0; i<B.M(); ++i)
        if (Bt[i].empty())
          uncon.push_back(i);

      return make_tuple(dofSets,conSets,uncon);
    }
  }


  using ProlongationMatrix = NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  namespace
  {

    // --------------------------------------------------------------------------------------------

   template <class Prolongation>
    Prolongation const& asMatrix(Prolongation const& p)
    {
      return p;
    }

    ProlongationMatrix asMatrix(MGProlongation const& p)
    {
      return p.template asMatrix<double>();
    }

  }

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <int d, class QPSolver, class Real>
  class QPSmootherBase
  {
  public:
    using VectorX = Dune::BlockVector<Dune::FieldVector<Real,d>>;
    using VectorB = Dune::BlockVector<Dune::FieldVector<Real,1>>;
    using EntryA  = Dune::FieldMatrix<Real,d,d>;
    using EntryB  = Dune::FieldMatrix<Real,1,d>;
    using MatrixA = NumaBCRSMatrix<EntryA>;
    using MatrixB = NumaBCRSMatrix<EntryB>;  
    // define type of local matrices in dependence of arithmetic 
    // used in unerlying solver (DENSE or SPARSE)
    using MatrixAloc = std::conditional_t<QPSolver::sparse==QPStructure::DENSE,
                                            DynamicMatrix<EntryA>,
                                            NumaBCRSMatrix<EntryA>>;
    using MatrixBloc = std::conditional_t<QPSolver::sparse==QPStructure::DENSE,
                                            DynamicMatrix<EntryB>,
                                            NumaBCRSMatrix<EntryB>>;
    using Self = QPSmootherBase<d,QPSolver,Real>;

    /**
     * \brief Constructor.
     *
     * The argument matrices A, B, and BtB are referenced internally and need to exist during the life time of the
     * smoother.
     */
    QPSmootherBase(MatrixA const& A_, MatrixB const& B_, bool direct=false, Real regularization=0, bool blocks=true)
    : A(&A_), B(&B_)
    {
      if (direct)                                         // We shall use a direct solver for the complete problem.
      {
        dofSets.push_back(DofSet(A->N()));                // That means there is just one dof set, containing
        std::iota(begin(dofSets[0]),end(dofSets[0]),0);   // all the dofs from 0 on.
        conSets.push_back(DofSet(B->N()));                // And there's one constraint set, containing all
        std::iota(begin(conSets[0]),end(conSets[0]),0);   // the constraints.
      }                                                   // uncon remains empty, as there are no "free" dofs
      else
        std::tie(dofSets,conSets,uncon) = getDofSets(*B,blocks);


      // Here we do the splitting of A and B into sub-blocks.
      for (int i=0; i<dofSets.size(); ++i)
      {
        MatrixAloc Aloc;
        if constexpr(QPSolver::sparse == QPStructure::DENSE)
          Aloc = full(*A,dofSets[i],dofSets[i]);
        else
          Aloc = submatrix<MatrixA>(*A,dofSets[i],dofSets[i]);
        
        // The matrix A may be *very* ill-conditioned, e.g. if there are low-energy modes in A like
        // rigid body motions in solid mechanics, and these low-energy modes barely affect constraints,
        // e.g., a rotating circular disc of stiff material inside a surrounding circle. In these cases,
        // rounding errors may lead to a semidefinite or indefinite A, and an ubounded QP. Then the QP
        // solver will fail. We therefore add positive values to the diagonal entries such that the
        // condition number will be bounded well below 1/eps.
        // WARNING: The drawback is that large penalties (e.g. for directional Dirichlet boundary conditions)
        // can prevent sufficient progress in non-penalized directions.
        Real maxDiag = 0;
        for (int k=0; k<Aloc.N(); ++k)
          maxDiag = std::max(maxDiag,Aloc[k][k].frobenius_norm());
        for (int k=0; k<Aloc.N(); ++k)
          Aloc[k][k] += regularization*unitMatrix<Real,d>();

        MatrixBloc Bloc;
        if constexpr(QPSolver::sparse == QPStructure::DENSE)
          Bloc = full(*B,conSets[i],dofSets[i]);
        else
          Bloc = submatrix<MatrixB>(*B,conSets[i],dofSets[i]);
          
        localQPs.push_back(QPSolver(Aloc,Bloc));
      }

      // For the unconstrained variables, we perform a plain Jacobi step. For that we need to
      // solve with the diagonal. Since accessing the diagonal means we have to do a binary
      // search in each row, and then invert the diagonal entry (which may be a small matrix itself),
      // this is expensive. We extract and invert the diagonal once beforehand.
      unconDiagAinv.resize(uncon.size());
      for (size_t i=0; i<uncon.size(); ++i)
      {
        unconDiagAinv[i] = (*A)[uncon[i]][uncon[i]];
        unconDiagAinv[i].invert();
      }
    }

    virtual ~QPSmootherBase()
    {}

    Self& setBulkMode(ParallelMode m)
    {
      if (m==ParallelMode::BLOCKPARALLEL)
        m = ParallelMode::PARALLEL;
      bulkMode = m;
      return *this;
    }

    Self& setConstraintsMode(ParallelMode m)
    {
      constraintsMode = m;
      return *this;
    }

    Self& setGlobalMode(ParallelMode m)
    {
      // TODO: implement block parallel mode
      if (m==ParallelMode::BLOCKPARALLEL)
        m = ParallelMode::PARALLEL;
      globalMode = m;
      return *this;
    }

    Self& setQPtol(Real tol)
    {
      qpTol = std::max((Real)0,tol);
      return *this;
    }



    /**
     * \brief Smoother, solves the QP approximately by solving many local QP problems in a Gauss-Seidel manner.
     *
     * This solves \f[ x^* = \mathrm{arg}\min \frac{1}{2} x^T A x + c^T x + \frac{gamma}{2} \| Bx - b\|_+^2 \f]
     * approximately, with initial guess zero, and returns the correction \f$ \delta x \f$.
     *
     * \param[out] dx the correction \f$ \delta x \f$. It is resized to the correct length.
     */
    void applyCorrection(VectorX c, VectorB b, VectorX& dx) const
    {
      ScopedTimingSection timerSection("QPMultiGrid smoother");
      auto& timer = Timings::instance();

      // We want to minimize the penalized augmented problem
      // 1/2 dx^T Adx + c^T dx + gamma/2 |Bdx-b|_+^2.

      timer.start("init");
      dx.resize(c.N());
      dx = 0;
      timer.stop("init");

      // Perform standard Jacobi on the unconstrained primal variables. Note that on this
      // subspace, the constraints and hence B play no role at all.
      // TODO: This is memory bound. For parallelization, consider to adapt JacobiSmoother on Numa machines.
      timer.start("bulk smoother");
      for (int i=0; i<uncon.size(); ++i)
      {
        int k = uncon[i];
        dx[k] = - (unconDiagAinv[i]*c[k]);

        // This is Gauss-Seidel. Actually it improves the contraction rate, but apparently
        // the (significant) increase in computational work leads to increased run times nevertheless.
        if (bulkMode==ParallelMode::SEQUENTIAL)
          for (auto ci=(*A)[k].begin(); ci!=(*A)[k].end(); ++ci)
            ci->umtv(dx[k],c[ci.index()]);
      }
      timer.stop("bulk smoother");

      // If we want to apply bulk and constrained dof sets sequentially, we need to
      // update the right hand side. If we have done the bulk sequentially, this has
      // already been done. Otherwise, we need additionally to make sure that we have
      // energy decrease in the (parallel) bulk step.
      if (globalMode==ParallelMode::SEQUENTIAL && bulkMode==ParallelMode::PARALLEL)
      {
        Real t = qpLinesearch(*A,*B,0.0,c,b,dx);// ensure energy decrease (TODO: we could use simpler and more efficient linesearch here)
                                                // Note that as in the bulk we have only unconstrained dofs, the value
                                                // of gamma is irrelevant. Set to 0 here.
        dx *= t;
        A->usmv(1.0,dx,c);                      // modify rhs for boundary smoother step (sequential, after free dof smoother)
      }

      // Perform overlapping block Jacobi/Gauss-Seidel QP on the constrained variables
      {
        ScopedTimingSection tim("QP smoother");

        // The number of parallel tasks to run depends on the requested mode of operation for the constrained dofs.
        // Sequential execution is enforced by taking a single task to execute all local PQPs.
        // TODO: use a coloring of the blocks and perform non-overlapping ones in parallel.
        int nTasks = constraintsMode==ParallelMode::SEQUENTIAL? 1: 1024;
        std::mutex rhsMutex;    // for preventing concurrent access to solution and right hand side

        parallelFor(0,dofSets.size(),[&](int const i)
        {
          // We extract the block subvectors of c and b
          int n = dofSets[i].size();
          int m = conSets[i].size();

          VectorX cloc(n);
          VectorB bloc(m);

          for (int j=0; j<dofSets[i].size(); ++j)
              cloc[j] = c[dofSets[i][j]];

          for (int j=0; j<conSets[i].size(); ++j)
            bloc[j] = b[conSets[i][j]];

          // Solve the QP for correction dx. As we go for corrections, we start at 0
          // and thus do not provide an initial guess.
          auto dxloc = std::get<0>(localQPs[i].solve(cloc,bloc,qpTol));

          // Update the solution. Take care of proper locking (as we run potentially in parallel)
          // and update the solution vector with a damping factor of 0.5 (see below).
          std::lock_guard<std::mutex> rhsLock(rhsMutex);
          for (int j=0; j<dofSets[i].size(); ++j)
            dx[dofSets[i][j]] += dxloc[j];

          // If we work sequentially or at least block-parallel, we must update the right hand
          // side and bound for the next dof sets to work on.
          if (constraintsMode != ParallelMode::PARALLEL)
          {
            // Update right hand side and constraints c <- c+A*dx, b <- b-B*dx
            for (int j=0; j<dofSets[i].size(); ++j)
            {
              int dof = dofSets[i][j];
              for (auto ci=(*A)[dof].begin(); ci!=(*A)[dof].end(); ++ci)
                ci->umtv(dxloc[j],c[ci.index()]);
            }
            for (int k=0; k<conSets[i].size(); ++k)
            {
              int con = conSets[i][k];
              for (auto ci=(*B)[con].begin(); ci!=(*B)[con].end(); ++ci)
              {
                auto it = std::find(dofSets[i].begin(),dofSets[i].end(),ci.index());
                if (it!=dofSets[i].end())
                {
                  int jloc = it-dofSets[i].begin();
                  b[con] -= *ci * dxloc[jloc];
                }
              }
            }
          }
        },nTasks);
      }
    }

  private:
    MatrixA const* A;
    MatrixB const* B;

    std::vector<typename MatrixA::value_type> unconDiagAinv;

    using DofSet = std::vector<size_t>;
    std::vector<DofSet> dofSets;
    std::vector<DofSet> conSets;
    DofSet              uncon;

  protected:
    std::vector<QPSolver> localQPs;


  private:
    ParallelMode bulkMode        = ParallelMode::SEQUENTIAL;
    ParallelMode constraintsMode = ParallelMode::SEQUENTIAL;
    ParallelMode globalMode      = ParallelMode::PARALLEL;    // only relevant if bulk mode is parallel
    Real qpTol = 1e-8;
  };



  //--------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::QPMultiGridBase(MatrixA A_, MatrixB const& B,
                                                                                std::vector<Prolongation>&& prolongations,
                                                                                Real smootherRegularization, 
                                                                                bool blocks, bool directOnCoarse)
  : mgStackA(std::move(prolongations),std::move(A_),false)
  , dummyLogger(new MGSolverStatistics<d,Real>())
  , logger(dummyLogger.get())
  {
    ScopedTimingSection timer("QPMultiGridBase construction");

    mgStackB.push_back(B);
    for (int l=mgStackA.levels()-1; l>0; --l)
      mgStackB.insert(mgStackB.begin(),mgStackB.front()*asMatrix(mgStackA.p(l-1)));

    // Create smoothers on all levels. On the coarse grid, we use a direct solver as smoother,
    // but on the higher levels a GS/Jacobi smoother.
    for (int l=mgStackA.levels()-1; l>0; --l)
      smoothers.insert(smoothers.begin(),Smoother(mgStackA.a(l),mgStackB[l],false,smootherRegularization,blocks));

    smoothers.insert(smoothers.begin(),CoarseSmoother(mgStackA.coarseGridMatrix(),mgStackB[0],directOnCoarse,smootherRegularization,blocks));

    correctionStack.resize(mgStackA.levels()-1);
  }

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>& 
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::setCoarseLevel(int cl)
  {
    coarseLevel = std::clamp(cl,0,mgStackA.levels()-1);
    return *this;
  }

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>& 
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::setSmoothings(int pr, int po)
  {
    assert(0<=pr && 0<=po && pr+po>0);
    pre = pr;
    post = po;
    return *this;
  }

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>& 
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::setCoarseCorrections(int n)
  {
    assert(n>=1);
    mid = n;
    return *this;
  }

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>&
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::setBulkMode(ParallelMode m)
  {
    for (auto& s: smoothers)
    {
      if constexpr(smoothersEqual)        // if all s in vector are of Smoother type 
        s.setBulkMode(m);
      else        
        std::visit([&m](auto&& s_val){ s_val.setBulkMode(m);},s);     // need std::visit here, since s is of type std::variant<Smoother,CoarseSmoother>
    }
    return *this;
  }

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>& 
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::setLogger(MGSolverStatistics<d,Real>* log)
  {
    if (log)
      logger = log;
    else
      logger = dummyLogger.get();
    return *this;
  }

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::~QPMultiGridBase()
  {
    // There is nothing special here. But we cannot use the default destructor, because it is inline,
    // and this collides with QPPJacobiSmoother2 being defined only in the cpp file.
  }

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  typename QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::MatrixB const&
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::B() const
  {
    return mgStackB[mgStackA.levels()-1];
  }

  // ----------------------------------------------------------------------------------------------

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  std::tuple<typename QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::VectorX,int>
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::solve(VectorX x, VectorX const& c, VectorB const& b,
                                                       double tol, int iter) const
  {
    MatrixA const& A = mgStackA.a(mgStackA.levels()-1);
    MatrixB const& B = mgStackB[mgStackA.levels()-1];
    
    coarseResidual.resize(mgStackA.levels());

    // We want to minimize
    // 1/2 x^T Ax + c^T x s.t. Bx <= b

    std::vector<char> activity(B.N(),'i');

    // TODO: use an adaptive truncation with tolerance-based termination
    for (int i=0; i<iter; ++i) 
    {
      // perform one MG cycle
      auto dx = step(x,c,b);
      x += dx;

      // compute constraints and active set
      // and changes in the active set
      auto con = b;
      B.usmv(-1.0,x,con); // b-B*x
      int activeSetChanges = 0;
      for (int j=0; j<con.N(); ++j)
      {
        char active = con[j][0]<=0? 'a': 'i';
        if (activity[j] != active)
          ++activeSetChanges;
        activity[j] = active;
      }
      
      // compute residual (res1 = energy component, res2 = penalty component)
      auto [res1, res2, cons, nActive] = residual(A,B,c,b,x);
      
      // pass correction info to the logger 
      logger->enterVCycleCorrection(i,dx,computeEnergy(A,B,c,b,x),res1,res2,cons,activeSetChanges,correctionStack,coarseResidual);  //activeSetChanges insted of nActive
    }

    return std::make_tuple(x,iter);
  }

  // ----------------------------------------------------------------------------------------------

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  typename QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::VectorX
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::step(VectorX const& x, VectorX c, VectorB b) const
  {
    ScopedTimingSection timer("QPMultiGridBase solve");
    MatrixA const& A = mgStackA.a(mgStackA.levels()-1);
    MatrixB const& B = mgStackB[mgStackA.levels()-1];

    // In terms of corrections dx, we aim at
    // 1/2 (x+dx)^T A (x+dx) + c^T(x+dx) s.t. B (x+dx) <= b or, equivalently, minimize
    // 1/2 dx^T A dx + (Ax+c)^T dx s.t. Bdx <= (b-Bx).
    // First compute the new vector quantities c <- Ax+c, b <- b-Bx
    A.umv(x,c);
    B.usmv(-1.0,x,b);

    // perform V-cycle
    VectorX dx = mgRecursion(mgStackA.levels()-1,pre,post,c,b);
    double t = qpLinesearch(A,B,c,b,dx);
    dx *= t;
    
    return dx;
  }

  // ----------------------------------------------------------------------------------------------

  template <int d, class Prolongation, class Smoother, class CoarseSmoother, class Real>
  typename QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::VectorX
  QPMultiGridBase<d,Prolongation,Smoother,CoarseSmoother,Real>::mgRecursion(int level, int lpre, int lpost,
                                                             VectorX c, VectorB b) const
  {
    MatrixA const& A = mgStackA.a(level);
    MatrixB const& B = mgStackB[level];

    VectorX x(A.N()), dx(A.N());

    // Choose initial solution. Since we compute a correction here, we start with zero
    x = 0;

    std::string timerlevel = "level " + std::to_string(level);

    auto update = [&]()
    {
      double t = qpLinesearch(A,B,c,b,dx);
      x.axpy(t,dx);
      A.usmv(t,dx,c);
      B.usmv(-t,dx,b);
    };
    
    // prolongate a correction up to
    // finest level (for output)
    auto prolongate = [&](VectorX const& dx)
    {
      // if on finest level, return correction
      if(level == mgStackA.levels()-1)
        return dx;
      else
      {
        auto corr = dx;
        int currLevel = level;
        while(currLevel < mgStackA.levels()-1)   
        {
          // prolongate correction from currLevel to currLevel+1
          auto const& P = mgStackA.p(currLevel);
          VectorX Pcorr(P.N());
          P.mv(corr,Pcorr);
          corr = Pcorr;
          ++currLevel;
        }
        return corr;
      }
    };


    // Presmoothing. 
    // On the coarse level, we use a direct solver (see QPMultiGridBase constructor).
    // Make sure it is actually called (at least once)
    if (level==coarseLevel)
      lpre = std::max(lpre,1);
      
    for (int j=0; j<lpre; ++j)
    {
      ScopedTimingSection tim(timerlevel+" smoother");
      // Apply (block) Jacobi smoother and update (dx is initialized in smoother).
      try
      {
        if constexpr(smoothersEqual)                        // if all smoothers are of type Smoother
          smoothers[level].applyCorrection(c,b,dx);
        else
          std::visit([&c,&b,&dx](auto&& s_val){s_val.applyCorrection(c,b,dx);},smoothers[level]);          // need std::visit here, since smoothers[level] is of type std::variant<Smoother,CoarseSmoother>
      }
      catch(IterativeSolverException const& ex)
      {
        std::string message = "Presmoother failed on MG level " + std::to_string(level) + ".\n" + ex.what();
        throw IterativeSolverException(message,__FILE__,__LINE__);
      }
      logger->enterPreSmoothingCorrection(level,dx);
// compute residual on coarselevel
auto [res1, res2, cons, nActive] = residual(A,B,c,b,dx);
auto res = res1; res +=res2;
coarseResidual.at(level) = std::make_tuple(res.two_norm(),res1.two_norm(),res2.two_norm());
      update();
    }

    if (level==coarseLevel)                 // no multigrid recursion below coarse grid
      return x;
    
    // Coarse grid correction. With dx = Pdx we minimize
    // 1/2 dx^T (P^TAP)dx + (P^Tc)^T dx s.t. BP dx <= b
    // P^TAP and BP are available in the multigrid stacks. So here we only need to
    // compute P^T c.
    auto const& P = mgStackA.p(level-1);
    VectorX ptc(P.M());

    for (int i=0; i<mid; ++i)
    {
      P.mtv(c,ptc);
      auto dxCoarse = mgRecursion(level-1,lpre,lpost,ptc,b);
      P.mv(dxCoarse,dx);
      update();
    }
    
    // prolongate coarse grid correction up to fine level for output
    correctionStack.at(level-1) = prolongate(dx);
    //  logger->enterCoarseGridCorrection(level,dx));

    // Postsmoothing
    for (int j=0; j<lpost; ++j)
    {
      ScopedTimingSection tim(timerlevel+" smoother");

      // Apply (block) Jacobi smoother and update (dx is initialized in smoother).
      try
      {
        if constexpr(smoothersEqual)
          smoothers[level].applyCorrection(c,b,dx);                        // if all smoothers are of type Smoother
        else
          std::visit([&c,&b,&dx](auto&& s_val){s_val.applyCorrection(c,b,dx);},smoothers[level]);          // need std::visit here, since s is of type std::variant<smoother,CoarseSmoother>     
      }
      catch(IterativeSolverException const& ex)
      {
        std::string message = "Postsmoother failed on MG level " + std::to_string(level) + ".\n" + ex.what();
        throw IterativeSolverException(message,__FILE__,__LINE__);
      }
      
// compute residual on coarselevel
auto [res1, res2, cons, nActive] = residual(A,B,c,b,dx);
auto res = res1; res +=res2;
coarseResidual.at(level) = std::make_tuple(res.two_norm(),res1.two_norm(),res2.two_norm());
      
      logger->enterPostSmoothingCorrection(level,dx);
      update();
    }

    return x;
  }

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------


  template <int d, class Real, QPStructure sparsity>
  class QPPenaltySmoother: public QPSmootherBase<d,QPPenalizedSolver<d,sparsity,Real>,Real>
  {
    using Base = QPSmootherBase<d,QPPenalizedSolver<d,sparsity,Real>,Real>;
  public:
    using typename Base::MatrixA;
    using typename Base::MatrixB;
    using Self = QPPenaltySmoother<d,Real,sparsity>;

    QPPenaltySmoother(MatrixA const& A_, MatrixB const& B_, bool direct=false, Real regularization=0, bool blocks=true)
    : Base(A_,B_,direct,regularization,blocks)
    {
      // The penalized solver needs to make at least *some* progress.
      // Thus, we force it to take at least one step, even if the tolerance
      // is already satisfied.
      for (auto& lqp: this->localQPs)
        lqp.setMinSteps(1);

      this->setQPtol(0.0); // perform an exact linesearch on penalized QPs.
    }

    Self& setPenalty(Real g)
    {
      assert(g>=0);
      for (auto& lqp: this->localQPs)
        lqp.setPenalty(g);
      return *this;
    }
    
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::QPMultiGrid(MatrixA const& A, MatrixB const& B,
                                                                        std::vector<Prolongation>&& prolongations,
                                                                        Real smootherRegularization, 
                                                                        bool blocks, bool directOnCoarse)
  : Base(A,B,std::move(prolongations),smootherRegularization,blocks,directOnCoarse)
  {}

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::~QPMultiGrid()
  {
    // There is nothing special here. But we cannot use the default destructor, because it is inline,
    // and this collides with the Smoother being defined only in the cpp file.
  }

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>& QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::setPenalty(Real g, bool levelDependent)
  {
    assert(g>=0);
    gamma = g;
    int levels = (this->smoothers).size();
    for (size_t i=0; i<levels; ++i)
    {
      if(levelDependent && i>0)
        g = g/4;
      std::cout << "set penalty = " << g << " on level " << levels-1-i << "\n";
      if constexpr (this->smoothersEqual)                        // if all smoothers are of type Smoother
        (this->smoothers)[levels-1-i].setPenalty(g);
      else
        std::visit([&g](auto&& s_val){s_val.setPenalty(g);},(this->smoothers)[levels-1-i]);          // need std::visit here, since s is of type std::variant<Smoother,CoarseSmoother>
    }
    return *this;
  }
  
  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  double QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::getPenalty()
  {
    return gamma;
  }

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>& QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::setConstraintsMode(ParallelMode m)
  {
    for (auto& s: this->smoothers)
      if constexpr (this->smoothersEqual)                    // if all s are of type Smoother
        s.setConstraintsMode(m);
      else
        std::visit([&m](auto&& s_val){s_val.setConstraintsMode(m);},s);          // need std::visit here, since s is of type std::variant<Smoother,CoarseSmoother>
    return *this;
  }

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>& QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::setGlobalMode(ParallelMode m)
  {
    for (auto& s: this->smoothers)
      if constexpr (this->smoothersEqual)             // if all smoothers are of type Smoother
        s.setGlobalMode(m);
      else
        std::visit([&m](auto&& s_val){s_val.setGlobalMode(m);},s);          // need std::visit here, since s is of type std::variant<Smoother,CoarseSmoother>
    return *this;
  }


  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  std::tuple<typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorB,int>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::solve(VectorX x, VectorX const& c, VectorB const& b, VectorB const& lambda,
                                           double tol, int iter) const
  {
    // We want to minimize
    // 1/2 x^T Ax + c^T x + gamma/2 |Bx-(b-lambda/gamma)|_+^2
    // and set bl <- b-lambda/gamma
    auto bl = b;
    bl.axpy(-1.0/gamma,lambda);

    // solve and multiplier update
    std::tie(x,iter) = Base::solve(x,c,bl,tol,iter);
    auto dLam = firstOrderUpdate(x,b,lambda);

    return std::tuple(x,dLam,iter);
  }
  
  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorB 
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::firstOrderUpdate(VectorX const& x, VectorB const& b, VectorB const& lambda) const
  {
    MatrixB const& B = this->B();

    VectorB dLambda(b.N());
    B.mv(x,dLambda);                                   // B*x
    dLambda.axpy(-1.0,b);                              // B*x-b
    for (int i=0; i<b.N(); ++i)                        
      dLambda[i] = std::max(-1.0*lambda[i][0],gamma*dLambda[i][0]);             // dLambda = max(-lambda, gamma*(B*x-b))

    return dLambda;
  }

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  double QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::computeEnergy(MatrixA const& A, MatrixB const& B,
                                                         VectorX const& c, VectorB const& b, VectorX const& x) const
  {
    VectorX Ax(x);
    A.mv(x,Ax);
    
    auto [grad,con,nActive] = gradient(A,B,c,b,x);

    Real energy1 = 0.5*(Ax*x) + c*x;
    Real energy2 = gamma*(con*con)/2;

    std::cout.precision(16);
    std::cout << "energy = " << energy1+energy2 << ", |grad| = " << grad.two_norm() << ", |x| = " << x.two_norm() << ", active=" << nActive << ", penalty = " << energy2  <<"\n";

    return energy1+energy2;
  }
  
  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  std::tuple<typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorB,int>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::gradient(MatrixA const& A, MatrixB const& B,
                                             VectorX const& c, VectorB const& b, VectorX const& x) const
  {
    // the gradient of objective function
    // 1/2 x^T Ax + c^T x + gamma/2 |Bx-b|_+^2
    // is given by Ax + c + gamma* B^T(Bx-b)_+
    VectorX Ax(x);
    A.mv(x,Ax);
    VectorB Bx(b);
    B.mv(x,Bx);
    
    // con = (Bx-b)_+
    VectorB con = Bx;
    con -= b;
    int nActive = 0;
    for (int i=0; i<con.N(); ++i)
    {
      if (con[i][0]>0)
        ++nActive;
      con[i] = std::max(0.0,con[i][0]);
    }
    
    VectorX grad = Ax;
    grad += c;
    B.usmtv(gamma,con,grad);
    
    return std::make_tuple(grad,con,nActive);
  }
  
  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  std::tuple<typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,
             typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,
             typename QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorB,int>
  QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::residual(MatrixA const& A, MatrixB const& B,
                                             VectorX const& c, VectorB const& b, VectorX const& x) const
  {
    // the gradient of objective function
    // 1/2 x^T Ax + c^T x + gamma/2 |Bx-b|_+^2
    // is given by Ax + c + gamma* B^T(Bx-b)_+
    VectorX Ax(x);
    A.mv(x,Ax);
    VectorB Bx(b);
    B.mv(x,Bx);
    
    // con = (Bx-b)_+
    VectorB con = Bx;
    con -= b;
    int nActive = 0;
    for (int i=0; i<con.N(); ++i)
    {
      if (con[i][0]>0)
        ++nActive;
      //  else
        //  con[i][0] = 0;
      con[i] = std::max(0.0,con[i][0]);
    }
    
    VectorX grad1 = Ax;
    grad1 += c;
    VectorX grad2(x);
    B.smtv(gamma,con,grad2);
    
    return std::make_tuple(grad1,grad2,con,nActive);
  }


  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  double QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::qpLinesearch(MatrixA const& A, MatrixB const& B, VectorX const& c,
                                                         VectorB const& b, VectorX const& dx) const
  {
    return Kaskade::qpLinesearch(A,B,gamma,c,b,dx);
  }


//   template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
//   void QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::leastSquaresMultiplierUpdate(VectorX const& c, VectorB const& b, VectorX const& x, VectorB& lambda) const
//   {
//     MatrixA const& A = mgStackA.a(mgStackA.levels()-1);
//     MatrixB const& B = mgStackB[mgStackA.levels()-1];
//
//     // The KKT adjoint equation for the augmented Lagrangian problem
//     // min 1/2 x^T Ax + c^T x + gamma/2 (Bx-b-s)^T (Bx-b-s) - lambda^T (Bx-b-s)
//     // is Ax+c + gamma*B^T(Bx-b-s) - B^T*lambda = 0.
//     // Let r = Ax+c + gamma*B^T(Bx-b-s). Choosing lambda such that
//     // this is "minimal" (in which norm?) amounts to solving the least squares
//     // system min |B^T lambda - r|, e.g., by solving the normal equations
//     // B B^T lambda = B r. This is what we do here.
//
//     // First compute r = Ax+c
//     VectorX r = c;
//     A.umv(x,r);
//
//     // s is not provided, but we can choose the optimal s = min(0,Bx-b) on the fly:
//     // Then q = Bx-b-s = max(0,Bx-b).
//     VectorB q(b.N());
//     B.mv(x,q);
//     q -= b;
//     for (size_t i=0; i<b.N(); ++i)
//       q[i] = std::max(Real(0),q[i][0]);
//
//     B.usmtv(gamma,q,r);         // now r = Ax+c + gamma*B^T(Bx-b-s)
//
//     auto tmp = r;
//     B.usmtv(-1.0,lambda,tmp); // residual
//     std::cout << "least squares multiplier update: residual is " << tmp.two_norm() << "\n";
//
//     // Now we solve B B^T lambda = B r with conjugate gradient method.
//     B.mv(r,q);      // q = Br
//
//     DefaultDualPairing<VectorB,VectorB> dp;
//     NumaBCRSMatrix<Dune::FieldMatrix<Real,d,1>> Bt(B,false,true); // transpose constructor
//     using BBtMatrix = NumaBCRSMatrix<Dune::FieldMatrix<Real,1,1>>;
//     BBtMatrix BBtm = B*Bt;
//     SymmetricMatrixOperator<VectorB,BBtMatrix> BBt(BBtm);
//     IdentityPreconditioner<VectorB> preco;
//     PCGCountingTerminationCriterion<Real> terminate(10);
//     Pcg<VectorB,VectorB> pcg(BBt,preco,terminate);
//
//     Dune::InverseOperatorResult res;
//     pcg.apply(lambda,q,res);
//
//     tmp = r;
//     B.usmtv(-1.0,lambda,tmp); // residual
//     std::cout << "least squares multiplier update: residual is now " << tmp.two_norm() << "\n";
//   }
//
//   template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
//   void QPMultiGrid<d,Prolongation,Real,Smoother,CoarseSmoother>::firstOrderMultiplierUpdate(VectorX const& c, VectorB const& b,
//                                                                     VectorX const& x, VectorB& lambda) const
//   {
//     MatrixB const& B = mgStackB[mgStackA.levels()-1];
//     VectorB tmp(b.N());
//     B.mv(x,tmp); tmp -= b;
//     VectorB s(b.N());
//     for (size_t i=0; i<s.N(); ++i)
//       s[i] = std::min(Real(0),tmp[i][0]-lambda[i][0]/gamma);
//     std::cout << "first order update: |s| = " << s.two_norm() << "\n";
//     tmp -= s;
//     std::cout << "first order update: constraint violation = " << tmp.two_norm() << "\n";
//     tmp *= gamma;
//     std::cout << "lambda update: " << tmp.two_norm() << "\n";
//     lambda -= tmp;
//   }

  template <int d, QPStructure sparsity>
  using QPSv = QPPenaltySmoother<d,double,sparsity>;


  template class QPMultiGridBase<1,ProlongationMatrix,QPSv<1,QPStructure::DENSE>,QPSv<1,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<1,ProlongationMatrix,QPSv<1,QPStructure::DENSE>,QPSv<1,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<2,ProlongationMatrix,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<2,ProlongationMatrix,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<3,ProlongationMatrix,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<3,ProlongationMatrix,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<2,MGProlongation,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<2,MGProlongation,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<3,MGProlongation,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<3,MGProlongation,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::SPARSE>,double>;


  template class QPMultiGrid<1,ProlongationMatrix,double,QPSv<1,QPStructure::DENSE>,QPSv<1,QPStructure::DENSE>>;
  template class QPMultiGrid<1,ProlongationMatrix,double,QPSv<1,QPStructure::DENSE>,QPSv<1,QPStructure::SPARSE>>;
  template class QPMultiGrid<2,ProlongationMatrix,double,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::DENSE>>;
  template class QPMultiGrid<2,ProlongationMatrix,double,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::SPARSE>>;
  template class QPMultiGrid<3,ProlongationMatrix,double,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::DENSE>>;
  template class QPMultiGrid<3,ProlongationMatrix,double,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::SPARSE>>;
  template class QPMultiGrid<2,MGProlongation,double,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::DENSE>>;
  template class QPMultiGrid<2,MGProlongation,double,QPSv<2,QPStructure::DENSE>,QPSv<2,QPStructure::SPARSE>>;
  template class QPMultiGrid<3,MGProlongation,double,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::DENSE>>;
  template class QPMultiGrid<3,MGProlongation,double,QPSv<3,QPStructure::DENSE>,QPSv<3,QPStructure::SPARSE>>;

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <int d, class Real, QPStructure sparsity>
  class QPSmoother: public QPSmootherBase<d,QPALSolver<d,sparsity,Real>,Real>
  {
    using Base = QPSmootherBase<d,QPALSolver<d,sparsity,Real>,Real>;
  public:
    using typename Base::MatrixA;
    using typename Base::MatrixB;
    using Self = QPSmoother<d,Real,sparsity>;

    QPSmoother(MatrixA const& A_, MatrixB const& B_, bool direct=false, Real regularization=0, bool blocks=true)
    : Base(A_,B_,direct,regularization,blocks)
    {
      // make sure that all dofs affecting constraints are processed sequentially
      this->setConstraintsMode(ParallelMode::SEQUENTIAL);
      this->setGlobalMode(ParallelMode::SEQUENTIAL); // maybe this is too restrictive? Can we relax this?

      for (auto& lqp: this->localQPs)
        lqp.setMaxSteps(200);

      this->setQPtol(1e-8);
    }
    
    Self& setPenalty(Real g)
    {
      assert(g>=0);
      for (auto& lqp: this->localQPs)
        lqp.setPenalty(g);
      return *this;
    }
    
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::QPMultiGridStrict(MatrixA const& A, MatrixB const& B,
                                                            std::vector<Prolongation>&& prolongations,
                                                            Real smootherRegularization, 
                                                            bool blocks, bool directOnCoarse)
  : Base(A,B,std::move(prolongations),smootherRegularization,blocks,directOnCoarse)
  {}

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::~QPMultiGridStrict()
  {
    // There is nothing special here. But we cannot use the default destructor, because it is inline,
    // and this collides with QPPJacobiSmoother2 being defined only in the cpp file.
  }
  
  
  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  std::tuple<typename QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,int>
  QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::solve(VectorX x, VectorX const& c, VectorB const& b, double tol, int iter) const
  {
    return Base::solve(x,c,b,tol,iter);
  }
  

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  double QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::computeEnergy(MatrixA const& A, MatrixB const& B,
                                                                                       VectorX const& c, VectorB const& b, VectorX const& x) const
  {
    VectorX Ax(x);
    A.mv(x,Ax);

    auto [grad,con,nActive] = gradient(A,B,c,b,x);

    Real energy1 = 0.5*(Ax*x) + c*x;
    
    double violation = con[0][0];
    for (int i=1; i<con.N(); ++i)
    {
      if (con[i][0]>violation)
        violation = con[i][0];
    }

    std::cout.precision(16);
    std::cout << "energy = " << energy1 << ", max Bx-b = " << violation
              << ", |x| = " << x.two_norm() << ", active=" << nActive << "\n";

    return energy1;
  }
  
  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  std::tuple<typename QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,typename QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorB,int>
  QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::gradient(MatrixA const& A, MatrixB const& B,
                                                   VectorX const& c, VectorB const& b, VectorX const& x) const
  {
    // the gradient of objective function
    // 1/2 x^T Ax + c^T x
    // is given by Ax + c
    VectorX Ax(x);
    A.mv(x,Ax);
    
    VectorB Bx(b);
    B.mv(x,Bx);
    
    // con = (Bx-b)_+
    VectorB con = Bx;
    con -= b;
    int nActive = 0;
    for (int i=0; i<con.N(); ++i)
    {
      if (con[i][0]>0)
        ++nActive;
      con[i] = std::max(0.0,con[i][0]);
    }
    
    VectorX grad = Ax;
    
    return std::make_tuple(grad,con,nActive);
  }
  
  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  std::tuple<typename QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,
             typename QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorX,
             typename QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::VectorB,int>
  QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::residual(MatrixA const& A, MatrixB const& B,
                                                   VectorX const& c, VectorB const& b, VectorX const& x) const
  {
    return std::make_tuple(x,x,b,0);
  }

  template <int d, class Prolongation, class Real, class Smoother, class CoarseSmoother>
  double QPMultiGridStrict<d,Prolongation,Real,Smoother,CoarseSmoother>::qpLinesearch(MatrixA const& A, MatrixB const& B, VectorX const& c,
                                                              VectorB const& b, VectorX const& dx) const
  {
    // We do a Gauss-Seidel scheme here, with (almost) exact preservation of feasibility
    // by the smoother steps. Thus, we can always accept the smoother result.
    return 1;
  }

  template <int d,QPStructure sparsity>
  using QPSm = QPSmoother<d,double,sparsity>;

  template class QPMultiGridBase<1,ProlongationMatrix,QPSm<1,QPStructure::DENSE>,QPSm<1,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<1,ProlongationMatrix,QPSm<1,QPStructure::DENSE>,QPSm<1,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<2,ProlongationMatrix,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<2,ProlongationMatrix,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<3,ProlongationMatrix,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<3,ProlongationMatrix,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<2,MGProlongation,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<2,MGProlongation,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::SPARSE>,double>;
  template class QPMultiGridBase<3,MGProlongation,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::DENSE>,double>;
  template class QPMultiGridBase<3,MGProlongation,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::SPARSE>,double>;


  template class QPMultiGridStrict<1,ProlongationMatrix,double,QPSm<1,QPStructure::DENSE>,QPSm<1,QPStructure::DENSE>>;
  template class QPMultiGridStrict<1,ProlongationMatrix,double,QPSm<1,QPStructure::DENSE>,QPSm<1,QPStructure::SPARSE>>;
  template class QPMultiGridStrict<2,ProlongationMatrix,double,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::DENSE>>;
  template class QPMultiGridStrict<2,ProlongationMatrix,double,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::SPARSE>>;
  template class QPMultiGridStrict<3,ProlongationMatrix,double,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::DENSE>>;
  template class QPMultiGridStrict<3,ProlongationMatrix,double,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::SPARSE>>;
  template class QPMultiGridStrict<2,MGProlongation,double,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::DENSE>>;
  template class QPMultiGridStrict<2,MGProlongation,double,QPSm<2,QPStructure::DENSE>,QPSm<2,QPStructure::SPARSE>>;
  template class QPMultiGridStrict<3,MGProlongation,double,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::DENSE>>;
  template class QPMultiGridStrict<3,MGProlongation,double,QPSm<3,QPStructure::DENSE>,QPSm<3,QPStructure::SPARSE>>;

}
