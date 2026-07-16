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

#ifndef LINALG_PATCHDOMAINDECOMPOSITIONPPRECONDITIONER_HH
#define LINALG_PATCHDOMAINDECOMPOSITIONPPRECONDITIONER_HH

#include <vector>

#include <boost/timer/timer.hpp>

#include <dune/common/dynvector.hh>
#include "dune/istl/preconditioner.hh"

#include "linalg/dynamicMatrix.hh"
#include "linalg/localMatrices.hh"
#include "linalg/symmetricOperators.hh"
#include "linalg/threadedMatrix.hh"
#include "fem/supports.hh"
#include "utilities/scalar.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"


// #define INEXACTSTORAGE
// #include <random> // for mixed precision *testing*
// extern double epsilon;



namespace Kaskade
{
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup preconditioners
   * \brief A base class providing common functionality for storage of inverse local matrices.
   */
  class InverseLocalMatrixStorageBase
  {
  public:
    InverseLocalMatrixStorageBase(std::vector<size_t> const& ap)
    : dofs(ap)
    {
    }

    /**
     * \brief Returns the global degrees of freedom defining the submatrix treated as a block here.
     */
    std::vector<size_t> const& indices() const
    {
      return dofs;
    }

    /**
     * \brief Reports the storage size (in bytes) of the inverse block matrix.
     *
     * Minor metadata that does not scale with the matrix size is neglected.
     */
    size_t storage() const
    {
      return storageSize;
    }

    /**
     * \brief Reports the number of floating point operations per preconditioner application.
     *
     * This is a (rather accurate) estimate. Depending on BLAS implementation, lower order contributions
     * may be inaccurate.
     */
    size_t flops() const
    {
      return computeFlops;
    }

  protected:
    size_t storageSize;
    size_t computeFlops;

    template <int m, class LocalMatrix, class Real>
    void regularize(LocalMatrix& aloc, Real r) const
    {
      if (r > 0)
      {
        Real Amax = 0;
        for (int i=0; i<aloc.N(); ++i)
          Amax = std::max(static_cast<Real>(std::sqrt(frobenius_norm2(aloc[i][i]))),Amax);
        for (int i=0; i<aloc.N(); ++i)
          aloc[i][i] += Amax/r * unitMatrix<Real,m>();
      }
    }

  private:
    std::vector<size_t> dofs;
  };

  // ----------------------------------------------------------------------------------------------

  /**
   * \brief A parametrized tag type for selecting dense storage of inverses with floating point representation of type Scalar.
   * \relates DenseInverseLocalMatrixStorage
   * \relates PatchDomainDecompositionPreconditioner
   */
  template <class Scalar>
  class DenseInverseStorageTag {};

  /**
   * \ingroup preconditioners
   * \brief An inverse local matrix representation that simply stores the dense inverse.
   * \relates PatchDomainDecompositionPreconditioner
   */
  template <int m, class Scalar=float>
  class DenseInverseLocalMatrixStorage: public InverseLocalMatrixStorageBase
  {
    using StorageEntry = Dune::FieldMatrix<Scalar,m,m>;

  public:

    using field_type = Scalar;

    template <class T, class Index>
    DenseInverseLocalMatrixStorage(std::vector<size_t> const& ap,
                                   NumaBCRSMatrix<Dune::FieldMatrix<T,m,m>,Index> const& A,
                                   field_type regularizeForCondition, DenseInverseStorageTag<Scalar>)
    : InverseLocalMatrixStorageBase(ap)
    {
      alocinv = full(A,indices(),indices());

#ifdef INEXACTSTORAGE
      using Real = typename ScalarTraits<Scalar>::Real;
      Real anorm = two_norm(alocinv);
#endif

      regularize<m>(alocinv,regularizeForCondition);

      invertSpd(alocinv); // invert the local block

#ifdef INEXACTSTORAGE
      DynamicMatrix<Dune::FieldMatrix<Scalar,m,m>> delta(alocinv.N(),alocinv.M());

      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dis(-1.0,1.0);
      for (int i=0; i<alocinv.N(); ++i)
        for (int j=0; j<alocinv.M(); ++j)
          for (int k=0; k<m; ++k)
            for (int l=0; l<m; ++l)
              delta[i][j][k][l] = dis(gen);
      Real dnorm = two_norm(delta);

      Real epsilon = ::epsilon;
      delta *= epsilon / anorm / dnorm;
      alocinv += delta;
#endif

      storageSize = sizeof(StorageEntry)*alocinv.N()*alocinv.M() + sizeof(size_t)*indices().size();
      auto n = alocinv.N() * m;
      computeFlops = n*(2*n-1); // in every row n multiplications and n-1 additions
    }


    template <class Vector>
    void mv(Vector const& x, Vector& y) const
    {
      alocinv.mv(x,y);
    }


  private:
    DynamicMatrix<StorageEntry> alocinv;
  };

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup solvers
   * \brief Computes a nested dissection of the given graph.
   * \param edges a sequence \f$ E \f$ of edges \f$ (v,w) \f$, given as pairs of vertices, sorted lexicographically.
   *              The edges are assumed to be symmetric, i.e. \f$ (v,w)\in E \Rightarrow (w,v)\in E \f$, and shall
   *              not contain the diagonal, i.e. \f$ (v,v) \not\in E \f$.
   * \param edgeStart an index sequence such that all edges \f$ (v,w) \f$ for some particular \f$ v \f$
   *                  are in the half-open range [edgeStart(v),edgeStart(v+1)[ of edges.
   *
   * The vertices are assumed to be a contiguous range \f$ 0, \dots, n-1 \f$ of integers.
   *
   * \return a tuple \f$ (v_1, v_2, s) \f$ of vertices, forming the two sets \f$ v_1, v_2 \f$ separated by \f$ s \f$.
   *
   * If no proper partitioning can be found, both \f$ v_1, v_2 \f$ are empty and \f$ s \f$ contains all vertices.
   * If the graph is disconnected, \f$ s \f$ is empty.
   */
  std::array<std::vector<size_t>,3> nestedDissection(std::vector<std::pair<size_t,size_t>> const& edges,
                                                     std::vector<size_t> const& edgeStart);

  /**
   * \ingroup solver
   * \brief a nested dissection factorization for symmetric positive definite sparse matrices
   * \tparam StorageScalar the floating point type used for storage of Cholesky factors
   */
  template <class StorageScalar>
  class NestedDissection
  {
  public:
    NestedDissection() = default;

    /**
     * \brief Constructor.
     */
    template <class SparseMatrix>
    NestedDissection(SparseMatrix const& A, size_t skipEntries, size_t skipDimension);

    /**
     * \brief move assignment from different storage scalar type
     *
     * Use this for compression after construction. E.g., the construction (including factorization)
     * can be done with double precision, and the Cholesky factor be converted to single precision
     * afterwards using this move assignment.
     */
    template <class FactorizationScalar>
    NestedDissection<StorageScalar>& operator=(NestedDissection<FactorizationScalar>&& nd);

    template <class T>
    void solve(T& y) const;

    size_t storageSize() const;
    size_t flops() const;

    /**
     * \brief reports the total number of separators in this nested dissection
     */
    size_t size()
    {
      return sz;
    }

  private:
    template <class FactorizationScalar>
    friend class NestedDissection;

    std::vector<size_t> v1, v2, s;
    std::unique_ptr<NestedDissection<StorageScalar>> l1FactorND, l2FactorND;
    DynamicMatrix<StorageScalar> r1, r2;
    std::vector<StorageScalar> sFactor;
    size_t sz;

    /**
     * \brief in-place solve of \f$ y \gets L^{-1} y \f$
     * \param trans if 'n', \f$ L \f$ is used, 'c' uses \f$ L^H \f$, and
     *              'b' performs the complete Cholesky solve \f$ y \gets (LL^H)^{-1} y \f$, which
     *              is equivalent to but (slightly) more efficient than calling
     *              trisolve(y,'n'); trisolve(y,'c').
     */
    template <int m>
    void trisolve(DynamicMatrix<StorageScalar>& y, char trans) const;
  };

  /**
   * \ingroup preconditioners
   * \brief A parametrized tag type for selecting storage of inverses by nested dissection
   *        with floating point representation of type Scalar.
   * \relates NestedDissectionCholeskyLocalMatrixStorage
   * \relates PatchDomainDecompositionPreconditioner
   */
  template <class StorageScalar, class FactorizationScalar>
  class NestedDissectionStorageTag
  {
  public:

    NestedDissectionStorageTag(size_t skipEntries_=100, size_t skipDimension_=32)
    : skipEntries(skipEntries_), skipDimension(skipDimension_)
    {}

    /**
     * \brief omit further dissection if the size of the zero-block (in terms of entries) is below this number
     */
    size_t skipEntries;

    /**
     * \brief omit further dissection if the size of the system is below this number
     */
    size_t skipDimension;
  };

  /**
   * \ingroup preconditioners
   * \brief An inverse local matrix representation that stores a sparse Cholesky factor using nested dissection.
   * \relates PatchDomainDecompositionPreconditioner
   */
  template <int m, class StorageScalar, class FactorizationScalar>
  class NestedDissectionCholeskyLocalMatrixStorage: public InverseLocalMatrixStorageBase
  {
  public:
    using field_type = StorageScalar;

    template <class T, class Index>
    NestedDissectionCholeskyLocalMatrixStorage(std::vector<size_t> const& ap,
                                               NumaBCRSMatrix<Dune::FieldMatrix<T,m,m>,Index> const& A,
                                               FactorizationScalar regularizeForCondition,
                                               NestedDissectionStorageTag<StorageScalar,FactorizationScalar> storageTag)
    : InverseLocalMatrixStorageBase(ap)
    {
      using LocalMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<T,m,m>>;

      LocalMatrix aloc = submatrix<LocalMatrix>(A,indices(),indices());

      regularize<m>(aloc,regularizeForCondition);

      choleskyFactor = NestedDissection<FactorizationScalar>(aloc,storageTag.skipEntries,
                                                             storageTag.skipDimension);  // factorize the local block and compress

      storageSize = choleskyFactor.storageSize();
      computeFlops = choleskyFactor.flops();
    }


    /**
     * \tparam Vector a vector class with entry type Dune::FieldVector<StorageScalar,m> (or compatible)
     */
    template <class Vector>
    void mv(Vector const& x, Vector& y) const
    {
      y = x;
      choleskyFactor.solve(y);
    }

  private:
    NestedDissection<StorageScalar> choleskyFactor;
  };


  // ----------------------------------------------------------------------------------------------

  namespace DomainDecompositionPreconditionerDetail
  {
    template <class Scalar>
    std::vector<Scalar> choleskyFactor(DynamicMatrix<Scalar> const& A);

    template <class Scalar>
    void choleskySolve(int n, int nrhs, std::vector<Scalar> const& L, Scalar* b);
  }

  /**
   * \brief A parametrized tag type for selecting storage of inverses by dense Cholesky factors
   *        with floating point representation of type Scalar.
   * \relates DenseCholeskyLocalMatrixStorage
   * \relates PatchDomainDecompositionPreconditioner
   */
  template <class Scalar>
  class DenseCholeskyStorageTag {};

  /**
   * \ingroup preconditioners
   * \brief An inverse local matrix representation that simply stores a dense Cholesky factor.
   * \relates PatchDomainDecompositionPreconditioner
   */
  template <int m, class Scalar=float>
  class DenseCholeskyLocalMatrixStorage: public InverseLocalMatrixStorageBase
  {
    using StorageEntry = Dune::FieldMatrix<Scalar,m,m>;

  public:

    using field_type = Scalar;

    template <class T, class Index>
    DenseCholeskyLocalMatrixStorage(std::vector<size_t> const& ap,
                                    NumaBCRSMatrix<Dune::FieldMatrix<T,m,m>,Index> const& A,
                                    field_type regularizeForCondition, DenseCholeskyStorageTag<Scalar>)
    : InverseLocalMatrixStorageBase(ap)
    {
      DynamicMatrix<Scalar> aloc = DynamicMatrixDetail::flatMatrix(full(A,indices(),indices()));
      n = aloc.N();

      regularize<1>(aloc,regularizeForCondition);

      choleskyFactor = DomainDecompositionPreconditionerDetail::choleskyFactor(aloc); // factorize the local block
      storageSize = sizeof(Scalar)*choleskyFactor.size() + sizeof(size_t)*indices().size();
      computeFlops =  2*n*n;

#ifdef INEXACTSTORAGE
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dis(-1.0,0.0);

      DynamicMatrix<Scalar> L(n,n), delta(n,n);
      L = 0;
      delta = 0;
      auto lpos = begin(choleskyFactor);
      for (int j=0; j<n; ++j)
        for (int i=j; i<n; ++i)
        {
          L[i][j] = *lpos;
          ++lpos;
          delta[i][j] = dis(gen) / (i-j+1);
        }

      invert(L);

      Scalar epsilon = ::epsilon;
      delta *= epsilon/two_norm(delta)/two_norm(L);

      lpos = begin(choleskyFactor);
      for (int j=0; j<n; ++j)
        for (int i=j; i<n; ++i)
        {
          *lpos += delta[i][j];
          ++lpos;
        }

#endif
    }


    template <class Vector>
    void mv(Vector const& x, Vector& y) const
    {
      y = x;
      DomainDecompositionPreconditionerDetail::choleskySolve(n,1,choleskyFactor,&y[0][0]);
    }


  private:
    std::vector<Scalar> choleskyFactor; // L in LAPACK packed storage
    int n;                              // size of L (with scalar entries, not blocks of size m)
  };



  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup preconditioners
   * \brief Provides the actual type of local matrix storage based on a tag and the number of components.
   *
   * In order to make a new local matrix storage method available, create a unique tag type and a partial
   * specialization of this class.
   *
   * \relates PatchDomainDecompositionPreconditioner
   */
  template <class StorageTag, int components>
  struct PatchDomainDecompositionPreconditionerTraits
  {
    /**
     * \brief The actual type used for storing local matrices.
     *
     * This has to be redefined in partial specializations.
     */
    using type = void;
  };

  template <class Scalar, int components>
  struct PatchDomainDecompositionPreconditionerTraits<DenseInverseStorageTag<Scalar>,components>
  {
    using type = DenseInverseLocalMatrixStorage<components,Scalar>;
  };

  template <class Scalar, int components>
  struct PatchDomainDecompositionPreconditionerTraits<DenseCholeskyStorageTag<Scalar>,components>
  {
    using type = DenseCholeskyLocalMatrixStorage<components,Scalar>;
  };

  template <class StorageScalar, class FactorizationScalar, int components>
  struct PatchDomainDecompositionPreconditionerTraits<NestedDissectionStorageTag<StorageScalar,FactorizationScalar>,components>
  {
    using type = NestedDissectionCholeskyLocalMatrixStorage<components,StorageScalar,FactorizationScalar>;
  };

  /**
   * \ingroup preconditioners
   * \brief An additive overlapping domain decomposition type preconditioner for higher order
   *        finite elements applied to elliptic equations.
   *
   * The subdomains are the stars around grid vertices. The ansatz functions with a support contained in a subdomain
   * are treated as a block, and the corresponding submatrix is inverted using dense linear algebra. The subdomain
   * corrections are added up.
   *
   * \tparam Space a finite element space (FEFunctionSpace)
   * \tparam m the number of components in the equation
   * \tparam StorageTag a tag class for selecting the method of computing and storing approximate inverse local matrices.
   *                    Common tags are, e.g., DenseInverseStorageTag<double>, DenseCholeskyStorageTag<float>,
   *                    NestedDissectionStorageTag<float,double>.
   */
  template <class Space, int m, class StorageTag,
            class SparseMatrixIndex = std::size_t>
  class PatchDomainDecompositionPreconditioner: public SymmetricPreconditioner<typename Space::template Element<m>::type::StorageType,
                                                                               typename Space::template Element<m>::type::StorageType>
  {
    using InverseLocalMatrix = typename PatchDomainDecompositionPreconditionerTraits<StorageTag,m>::type;

  public:
    using domain_type = typename Space::template Element<m>::type::StorageType;
    using range_type = domain_type;

    /**
     * \brief Constructor 
     * \tparam Matrix the sparse matrix type of the global stiffness matrix
     * 
     * \param A matrix to be preconditioned (symmetric and positive definite)
     * \param regularizeForCondition aim at the given condition number. No regularization if negative values are provided.
     * \param nTasks the number of groups in which to split the patches (each group is processed in one task,
     *               tasks can be executed in parallel). Nonpositive values mean to choos an appropriate default value.
     * 
     * For the patch matrices, the preconditioner can perform regularization that can help if rounding errors
     * interfere with badly conditioned patch matrices, which may cause failure of the Cholesky inversion.
     * If regularizeForCondition is positive, we regularize each patch matrix \f$ A \f$ as follows:
     * \f[ A \leftarrow A + \alpha I, \quad \alpha = \max_{ij} A_{ij} / \mathrm{regularizeForCondition} \f]
     * Reasonable values for regularizeForCondition are above 1000.
     */
    template <class Matrix>
    PatchDomainDecompositionPreconditioner(Space const& space, Matrix const& A, 
                                           typename Matrix::field_type regularizeForCondition = -1.0,
                                           StorageTag storageTag = StorageTag(),
                                           int nTasks = 0)
    {
      assert(A.N()>0);
      assert(A.N()==A.M());
      Timings& timer = Timings::instance();
      timer.start("make Schwarz preconditioner");

      // for each patch around a vertex, store the global dof indices
      std::vector<std::vector<size_t>> patchDofs(space.gridView().size(space.gridView().dimension));
      
      auto patches = computePatches(space.gridView());
      assert(patches.size()==patchDofs.size());
      
      int nThreads = nTasks<=0? 4*NumaThreadPool::instance().cpus()
                              : std::min(nTasks,(int)patches.size());
      inverseLocalMatrices.resize(nThreads);

      TaskTiming taskTiming(nThreads);
      
      // Compute the inverse of each block's local stiffness matrix in parallel.
      parallelFor([&patches,&A,&space,regularizeForCondition,storageTag,&taskTiming,this](size_t j, size_t s)
      {
        taskTiming.start(j);
        // For all patches (in our parallel loop range)
        for (size_t i=uniformWeightRangeStart(j,s,patches.size()); i<uniformWeightRangeStart(j+1,s,patches.size()); ++i)
        {
          auto const& ap = algebraicPatch(space,patches[i]);
          // todo: (catch possible exception, store by std::exception_ptr (which needs to be passed from main thread and locked), then rethrow in main thread
          try 
          {
            InverseLocalMatrix aloc = InverseLocalMatrix(ap,A,regularizeForCondition,storageTag);
            this->inverseLocalMatrices[j].push_back(std::move(aloc));
          } 
          catch (DetailedException const& ex)
          {
            std::cout << "exception in inverting block " << i << " in thread " << j << " of " << s << ":\n";
            std::cout << ex.what() << "\n";
          }
        }
        taskTiming.stop(j);
      },nThreads);

      timer.stop("make Schwarz preconditioner");

std::cerr << "Schwarz preconditioner storage size: " << storage()/(1024*1024) << " MB\n"
          << "                       flop count:   " << flops()/(1024*1024) << " M\n";

      std::ofstream f("out.gnu");
      f << taskTiming;
    }

    virtual typename Space::Scalar applyDp(domain_type& x, range_type const& y)
    {
      apply(x,y);
      return x*y;
    }

    /**
     * \brief Applys the preconditioner as \f$ x \leftarrow P^{-1} y \f$.
     * \param[in] y the right hand side
     * \param[out] x the (approximate) solution
     */
    virtual void apply(domain_type& x, range_type const& y)
    {      
      using StorageScalar = typename InverseLocalMatrix::field_type;
      using LocalVector = Dune::BlockVector<Dune::FieldVector<StorageScalar,m>>;

      // TODO: consider scattering into local copies of global vector? Maybe dependent on the actual
      // local matrix size? Or range locking?
      // But for 3D order 3&4 locking appears not to be a performance bottleneck on
      // a quadcore (as of 2018-11).

      ScopedTimingSection timer("Schwarz preconditioner apply");
      
      std::mutex rhsAccess;
      parallelFor([&rhsAccess,&x,&y,this](int j, int s)
      {
        for (auto const& alocinv: this->inverseLocalMatrices[j])
        {
          auto const& indices = alocinv.indices();
          int const n = indices.size();

          LocalVector xloc(n);
          LocalVector yloc(n);
          
          for (int i=0; i<n; ++i)               // gather relevant vector entries from x
            yloc[i] = y[indices[i]];
          
          alocinv.mv(yloc,xloc);                // perform local multiplication
          
          std::lock_guard<std::mutex> lock(rhsAccess);
          for (int i=0; i<n; ++i)               // scatter local result to global vector y
            x[indices[i]] += xloc[i];
        }
      },inverseLocalMatrices.size());
    }

    virtual bool requiresInitializedInput() const
    {
      return true;
    }

    /**
     * \brief Reports the storage size (in bytes) of the preconditioner.
     *
     * Minor metadata that does not scale with the matrix size is neglected.
     */
    size_t storage() const
    {
      size_t s = 0;
      for (auto const& as: inverseLocalMatrices)
        for (auto const& a: as)
          s += a.storage();
      return s;
    }

    /**
     * \brief reports the number of floating point operations when applying the preconditioner
     */
    size_t flops() const
    {
      size_t f = 0;
      for (auto const& as: inverseLocalMatrices)
      {
        for (auto const& a: as)
          f += a.flops();
      }
      return f;
    }

  private:
    std::vector<std::vector<InverseLocalMatrix>> inverseLocalMatrices;
  };
}

#endif
