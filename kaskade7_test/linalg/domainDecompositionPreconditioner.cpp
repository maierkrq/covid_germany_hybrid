/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2018-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cassert>
#include <limits>
#include <queue>
#include <type_traits>

#include "linalg/domainDecompositionPreconditioner.hh"
#include "linalg/dynamicMatrixOps.hh"
#include "linalg/lapack/lapacke.h"
#include "utilities/detailed_exception.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"


namespace Kaskade
{
  namespace
  {
    int pptrf(int n, double* A)
    {
      return LAPACKE_dpptrf(LAPACK_COL_MAJOR,'L',n,A);
    }

    int pptrf(int n, float* A)
    {
      return LAPACKE_spptrf(LAPACK_COL_MAJOR,'L',n,A);
    }

    int pptrs(int n, int nrhs, double const* L, double* b)
    {
      return LAPACKE_dpptrs(LAPACK_COL_MAJOR,'L',n,nrhs,L,b,n);
    }

    int pptrs(int n, int nrhs, float const* L, float* b)
    {
      return LAPACKE_spptrs(LAPACK_COL_MAJOR,'L',n,nrhs,L,b,n);
    }

    int tptrs(char trans, int n, int nrhs, double const* L, double* b)
    {
      return LAPACKE_dtptrs(LAPACK_COL_MAJOR,'L',trans,'N',n,nrhs,L,b,n);
    }

    int tptrs(char trans, int n, int nrhs, float const* L, float* b)
    {
      int info = LAPACKE_stptrs(LAPACK_COL_MAJOR,'L',trans,'N',n,nrhs,L,b,n);
//       if (true || info!=0)
//       {
//         std::cout << "n = " << n << "\n";
//         std::cout << "ldb = " << n << "\n";
//         std::cout << "nrhs = " << nrhs << "\n";
//         std::cout << "info = " << info << "\n";
//       }
      return info;
    }

  }

  // ----------------------------------------------------------------------------------------------


  namespace DomainDecompositionPreconditionerDetail
  {
    template <class Scalar>
    std::vector<Scalar> choleskyFactor(DynamicMatrix<Scalar> const& A)
    {
      assert(A.N()==A.M());

      std::vector<Scalar> l;
      l.reserve(A.N()*(A.N()+1)/2);

      for (int j=0; j<A.M(); ++j)     // store A in packed storage
        for (int i=j; i<A.N(); ++i)
          l.push_back(A[i][j]);

      int info = pptrf(A.N(),&l[0]);  // factorize
      assert(info >= 0);
      if (info>0)
        throw NonpositiveMatrixException("Nonpositive matrix for Cholesky factorization.",__FILE__,__LINE__);

      return l;
    }

    // explicit instantiation
    template std::vector<double> choleskyFactor(DynamicMatrix<double> const& A);
    template std::vector<float> choleskyFactor(DynamicMatrix<float> const& A);

    // --------------------------------------------------------------------------------------------

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
    template <class Scalar>
    void choleskySolve(int n, int nrhs, std::vector<Scalar> const& L, Scalar* b)
    {
      static_assert(std::is_floating_point<Scalar>::value,"choleskySolve only defined for floating point entries");
      int info = pptrs(n,nrhs,&L[0],b);
      assert(info>=0);
    }
#pragma GCC diagnostic pop

    // explicit instantiation
    template void choleskySolve(int n, int nrhs, std::vector<double> const& L, double* b);
    template void choleskySolve(int n, int nrhs, std::vector<float> const& L, float* b);

    template <class Scalar>
    void choleskySolve(std::vector<Scalar> const& L, DynamicMatrix<Scalar>& b)
    {
      assert(L.size()==b.N()*(b.N()+1)/2);
      choleskySolve(b.N(),b.M(),L,&b[0][0]);
    }

    // explicit instantiation
    template void choleskySolve(std::vector<double> const& L, DynamicMatrix<double>& b);
    template void choleskySolve(std::vector<float> const& L,  DynamicMatrix<float>& b);

    template <class Scalar>
    void triSolve(int n, int nrhs, char trans, std::vector<Scalar> const& L, Scalar* b)
    {
      static_assert(std::is_floating_point<Scalar>::value,"choleskySolve only defined for floating point entries");
      int info = tptrs(trans,n,nrhs,&L[0],b);
      if (info!=0)
      {
        std::cout << "info = " << info << "\n";
      }
      assert(info==0);
    }

    template <class Scalar>
    void triSolve(char trans, std::vector<Scalar> const& L,  DynamicMatrix<Scalar>& b)
    {
      assert(b.N()>0);
      assert(L.size()==b.N()*(b.N()+1)/2);
      triSolve(b.N(),b.M(),trans,L,&b[0][0]);
    }

    // explicit instantiation
    template void triSolve(char trans, std::vector<double> const& L,  DynamicMatrix<double>& b);
    template void triSolve(char trans, std::vector<float> const& L,  DynamicMatrix<float>& b);
  }

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  std::pair<size_t,size_t> farVertex(size_t const v, std::vector<std::pair<size_t,size_t>> const& edges,
                                     std::vector<size_t> const& edgeStart)
  {
    using namespace std;

    // We use a simplified version of Dijkstra, specialized for identical edge weights.

    // Compute the number of vertices.
    size_t const n = edgeStart.size()-1;
    assert(v<n);

    // initialize distances
    vector<size_t> dist(n,numeric_limits<size_t>::max());
    dist[v] = 0;

    // initialize vertices with known distance (sorted by distance)
    queue<size_t> fringe;       // todo: consider vector as basic container?
    fringe.push(v);

    // compute distances
    while (!fringe.empty())
    {
      size_t i = fringe.front();                                                  // get the nearest vertex
      fringe.pop();                                                               // and take it out
      for (size_t p=edgeStart[i]; p<edgeStart[i+1]; ++p)                          // step through all adjacent vertices
      {
        size_t j = edges[p].second;                                               // an adjacent vertex
        if (dist[j] > dist[i]+1)                                                  // if it's far away (in our case this means
        {                                                                         // we haven't seen it before),
          dist[j] = dist[i]+1;                                                    // we've found a shorter path
          fringe.push(j);                                                         // and consider it (again)
        }
      }
    }

    // Find a far vertex with minimal degree and report both the vertex and its distance to v
    size_t far = 0;
    size_t d = 0;
    int deg = std::numeric_limits<int>::max();
    for (size_t i=0; i<n; ++i)
      if (dist[i]>=d)
      {
        int degi = edgeStart[i+1]-edgeStart[i];

        if (dist[i]>d || degi<deg)
        {
          far = i;
          d = dist[i];
          deg = degi;
        }
      }

    return make_pair(far,d);
  }

  // ----------------------------------------------------------------------------------------------

  std::array<std::vector<size_t>,3> nestedDissection(std::vector<std::pair<size_t,size_t>> const& edges,
                                                     std::vector<size_t> const& edgeStart)
  {
    assert(std::is_sorted(begin(edges),end(edges),[](auto const& p1, auto const& p2)
                                           { return p1.first<p2.first || (p1.first==p2.first && p1.second < p2.second); }));

    // Find two vertices that are far from each other and use them to initialize the fringe sets f1 and f2.
    std::vector<size_t> f1{farVertex(0,edges,edgeStart).first};
    auto f2d = farVertex(f1[0],edges,edgeStart);
    std::vector<size_t> f2{f2d.first};

    // Compute number of vertices.
    size_t const n = edgeStart.size()-1;

    // If we don't find two vertices with a distance of at least 2, there is no proper partitioning.
    // We report empty sets 1 and 2 and a separator that consists of all vertices.
    if (f2d.second <= 1)
    {
      std::array<std::vector<size_t>,3> result;
      result[2].resize(n);
      std::iota(begin(result[2]),end(result[2]),0);
      return result;
    }

    // Compute the partitioning. This realizes a two-sided King algorithm.
    // We start with empty v1,v2 partitions and the two vertices far from each other
    // in the corresponding fringes.
    std::vector<char> v1Flag(n,false), v2Flag(n,false), sFlag(n,false), f1Flag(n,false), f2Flag(n,false);
    f1Flag[f1[0]] = true;
    f2Flag[f2[0]] = true;

    std::vector<size_t> tmpSep;                                                   // defined here to prevent frequent reallocations
    while (! (f1.empty() && f2.empty()))
    {
      std::swap(v1Flag,v2Flag); // treat both vertex sets alternatively in order to
      std::swap(f1,f2);         // get vertex groups of similar size
      std::swap(f1Flag,f2Flag);

      // Move all vertices from the fringe that are adjacent to both vertex sets to the separator.
      auto it = std::partition(begin(f1),end(f1),[&](size_t v)                    // partition f1 according to adjacency
      {                                                                           // such that vertices adjacent to V2 end up last.
        for (size_t i=edgeStart[v]; i!=edgeStart[v+1]; ++i)                       // Find all edges of vertex v in f1
          if (v2Flag[edges[i].second])                                            // if the edge goes into V2, then v is
            return false;                                                         // adjacent to V2 (and V1 by definition)
        return true;                                                              // - otherwise, it's a normal fringe vertex
      });
      std::for_each(it,end(f1),[&](size_t v) { sFlag[v] = true; f1Flag[v] = false; });  // move the bi-adjacent vertices into the separator
      f1.erase(it,end(f1));                                                       // and remove them from the fringe

      if (f1.empty())         // no candidate vertex available for this vertex set
        continue;             // -> switch to other set

      // Find the fringe vertex that has least edges outside V1/S and insert it into V1.
      int count = std::numeric_limits<int>::max();
      auto vIt = begin(f1);
      for (auto it=vIt; it!=end(f1); ++it)                                        // For each vertex w in f1
      {
        int c = std::count_if(begin(edges)+edgeStart[*it],                        // count the edges
                              begin(edges)+edgeStart[*it+1],                      // that do not go into V1 or Separator
                              [&](auto const& e) { return !(v1Flag[e.second] || sFlag[e.second]); });
        if (c < count)                                                            // and remember the vertex
        {                                                                         // with the least outside edges
          count = c;
          vIt = it;
        }
      }

      size_t v = *vIt;
      v1Flag[v] = true;                                                           // enter it into V1
      f1Flag[v] = false;
      std::swap(*vIt,f1.back());                                                  // and remove it from the fringe
      f1.pop_back();                                                              // (here with O(1) complexity)

      // Enter all neighbours of v into the fringe (unless they are already known to be in the separator or V1)
      for (size_t i=edgeStart[v]; i!=edgeStart[v+1]; ++i)
      {
        size_t w = edges[i].second;
        if ( !sFlag[w] && !v1Flag[w] && !f1Flag[w])
        {
          f1.push_back(w);
          f1Flag[w] = true;
        }
      }
    }

    // Extract the partitions. Vertices that could not be reached (e.g., if the graph is not
    // connected or if the fringes collide at a bottleneck behind which there are more vertices)
    // end up in v2.
    std::vector<size_t> v1, v2, separator;
    for (size_t v=0; v<n; ++v)
      if (sFlag[v])        separator.push_back(v);
      else if (v1Flag[v])  v1.push_back(v);
      else                 v2.push_back(v);

    return std::array<std::vector<size_t>,3>{v1,v2,separator};
  }

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  using namespace DomainDecompositionPreconditionerDetail;

  template <class StorageScalar>
  template <class SparseMatrix>
  NestedDissection<StorageScalar>::NestedDissection(SparseMatrix const& A, size_t skipEntries, size_t skipDimension)
  : sz(1)
  {
bool const report = NumaThreadPool::instance().isSequential();
Timings& timer=Timings::instance();

    // Get the size of matrix entries
    constexpr int m = SparseMatrix::block_type::rows;

    // Nested dissection is only beneficial if the overhead it induces pays off in reduction of
    // Cholesky factor entries. If the size of A is too small, we do not dissect (i.e. everything
    // is the separator), and perform dense factorization. Note that this can also result if we
    // fail to find a proper partitioning.
    if (A.N() >= skipDimension)
    {
if (report)  timer.start("edge generation");
      // Compute adjacency graph as edge list
      std::vector<std::pair<size_t,size_t>> edges;
      std::vector<size_t> edgeStart{0};

      for (size_t r=0; r<A.N(); ++r)
      {
        for (auto ci=A[r].begin(); ci!=A[r].end(); ++ci)
          if (ci.index() != r)                                  // omit diagonal
            edges.push_back(std::make_pair(r,ci.index()));
        edgeStart.push_back(edges.size());
      }
if (report)  timer.stop("edge generation");

if (report)  timer.start("dissection");
      // Compute row and column partitioning
      auto idx = nestedDissection(edges,edgeStart);
      v1 = std::move(idx[0]);
      v2 = std::move(idx[1]);
      s = std::move(idx[2]);

      // If the zero block generated by dissection is not large enough
      // to reward the overhead of dissection (non-local memory access,
      // small matrices), we omit dissection and perform a dense factorization.
      if (v1.size()*v2.size()*m*m < skipEntries)
      {
        v1.clear();
        v2.clear();
        s.resize(A.N());
        std::iota(begin(s),end(s),0);
      }
if (report)  timer.stop("dissection");
    }
    else
    {
      s.resize(A.N());
      std::iota(begin(s),end(s),0);
    }
    // ND Cholesky factor looks like
    // [L1     ]
    // [   L2  ]
    // [R1 R2 S]

    // In case the index sets 1 and 2 are empty, we just compute S.
    if (v1.empty())
    {
      assert(v2.empty());
if (report)  timer.start("extraction dense");
      DynamicMatrix<StorageScalar> ss = DynamicMatrixDetail::flatMatrix(full(A,s,s));
if (report)  timer.stop("extraction dense");

if (report)  timer.start("S factor");
      auto ls = choleskyFactor(ss);
      sFactor.reserve(ls.size());
      std::copy(begin(ls),end(ls),std::back_inserter(sFactor));
if (report)  timer.stop("S factor");
      return;
    }

    assert(!v2.empty());


    // Extract the relevant blocks
if (report)  timer.start("extraction sparse");
    auto a1 = submatrix<SparseMatrix>(A,v1,v1);
    auto a2 = submatrix<SparseMatrix>(A,v2,v2);
if (report)  timer.stop("extraction sparse");
if (report)  timer.start("extraction dense");
    DynamicMatrix<StorageScalar> ss = DynamicMatrixDetail::flatMatrix(full(A,s,s));
    DynamicMatrix<StorageScalar> r1t = DynamicMatrixDetail::flatMatrix(full(A,v1,s));
    DynamicMatrix<StorageScalar> r2t = DynamicMatrixDetail::flatMatrix(full(A,v2,s));
if (report)  timer.stop("extraction dense");

    // Factorize the diagonal blocks L1*L1' = A1
    l1FactorND.reset(new NestedDissection<StorageScalar>(a1,skipEntries,skipDimension));
    l2FactorND.reset(new NestedDissection<StorageScalar>(a2,skipEntries,skipDimension));
    sz += l1FactorND->size() + l2FactorND->size();


    // Compute subdiagonal blocks L1*R1' = A1s
if (report)  timer.start("L1/L2 blocksolve");
    l1FactorND->template trisolve<m>(r1t,'N');
    r1 = transpose(r1t);
    l2FactorND->template trisolve<m>(r2t,'N');
    r2 = transpose(r2t);
if (report)  timer.stop("L1/L2 blocksolve");

    // Compute and factorize separator block SS' = As - R1*R1' - R2*R2'
    if (! s.empty())
    {
if (report)  timer.start("S build");
    ss -= transpose(r1t)*r1t + transpose(r2t)*r2t;
if (report)  timer.stop("S build");
if (report)  timer.start("S factor");
      auto ls = choleskyFactor(ss);
      sFactor.reserve(ls.size());
      std::copy(begin(ls),end(ls),std::back_inserter(sFactor));
if (report)  timer.stop("S factor");
    }
  }

  template <class StorageScalar>
  template <class FactorizationScalar>
  NestedDissection<StorageScalar>& NestedDissection<StorageScalar>::operator=(NestedDissection<FactorizationScalar>&& nd)
  {
    v1 = std::move(nd.v1);
    v2 = std::move(nd.v2);
    s = std::move(nd.s);

    r1 = nd.r1;
    r2 = nd.r2;
    sFactor.reserve(nd.sFactor.size());
    std::copy(begin(nd.sFactor),end(nd.sFactor),back_inserter(sFactor));

    if (nd.l1FactorND.get())
    {
      assert(nd.l2FactorND.get());
      l1FactorND.reset(new NestedDissection<StorageScalar>());
      *l1FactorND = std::move(*nd.l1FactorND);
      l2FactorND.reset(new NestedDissection<StorageScalar>());
      *l2FactorND = std::move(*nd.l2FactorND);
    }

    sz = nd.sz;
    return *this;
  }


// USEFUL FOR DEBUGGING OUTPUT
template <class Scalar>
DynamicMatrix<Scalar> matrixOf(std::vector<Scalar> const& l)
{
  int n = (int) ((std::sqrt(1.0+8*l.size())-1)/2 + 0.4);
  assert(n*(n+1)/2 == l.size());
  DynamicMatrix<Scalar> L(n,n);
  auto it = begin(l);
  for (int j=0; j<n; ++j)
    for (int i=j; i<n; ++i)
      L[i][j] = *it++;
  return L;
}

  template <class StorageScalar>
  template <class Vector>
  void NestedDissection<StorageScalar>::solve(Vector& y) const
  {
    constexpr int m = Vector::value_type::dimension;
    DynamicMatrix<StorageScalar> ys(y.N()*m,1);

    for (size_t i=0; i<y.N(); ++i)
      for (int k=0; k<m; ++k)
        ys[m*i+k][0] = y[i][k];

    trisolve<m>(ys,'B');

    for (size_t i=0; i<y.N(); ++i)
      for (int k=0; k<m; ++k)
        y[i][k] = ys[m*i+k][0];
  }


  template <class StorageScalar>
  template <int m>
  void NestedDissection<StorageScalar>::trisolve(DynamicMatrix<StorageScalar>& y, char trans) const
  {
    trans = std::toupper(trans);

    auto gather = [](auto const& from, auto const& idx)
    {
      DynamicMatrix<StorageScalar> to(m*idx.size(),from.M());
      for (size_t j=0; j<from.M(); ++j)
        for (size_t i=0; i<idx.size(); ++i)
          for (int k=0; k<m; ++k)
            to[i*m+k][j] = from[m*idx[i]+k][j];
      return to;
    };

    DynamicMatrix<StorageScalar> y1 = gather(y,v1), y2 = gather(y,v2), ys = gather(y,s);

    if (trans=='N' || trans=='B' )
    {
      if (v1.size()>0)
      {
        l1FactorND->template trisolve<m>(y1,'N');          // y2 <- L2 \ y2
//         ys -= r1*y1;                                       // ys <- ys - R1*y1
        DynamicMatrixDetail::gemm(false,false,ys.N(),ys.M(),y1.N(),-1.0,&r1[0][0],r1.N(),&y1[0][0],y1.N(),1.0,&ys[0][0],ys.N());
        l2FactorND->template trisolve<m>(y2,'N');          // y2 <- L2 \ y2
//         ys -= r2*y2;                                       // ys <- ys - R2*y2
        DynamicMatrixDetail::gemm(false,false,ys.N(),ys.M(),y2.N(),-1.0,&r2[0][0],r2.N(),&y2[0][0],y2.N(),1.0,&ys[0][0],ys.N());
      }
      if (! s.empty())
        triSolve('N',sFactor,ys);                          // ys <- S \ ys
    }

    if (trans=='C' || trans=='B')
    {
      if (! s.empty())
        triSolve('C',sFactor,ys);                          // ys <- S' \ ys
      if (v2.size()>0)
      {
//         y2 -= transpose(r2)*ys;                            // y2 <- y2 - R2'*ys
// The above line creates *two* temproraries, while the direct BLAS level 3 call to GEMM creates none.
// The difference is negligible for factorization, but during solution, where the compute load is
// relatively smaller, this makes a difference.
        DynamicMatrixDetail::gemm(true,false,y2.N(),y2.M(),ys.N(),-1.0,&r2[0][0],r2.N(),&ys[0][0],ys.N(),1.0,&y2[0][0],y2.N());

        l2FactorND->template trisolve<m>(y2,'C');           // y2 <- L2' \ y2

//         y1 -= transpose(r1)*ys;                             // y1 <- y1 - R1'*ys
        DynamicMatrixDetail::gemm(true,false,y1.N(),y1.M(),ys.N(),-1.0,&r1[0][0],r1.N(),&ys[0][0],ys.N(),1.0,&y1[0][0],y1.N());

        l1FactorND->template trisolve<m>(y1,'C');           // y1 <- L1' \ y1
      }
    }

    auto scatter = [](auto const& from, auto const& idx, auto& to)
    {
      for (size_t j=0; j<from.M(); ++j)
        for (size_t i=0; i<idx.size(); ++i)
          for (int k=0; k<m; ++k)
            to[m*idx[i]+k][j] = from[i*m+k][j];
    };
    scatter(y1,v1,y);
    scatter(y2,v2,y);
    scatter(ys,s,y);
  }

  template <class StorageScalar>
  size_t NestedDissection<StorageScalar>::storageSize() const
  {
    size_t result = sizeof(StorageScalar)*sFactor.size() + sizeof(size_t)*(v1.size()+v2.size()+s.size());
    if (!v1.empty())
      result += l1FactorND->storageSize() + l2FactorND->storageSize() + sizeof(StorageScalar)*(r1.N()*r1.M() + r2.N()*r2.M());
    return result + sizeof(NestedDissection<StorageScalar>);
  }

  template <class StorageScalar>
  size_t NestedDissection<StorageScalar>::flops() const
  {
    size_t result = 2*sFactor.size();
    if (!v1.empty())
      result += l1FactorND->flops()+l2FactorND->flops()+2*(r1.N()*r1.M()+r2.N()*r2.M());
    return result;
  }

  // explicit instantiation
  template class NestedDissection<double>;
  template class NestedDissection<float>;
  template NestedDissection<double>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> const& A, size_t, size_t);
  template NestedDissection<double>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<double,2,2>> const& A, size_t, size_t);
  template NestedDissection<double>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<double,3,3>> const& A, size_t, size_t);
  template NestedDissection<float>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> const& A, size_t, size_t);
  template NestedDissection<float>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<double,2,2>> const& A, size_t, size_t);
  template NestedDissection<float>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<double,3,3>> const& A, size_t, size_t);
  template NestedDissection<float>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<float,1,1>> const& A, size_t, size_t);
  template NestedDissection<float>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<float,2,2>> const& A, size_t, size_t);
  template NestedDissection<float>::NestedDissection(Dune::BCRSMatrix<Dune::FieldMatrix<float,3,3>> const& A, size_t, size_t);
  template NestedDissection<double>& NestedDissection<double>::operator=(NestedDissection<double>&& nd);
  template NestedDissection<float>& NestedDissection<float>::operator=(NestedDissection<double>&& nd);
  template NestedDissection<float>& NestedDissection<float>::operator=(NestedDissection<float>&& nd);
  template void NestedDissection<double>::solve(Dune::BlockVector<Dune::FieldVector<double,1>>& y) const;
  template void NestedDissection<double>::solve(Dune::BlockVector<Dune::FieldVector<double,2>>& y) const;
  template void NestedDissection<double>::solve(Dune::BlockVector<Dune::FieldVector<double,3>>& y) const;
  template void NestedDissection<float>::solve(Dune::BlockVector<Dune::FieldVector<float,1>>& y) const;
  template void NestedDissection<float>::solve(Dune::BlockVector<Dune::FieldVector<float,2>>& y) const;
  template void NestedDissection<float>::solve(Dune::BlockVector<Dune::FieldVector<float,3>>& y) const;
}

#ifdef UNITTEST

int main(void)
{
  using namespace Kaskade;

  using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  Matrix A(5,5,5,0.5,Matrix::implicit);
  A.entry(0,0) = 1;                       A.entry(2,0) = -.1;                     A.entry(4,0) = -.1;
                      A.entry(1,1) = 1;   A.entry(2,1) = -.1;                     A.entry(4,1) = -.1;
  A.entry(0,2) = -.1; A.entry(1,2) = -.1; A.entry(2,2) = 1;                       A.entry(4,2) = -.1;
                                                              A.entry(3,3) = 1;   A.entry(4,3) = -.1;
  A.entry(0,4) = -.1; A.entry(1,4) = -.1; A.entry(2,4) = -.1; A.entry(3,4) = -.1; A.entry(4,4) = 1;
  A.compress();

  std::cout << A << "\n";
  std::cout << full(A) << "\n";

  NestedDissection<double> nd(A, (double)1.0);
  Dune::BlockVector<Dune::FieldVector<double,1>> r(5);
  for (int i=0; i<A.N(); ++i)
    r[i] = i+1;
  nd.solve(r);
  std::cout << r << "\n";

  return 0;
}

#endif
