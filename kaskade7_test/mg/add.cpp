/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2022-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "metis.h"

#include "mg/add.hh"


namespace Kaskade::BDDC
{
  template <class Index>
  std::vector<std::vector<Index>> algebraicDomainDecomposition(NumaCRSPattern<Index> const& A, int const n)
  {
    // We decompose the indices (i.e. dofs) using Metis. For that, we create a "mesh" data 
    // structure where each nonzero matrix entry (i,j) is an "element" connecting the 
    // "mesh nodes" i and j. We then partition the mesh elements into different groups.
    // Finally, the subdomain corresponding to such a group consists of all indices that
    // are incident to some element in the group.
    
    // First compute the number of "elements". We can restrict the attention to the lower
    // left triangular part of the matrix. We also store the row/column indices of 
    // relevant matrix entries for later use.
    std::vector<idx_t> eind;
    for (int i=0; i<A.nodes(); ++i)                               // step through all the NUMA chunks
    {                                                             // TODO: run this in parallel
      auto const& chunk = *A.pattern(i);                          // access the chunk
      
      auto rStart = chunk.first();
      for (size_t r=rStart; r<chunk.last(); ++r)                  // step through all rows
        for (auto cit=chunk.colStartIterator(r-rStart);           // and extract the entries
             cit!=chunk.colStartIterator(r-rStart+1) && *cit<r;   // in the lower left triangle
             ++cit)
        {
          eind.push_back(*cit);
          eind.push_back(r);
        }
    }
    idx_t m = eind.size()/2;                                      // number of "elements" in the "mesh"
    idx_t k = A.N();                                              // number of "nodes" in the "mesh"
    
    // Create the "mesh" data structure. This consists of a list of nodes for each "element".
    // In our case, all "elements" are just edges incident to exactly two nodes.
    std::vector<idx_t> eptr(m+1);           // start index of node list of element i in eind
    for (Index i=0; i<=m; ++i)              // this is trivial, as the node lists all have 
      eptr[i] = 2*i;                        // length 2
      
    // Now we observe that eind already contains the concatenated node lists.
    
    idx_t ncommon = 1;
    idx_t nparts = n;
    idx_t objval;
    std::vector<idx_t> epart(m), npart(k);
    int retcode = METIS_PartMeshDual(&m,&k,&eptr[0],&eind[0],
                                     nullptr,nullptr,&ncommon,&nparts,
                                     nullptr,nullptr,&objval,&epart[0],&npart[0]);
    
    if (retcode != METIS_OK)
      std::cerr << "Metis failed: " << retcode << "\n";
    
    // Now extract the partitioning.
    std::vector<std::vector<Index>> subdoms(n);
    for (size_t i=0; i<m; ++i)                          // First include the vertices (dofs) into the 
    {                                                   // partition of their incident edges
      assert(epart[i]<n);
      subdoms[epart[i]].push_back(eind[eptr[i]]);
      subdoms[epart[i]].push_back(eind[eptr[i]+1]);
    }
    for (int i=0; i<n; ++i)                             // Then remove duplicate dofs from the partitions.
    {
      auto& s = subdoms[i];
      std::sort(begin(s),end(s));
      s.erase(std::unique(begin(s),end(s)),end(s));
    }
    
    return subdoms;
  }
  
  // Explicit instantiation
  template std::vector<std::vector<int>> algebraicDomainDecomposition(NumaCRSPattern<int> const&, int);
  template std::vector<std::vector<size_t>> algebraicDomainDecomposition(NumaCRSPattern<size_t> const&, int);

  // ----------------------------------------------------------------------------------------------

  template <class Scalar>
  Scalar equal(Scalar a, Scalar b)
  {
    return std::abs(a-b) / ((std::abs(a)+std::abs(b))*std::numeric_limits<Scalar>::epsilon());
  }

  template <class Matrix>
  double checkSymmetry(Matrix const& A)
  {
    double maxRelError = 0;
    for (auto ri=A.begin(); ri!=A.end(); ++ri)
      for (auto ci=ri.begin(); ci!=ri.end(); ++ci)
        maxRelError = std::max(maxRelError,equal((double)*ci,(double)(A[ri.index()][ci.index()])));

    return maxRelError;
  }

  template <class Matrix, class Index>
  bool isSum(Matrix A, std::vector<Matrix> const& As, std::vector<std::vector<Index>> const& globIdx)
  {
    double maxEntry = 0 ;
    for (auto ri=A.begin(); ri!=A.end(); ++ri)
      for (auto ci=ri.begin(); ci!=ri.end(); ++ci)
        maxEntry = std::max(maxEntry,std::abs((double)*ci));

    for (int k=0; k<As.size(); ++k)
    {
      auto const& Ak = As[k];
      for (auto ri=Ak.begin(); ri!=Ak.end(); ++ri)
        for (auto ci=ri.begin(); ci!=ri.end(); ++ci)
          A[globIdx[k][ri.index()]][globIdx[k][ci.index()]] -= *ci;
    }
    double maxEntry2 = 0 ;
    for (auto ri=A.begin(); ri!=A.end(); ++ri)
      for (auto ci=ri.begin(); ci!=ri.end(); ++ci)
        maxEntry2 = std::max(maxEntry2,std::abs((double)*ci));

    std::cout << "sum deviation: " << maxEntry2 << " vs " << maxEntry << "\n";
    return maxEntry2 < 1e-13*maxEntry;
  }


  template <class Scalar, class Index>
  std::vector<NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Index>>
  algebraicMatrixDecomposition(NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Index> const& A,
                               std::vector<std::vector<Index>> const& subIndices,
                               Scalar tolerance, int maxIter)
  {
    using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Index>;
    ScopedTimingSection tsec("algebraicMatrixDecomposition");

    assert(A.N() == A.M());

    int const K = subIndices.size();

std::cout.precision(12);
// std::cout << "A = [\n" << full(A) << "];\n";
    if (Scalar s = checkSymmetry(A); s > 100)
      std::cout << "A is not symmetric in " << __LINE__ << " by relative deviation " << s << " eps\n";
auto Asave = A;

    // Get the subdomains in which any global dof is contained. Note that by construction
    // the lists of subdomains are sorted.
    std::vector<std::vector<Index>> subdomains(A.N());
    for (int s=0; s<K; ++s)
    {
      for (auto i: subIndices[s])
        subdomains[i].push_back(s);
    }

    std::vector<Matrix> As;
    std::vector<std::vector<std::tuple<Index,Index,Index,Index,int>>> bData(K);

    tsec.timer().start("initial decomposition");
    for (int k=0; k<K; ++k)
    {
      // Create submatrix.
      Matrix Ak(subIndices[k],A);
      Ak = 0;
      // We first transfer all the entries which are interior to a subdomain, i.e.
      // they belong only to *one* subdomain. These can be copied directly.
      for (auto ri=Ak.begin(); ri!=Ak.end(); ++ri)
        for (auto ci=ri.begin(); ci!=ri.end(); ++ci)
        {
          Index rglob = subIndices[k][ri.index()];
          Index cglob = subIndices[k][ci.index()];
          std::vector<Index> intersection;
          std::set_intersection(begin(subdomains[rglob]),end(subdomains[rglob]),
                                begin(subdomains[cglob]),end(subdomains[cglob]),
                                back_inserter(intersection));
          int nSubs = intersection.size();
          assert(nSubs >= 1);

          if (nSubs == 1)
          {
            *ci = A[rglob][cglob];
            A[rglob][cglob] = 0;
          }
          else
            bData[k].push_back({ri.index(),ci.index(),rglob,cglob,nSubs});
        }
      As.push_back(Ak);
//       if (Scalar s = checkSymmetry(A); s > 100)
//       {
//         std::cout << "k=" << k << ": submatrix not symmetric in " << __LINE__ << "\n";
//         std::cout << "A" << k << " = [\n" << full(Ak) << "];\n";
//       }
    }
    tsec.timer().stop("initial decomposition");

    // Now move as much of the diagonal entries over to the submatrices as is absolutely necessary
    // to have them weakly diagonally dominant.
    {
      ScopedTimingSection tsec("diagonal distribution");
      for (int k=0; k<K; ++k)
      {
        auto& Ak = As[k];
        for (auto& [rloc,cloc,rglob,cglob,ns]: bData[k])
          if (rloc==cloc)
          {
            Scalar rowSum = 0;
            for (auto ci=Ak[rloc].begin(); ci!=Ak[rloc].end(); ++ci)
              if (ci.index() != cloc)
                rowSum += std::abs(*ci);
              Ak[rloc][cloc] = rowSum;
            A[rglob][cglob] -= rowSum;
            assert(A[rglob][cglob] >= 0);
          }
          //       if (Scalar s = checkSymmetry(A); s > 100)
          //         std::cout << "k=" << k << ": submatrix not symmetric in " << __LINE__ << "\n";
      }
    }

    // Now we have to distribute the remaining parts of A to the Ak's, such that all Ak's
    // remain diagonally dominant. Thus, we need to find an admissible solution of the
    // following linear inequality system, where k denotes the subdomain, i row index, and
    // j column index, with some sloppy notation not distinguishing between local and global
    // indices):
    // 1) for all i,j: \sum_k a_{ij}^k = a_{ij}       -> submatrices sum up to the global matrix
    // 2) for all k,i: \sum_j |a_{ij}^k| \le a_{ii}^k -> submatrices are weakly diagonally dominant
    // 3) for all k,i,j: a_{ij}^k = a_{ji}^k          -> symmetry
    // This could be formulated as a linear program to be solved, but even though this problem is
    // defined only on the interface matrix entries, it is of considerable size. We thus try a
    // quick heuristic here instead of calling a general purpose LP solver.

    // We start by distributing the entries of A equally, i.e. just divided by the number of
    // subdomains containing the entry. This satisfies equations 1) and 3), but may violate 2).
    // In a Gauss-Seidel (or rather Southwell) manner, we look for (k,i) violating 2), and aim at
    // redistributing entries between submatrices such that 2) is locally satisfied and 1) and 3)
    // are retained. This will (most likely) converge to the desired decomposition, though no
    // convergence proof is available.

    for (int k=0; k<K; ++k)
    {
      auto& Ak = As[k];
      auto& b = bData[k];
      for (auto& [rloc,cloc,rglob,cglob,ns]: b)
        Ak[rloc][cloc] += A[rglob][cglob] / ns;

      // For the remaining operations, we need no longer the bData entries for off-diagonal
      // entries. So we remove them here once.
      auto last = std::remove_if(begin(b),end(b),
                                 [](auto const& bd) { return std::get<0>(bd) != std::get<1>(bd); });
      b.erase(last,end(b));
    }

//     isSum(Asave,As,subIndices);

    // For a global dof index, get its local index in the given subdomain k.
    // Precondition: the dof shall be contained in the subdomain.
    auto localIndex = [&](int k, Index globalIndex)
    {
      auto const& idxk = subIndices[k];
      auto it = std::find(idxk.begin(),idxk.end(),globalIndex);
      assert(it != idxk.end());
      return it-idxk.begin();
    };

    // Add the value a to all submatrix entries corresponding to the global dof given by
    // the indices (rglob,cglob), where the subdomain numbers are stored in subdoms. Omit
    // the subdomain with number k.
    auto addToOtherSubmatrices = [&](auto const& subdoms, int k, Index rglob, Index cglob, Scalar a)
    {
      for (int khat: subdoms)
        if (khat != k)
        {
          Index rhat = localIndex(khat,rglob);
          Index chat = localIndex(khat,cglob);
          As[khat][rhat][chat] += a;
          if (rhat != chat)
            As[khat][chat][rhat] += a;            // retain symmetry!
        }
    };

    {
      ScopedTimingSection tsec("GS");

      for (int gsiter=0; gsiter<600; ++gsiter)
      {
        Scalar maximumRelativeDefect = 0;
        std::cout << "--- GS iter " << gsiter << " ----\n";
        for (int k=0; k<K; ++k)
        {
          auto& Ak = As[k];
          for (auto& [rloc,cloc,rglob,_dummy1,_dummy2]: bData[k])
          {
            assert(rloc==cloc);                                          // Remember we erased all off-diagonal entries.
            // Now compute the defect \sum_j |a_{ij}^k| - a_{ii}^k
            // which shall be nonpositive for weakly diagonally
            // dominant matrices. Also compute the relative defect
            Scalar defect = 0;
            Scalar magnitude;
            for (auto ci=Ak[rloc].begin(); ci!=Ak[rloc].end(); ++ci)
              if (ci.index() == cloc)
              {
                defect -= *ci;
                magnitude = *ci;
              }
              else
                defect += std::abs(*ci);

              Scalar const relativeDefect = defect/magnitude;

            if (relativeDefect > 10*std::numeric_limits<Scalar>::epsilon())
            {
              maximumRelativeDefect = std::max(maximumRelativeDefect,relativeDefect);
              //             std::cout << "k=" << k << " i=" << rloc << " defect=" << defect << " magnitude=" << magnitude << "\n";

              int nEntries = Ak[rloc].size();
              Scalar diagonalCorrection = defect/2; // defect/nEntries;
              Scalar offDiagonalCorrection = defect/2/(nEntries-1);

              for (auto ci=Ak[rloc].begin(); ci!=Ak[rloc].end(); ++ci)
                if (ci.index() == rloc)
                {
                  *ci += diagonalCorrection;
                  addToOtherSubmatrices(subdomains[rglob],k,rglob,rglob,
                                        -diagonalCorrection/(subdomains[rglob].size()-1));
                }
                else
                {
                  Index cglob = subIndices[k][ci.index()];
                  std::vector<Index> intersection;
                  std::set_intersection(begin(subdomains[rglob]),end(subdomains[rglob]),
                                        begin(subdomains[cglob]),end(subdomains[cglob]),
                                        back_inserter(intersection));
                  if (intersection.size() > 1)
                  {
                    Scalar sign = *ci>0? 1: -1;
                    *ci -= sign*offDiagonalCorrection;
                    Ak[ci.index()][rloc] -= sign*offDiagonalCorrection;
                    addToOtherSubmatrices(intersection,k,rglob,cglob,
                                          sign*offDiagonalCorrection/(intersection.size()-1));
                  }
                }
            }
          }
          //         if (Scalar s = checkSymmetry(A); s > 100)
          //           std::cout << "k=" << k << ": submatrix not symmetric in " << __LINE__ << "\n";

        }
        std::cout << "max rel defect: " << maximumRelativeDefect << "\n";
        if (maximumRelativeDefect < tolerance)
          break;

        // isSum(Asave,As,subIndices);
      }
    }

//     for (int k=0; k<K; ++k)
//       std::cout << "A" << k+1 << " = [ \n" << full(As[k]) << "];\n";

    return As;

  }

  template std::vector<NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>,size_t>>
  algebraicMatrixDecomposition(NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>,size_t> const& A,
                               std::vector<std::vector<size_t>> const& subIndices,
                               double tolerance, int maxIter);
  
}
