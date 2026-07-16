 
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

#ifndef BDDC_HPP
#define BDDC_HPP

#include <algorithm>

#include "linalg/direct.hh"
#include "mg/bddc.hh"
#include "utilities/timing.hh"

namespace Kaskade::BDDC
{


  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <int m, class Index>
  InterfaceAverages<m,Index>::InterfaceAverages(std::vector<std::vector<LocalDof>> const& sharedDofs,
                                                std::vector<Index> const& subdomainSize_,
                                                int const types)
  : ifs(subdomainSize_.size())
  , subdomainSize(subdomainSize_)
  , cifs(subdomainSize.size())
  {
    for (auto const& localDofs: sharedDofs)
      for (auto [s1,n1]: localDofs)
        for (auto [s2,n2]: localDofs)
          if (s1 != s2)
          {
            auto i = std::find_if(ifs[s1].begin(),ifs[s1].end(),
                                  [=](auto const& neighbor) { return neighbor.first==s2; });
            if (i==ifs[s1].end())
              i = ifs[s1].insert(i,{s2,{}});
            i->second.push_back(n1);
          }

    for (auto& interface: ifs)                                    // Sort the interface definition of
      std::sort(interface.begin(),interface.end(),FirstLess());   // each subdomain by neigbor id.
                                                                  // Note the local dofs are sorted
                                                                  // according to global ordering.


    // Find the coarse interfaces, i.e. the sets of dofs shared by the same group [s1,...,sk] of
    // subdomains. We create a map [s1,...,sk] -> [(s1,i1),(s1,j1),...,(sk,lk)] mapping the
    // subdomain groups (sorted) to the list of subdomain-local dofs included in the coarse interface.
    // Each such coarse interface defines a single constraint: the average over all dofs in this
    // interface shall agree for all incident subdomains.
    std::map<std::vector<int>,std::vector<LocalDof>> coarseInterfaces; // TODO: consider sorted vector
    for (int i=0; i<sharedDofs.size(); ++i)                       // consider all global dofs
    {
      auto const& localDofs = sharedDofs[i];
      if (localDofs.size()==1)                                    // ignore (useless) dofs contained only
        continue;                                                 // in one subdomain - that's no interface

      std::vector<int> subs;
      for (auto const& [s,_]: localDofs)                          // For each global dof, find the
        subs.push_back(s);                                        // set of affected subdomains
      std::sort(subs.begin(),subs.end());                         // and sort them,

      auto& dofs = coarseInterfaces[subs];                        // then enter them without duplicates,
      dofs.insert(dofs.end(),localDofs.begin(),localDofs.end());  // collecting all affected local dofs.
    }


    // For each subdomain, extract the list of its local dofs associated to the coarse constraints
    // in which the subdomain participates. Simultaneously, we count the dimensionality, i.e. number
    // of scalar constraints, for each coarse constraint. In our case, we aim at enforcing
    // equality of the Dirichlet interface values, which means all constraint dimensions are
    // just the number of variable components.
    for (auto const& [subs,localDofs]: coarseInterfaces)
    {
      // We distinguish three types of interfaces: (i) faces: only two incident subdomains,
      // (ii) corners: only a single dof (iii) edges: all the rest. This distinction will usually
      // but need not agree with geometric notion of faces, corners, and edges.
      InterfaceType type = subs.size() == 2?                FACE
                         : subs.size() == localDofs.size()? CORNER
                         :                                  EDGE;

      if (type & types)                                           // Create a constraint only if
      {                                                           // it is of required type.
        ccn.push_back({subs,components});                         // Register this constraint.

        for (auto s: subs)                                        // For each affected subdomain,
          cifs[s].push_back({});                                  // extract the affected local dofs.
        for (auto [s,dof]: localDofs)
          cifs[s].back().push_back(dof);
      }
    }

    // If there are two or more subdomains, they should be coupled by some coarse constraints.
    if (ccn.empty() && ifs.size()>1)
      throw UnboundedProblemException("No BDDC coarse constraints found.",__FILE__,__LINE__);
  }


  template <int m, class Index>
  Interfaces const& InterfaceAverages<m,Index>::interfaces(int id) const
  {
    assert(0<=id && id<ifs.size());
    return ifs[id];
  }


  template <int m, class Index>
  template <class Scalar>
  NumaBCRSMatrix<Dune::FieldMatrix<Scalar,m,m>> InterfaceAverages<m,Index>::coarseConstraint(int id) const
  {
    Index const nCols = subdomainSize[id];            // as many columns as subdomain has dofs
    Index const nRows = cifs[id].size();              // and rows as many as coarse constraints

    NumaCRSPatternCreator<> creator(nRows,nCols);
    for (int r=0; r<nRows; ++r)
    {
      auto const& dofs = cifs[id][r];
      for (auto c: dofs)
        creator.addElement(r,c);
    }

    NumaBCRSMatrix<Dune::FieldMatrix<Scalar,m,m>> C(creator);
    C = unitMatrix<Scalar,m>();                       // each entry is 1 - this means in each row
                                                      // (constraint) we're summing up: averaging.

    return C;
  }


  template <int m, class Index>
  CoarseConstraints InterfaceAverages<m,Index>::coarseConstraints() const
  {
    return ccn;
  }


  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <int m, class Scalar_, class CoarseScalar>
  template <class Interface>
  Subdomain<m,Scalar_,CoarseScalar>::Subdomain(int id, AMatrix const& A, XVector const& f, Interface const& ifs)
  : Subdomain(id,A,f,ifs.interfaces(id),ifs.coarseConstraint(id))
  { }


  template <int m, class Scalar_, class CoarseScalar>
  template <int cm>
  Subdomain<m,Scalar_,CoarseScalar>::Subdomain(int id_, AMatrix const& A_, XVector const& f_,
                                               std::vector<std::pair<int,std::vector<int>>> const& interfaces_,
                                               NumaBCRSMatrix<Dune::FieldMatrix<Scalar,cm,components>> const& C_)
  : id(id_), nx(A_.N()*components), nc(C_.N()*components)
  , interfaces(interfaces_)
  , f(nx)
  , u(nx)
  , r(nx)
  , du(nx)
  , avg(nx), averagingFactor(nx,1)
  , rawCorrection(nx), residual(nx), restrictedResidual(nx)
  , X(nx,nc)
  {
    bool const doTimings = false;

    // Flatten out the provided right hand side vector and compute the residual
    // r = f - A*u (which, as u=0, is just r = f).
    for (int i=0; i<A_.N(); ++i)
      for (int j=0; j<components; ++j)
        f[i*components+j] = f_[i][j];

    auto& timer = Timings::instance();

    // When exchanging boundary values with neighboring subdomains across the shared
    // interface, we will add up all contributions from the neighbors, and average by
    // dividing by the number of contributions to each entry. Here we count the number
    // of contributions for each entry.
    for (auto const& [s,ids]: interfaces)
      for (int i: ids)
        for (int j=0; j<components; ++j)
          averagingFactor[i*components+j] += 1;


    // Create B = [A C']
    //            [C 0 ]. This could be done using reshapeBlocks, horzcat and vertcat
    // very concisely, but with a lot of overhead. Instead, we create the sparsity
    // pattern directly.
    if (doTimings) timer.start("matrix creation");

    A = reshapeBlocks<1,1>(A_);
    C = reshapeBlocks<1,1>(C_);

    if (doTimings) timer.start("entering A");
    std::vector<size_t> colIndices, rowIndices;
    NumaCRSPatternCreator<> Bcreator(nx+nc,nx+nc);
    for (int row=0; row<nx; ++row)
    {
      colIndices.clear();
      for (auto ci=A_[row].begin(); ci!=A_[row].end(); ++ci)
        for (int j=0; j<components; ++j)
          colIndices.push_back(ci.index()*components+j);
      Bcreator.addElements(&row,&row+1,colIndices.begin(),colIndices.end());
    }
    if (doTimings) timer.stop("entering A");

    // remember not to enter entries alternatingly into dense row (C) and sparse rows (C^T)
    if (doTimings) timer.start("entering C");
    for (int row=nx; row<nx+nc; ++row)
    {
      colIndices.clear();
      rowIndices.clear();
      for (auto ci=C_[row-nx].begin(); ci!=C_[row-nx].end(); ++ci)
        for (int j=0; j<components; ++j)
          colIndices.push_back(ci.index()*components+j);
      for (int i=0; i<cm; ++i)
        rowIndices.push_back(row*cm+i);

      Bcreator.addElements(rowIndices.begin(),rowIndices.end(),
                           colIndices.begin(),colIndices.end());
    }
    if (doTimings) timer.stop("entering C");
    if (doTimings) timer.start("entering C^T");
    for (int row=0; row<nc; ++row)
    {
      for (auto ci=C_[row].begin(); ci!=C_[row].end(); ++ci)
        Bcreator.addDenseBlock(ci.index()*components,(ci.index()+1)*components,
                               nx+row*cm,nx+(row+1)*cm);
    }
    if (doTimings) timer.stop("entering C^T");
    if (doTimings) timer.start("creating matrix");
    auto B = BMatrix(Bcreator);
    if (doTimings) timer.stop("creating matrix");

    if (doTimings) timer.start("filling matrix");
    for (int row=0; row<nx; ++row)
      for (auto ci=A_[row].begin(); ci!=A_[row].end(); ++ci)
        for (int i=0; i<components; ++i)
          for (int j=0; j<components; ++j)
            B[row*components+i][ci.index()*components+j] = (*ci)[i][j];

    for (int row=0; row<nc; ++row)
      for (auto ci=C_[row].begin(); ci!=C_[row].end(); ++ci)
        for (int i=0; i<cm; ++i)
          for (int j=0; j<components; ++j)
          {
            B[nx+row*cm+i][ci.index()+j] = (*ci)[i][j];   // C
            B[ci.index()+j][nx+row*cm+i] = (*ci)[i][j];   // C^T
        }
    if (doTimings) timer.stop("filling matrix");
    if (doTimings) timer.stop("matrix creation");


    // TODO: MUMPS provides a partial factorization, splitting interior and boundary nodes, and gives access
    // to the Schur complement. This might be of interest for not having both Aii and B factorized.
    // But thenn, MUMPS can apparently not be called from several threads simultaneously, since its memory
    // management is not thread-safe.

    if (doTimings) timer.start("matrix factorization");
    Binv = DirectSolver<Vector,Vector>(B,DirectType::UMFPACK,MatrixProperties::SYMMETRICSTRUCTURE);
    if (doTimings) timer.stop("matrix factorization");

    // Create factorization of Aii for solving Dirichlet problems on the interior nodes.
    for (int i=0; i<nx; ++i)
      if (averagingFactor[i]==1)
        interiorDofs.push_back(i);
      else
        interfaceDofs.push_back(i);

    Aii = BMatrix(interiorDofs,A_);
    if (doTimings) timer.start("matrix factorization");
    AiiInv = DirectSolver<Vector,Vector>(Aii,DirectType::UMFPACK,MatrixProperties::POSITIVEDEFINITE);
    if (doTimings) timer.stop("matrix factorization");

    // Solve Dirichlet problem on interior nodes, completely eliminating the interior residual.
    solveDirichlet(u,f);    // solve and update solution
  }


  template <int m, class Scalar_, class CoarseScalar>
  std::vector<int> Subdomain<m,Scalar_,CoarseScalar>::neighbors() const
  {
    std::vector<int> neigh;               // extract the neighboring subdomain
    for (auto const& [s,_]: interfaces)   // numbers s from the list of interfaces
      neigh.push_back(s);
    return neigh;
  }

  // ----------------------------------------------------------------------------------------------

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::preRestrict()
  {
    // First compute the residual r = f-Au
    // We know that the inner residual ri vanishes. Thus we only need to compute the
    // interface part rD, i.e. rD = f_D - AiD ui - ADD uD
    for (int i: interfaceDofs)
    {
      r[i] = f[i];
      for (auto ci=A[i].begin(); ci!=A[i].end(); ++ci)
        r[i] -= *ci * u[ci.index()];
    }

    // In principle, we need to compute
    // ri_new = ri, rD_new = H^T ri, rDbar_new = rD-H^T ri, where H is the harmonic extension.
    // Now, since we maintain an approximate solution that satisfies ri = 0 (since in the
    // prolongation we solve the interior Dirichlet problem exactly (up to rounding errors)),
    // we know that ri_new = 0, rD_new = 0, and rDbar_new = rD. Thus, the computation is trivial
    // and can be omitted completely.


residual = r; // for debugging purposes
 }


  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getInterfaceResidual(int s, XVector& res) const
  {
    // First find the subdomain for which the data shall provided, and the
    // ids of the shared degrees of freedom.
    auto const& [_,ids] = *find_if(interfaces.begin(),interfaces.end(),
                                   [s](auto const& n) { return n.first==s; });
    res.resize(ids.size());
    for (int i=0; i<ids.size(); ++i)
      for (int j=0; j<components; ++j)
        res[i][j] = r[ids[i]*components+j];
  }

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::setInterfaceResidual(int s, XVector const& res)
  {
    // First find the subdomain for which the data is provided, and the
    // ids of the shared degrees of freedom.
    auto neighbor = find_if(interfaces.begin(),interfaces.end(),
                            [s](auto const& n) { return n.first==s; });
    assert(neighbor != interfaces.end());       // make sure we are given actual neighbor data
    auto const& ids = neighbor->second;         // get the affected dof indices
    assert(res.N()==ids.size());                // and make sure their number is as expected

    // Accumulate the provided values for later averaging by deviding by the number of contributions.
    for (int i=0; i<ids.size(); ++i)
      for (int j=0; j<components; ++j)
        avg[ids[i]*components+j] += res[i][j];
  }

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::restrict()
  {
    for (int i: interfaceDofs)                      // perform averaging of accumulated
      r[i] = (r[i]+avg[i]) / averagingFactor[i];    // interface coarse residuals

    avg = 0;
restrictedResidual = r;
  }


  // ----------------------------------------------------------------------------------------------

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getS(SMatrix& S) const
  {
    S.resize(nc,nc);
    X.resize(nx,nc);
    Vector b(nx+nc), v(nx+nc);
    for (int c=0; c<nc; ++c)          // For each constraint, create a
    {                                 // corresponding unit vector, solve
      b[nx+c] = 1;                    // the system with B, and
      Binv.apply(v,b);                // extract all the multiplier values,
      for (int row=0; row<nc; ++row)  // storing them in the corresponding
        S[row][c] = v[nx+row][0];     // column of S.
      for (int row=0; row<nx; ++row)  // TODO. restrict this to the interface dofs,
        X[row][c] = v[row][0];        //       i.e. xD.
      b[nx+c] = 0;
    }
  }


  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getSchurResidual(Vector& lambda) const
  {
    // For existing residual r (which, actually, vanishes on the interior dofs),
    // we solve [A C^T][x]   [r]
    //          [C    ][l] = [0], and return the l component. Denoting the KKT
    // matrix by B and projectors on the the x and l components by Px and Pl,
    // respectively, we can write this operation as l = Pl B^{-1} Px^T r.
    // This can be written as l = X^T r with the matrix X that has been computed
    // in the getS() method.
    lambda.resize(nc);
    lambda = 0;                     // lambda = X^T * r
    for (int i=0; i<nc; ++i)
      for (int j=0; j<nx; ++j)
        lambda[i] += X[j][i]*r[j];
  }

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::setCorrection(Vector const& dc)
  {
    // In principle, we need to solve
    // [A C^T][x]   [r ]
    // [C    ][l] = [dc], with r (where ri=0, so only rD counts) computed 
    // already in the restriction phase and dc provided now. Then we extract
    // just the x component of the solution. 
    // The system can either be solved as stated, or it can be solved separately for 
    // r and for dc.
    //
    // The advantage of simultaneous solve is clear: just one system solve.
    // This is what we do here.
    Vector rtmp(nx+nc), dutmp(nx+nc);
    for (int i=0; i<nx; ++i)
      rtmp[i] = r[i];
    for (int i=0; i<nc; ++i)
      rtmp[nx+i] = dc[i];
    Binv.apply(dutmp,rtmp);               // solve system for the coarse grid correction du
    for (int i=0; i<nx; ++i)
      du[i] = dutmp[i];
    
    
    // The advantage of the alternative, separte solves, is less obvious:
    // (i) r is available at the end of the restriction phase, and the system for 
    //     r could therefore be solved in parallel to the sequential process of 
    //     computing the coarse correction.
    // (ii) We may apply different solution strategies for the two systems. While the 
    //      small dimension of dc can be exploited for precomputing the solution operator
    //      and having a fast update, we may employ an approximative inverse for the r
    //      part.
    // As an example, we implement this strategy using a simple Jacobi update for the residual
    // part. This is, however, rather ineffective as a preconditioner, such that we provide this
    // just as an (inefficient) example. Smarter solvers, e.g., an H-matrix approximation of the 
    // solution operator, should perform much better.
    if (false)
    {
      // We need to solve, for given constraints value c, the system
      // we solve [A C^T][x]   [0]
      //          [C    ][l] = [c], and return the x component. Denoting the 
      // KKT matrix by B and projectors on the the x and l components by Px and Pl,
      // respectively, we can write this operation as x = Px B^{-1} Pl^T c = X c. 
      X.mv(dc,du);
      
      // Take a single Jacobi step for the residual part. Note that this does not  
      // satisfy the constraints (which is actually not necessary), but is anyways a 
      // rather poor solver.
      for (int i: interfaceDofs)
        du[i] += r[i] / A[i][i];
    }
    
    rawCorrection = du;
  }

  // ----------------------------------------------------------------------------------------------


  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getInterfaceValues(int s, XVector& values) const
  {
    // First find the subdomain for which the data shall provided, and the
    // ids of the shared degrees of freedom.
    auto const& [_,ids] = *find_if(interfaces.begin(),interfaces.end(),
                                   [s](auto const& n) { return n.first==s; });
    values.resize(ids.size());

    // Extract the interface values of the shared nodes.
    for (int i=0; i<ids.size(); ++i)
      for (int j=0; j<components; ++j)
        values[i][j] = du[ids[i]*components+j];
  }

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::setInterfaceValues(int s, XVector const& values)
  {
    // First find the subdomain for which the data is provided, and the
    // ids of the shared degrees of freedom.
    auto neighbor = find_if(interfaces.begin(),interfaces.end(),
                            [s](auto const& n) { return n.first==s; });
    assert(neighbor != interfaces.end());       // make sure we are given actual neighbor data
    auto const& ids = neighbor->second;         // get the affected dof indices
    assert(values.N()==ids.size());             // and make sure their number is as expected

    // Accumulate the provided values for later averaging by deviding by the number of contributions.
    for (int i=0; i<ids.size(); ++i)
      for (int j=0; j<components; ++j)
        avg[ids[i]*components+j] += values[i][j];
  }

  template <int m, class Scalar_, class CoarseScalar>
  Dune::FieldVector<double,2> Subdomain<m,Scalar_,CoarseScalar>::prolongate()
  {
    for (int i: interfaceDofs)                                // perform averaging of summed
      avg[i] = (du[i]+avg[i]) / averagingFactor[i] - du[i];   // interface residuals and subtract du,
                                                              // yields correction to be added to du.

    // Solve Dirichlet problem Av = 0 s.t. v = avg on boundary, i.e.
    // A_ii v_i = -A_Di avg_D, and add this to du, i.e. du_i = du_i - A_ii^{-1}A_Di avg_D
    // and du_D = du_D + avg_D. This yields a minimum-energy prolongation.
    Vector rhs(nx);                  // Form right hand side A_Di avg_D
    Vector avg_D(nx);
    for (int i: interfaceDofs)
      avg_D[i] = avg[i];
    A.smv(-1,avg_D,rhs);

    solveDirichlet(avg,rhs);            // Solve A_ii rhs -> avg_i (note that the _D part is not touched)
    du += avg;
    avg = 0;

    // Compute products dx*r and dx*A*dx, which are necessary for selecting
    // step lengths in gradient and conjugate gradient methods. Since r and the "restricted"
    // residual agree (as linear functionals) on the "fine" (i.e. continuous) subspace,
    // and du is (by averaging) contained in that subspace, we can directly use
    // the restricted residual.
    Dune::FieldVector<double,2> p;
    Vector Adu(nx);
    A.mv(du,Adu);
    for (int i=0; i<nx; ++i)
    {
      p[0] += r[i]*du[i];
      p[1] += du[i] * Adu[i];
    }

    return p;
  }

  // ----------------------------------------------------------------------------------------------


  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::updateSolution(Scalar alpha)
  {
    u.axpy(alpha,du);    // update solution
  }



  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getSolution(XVector& xout) const
  {
    for (int i=0; i<nx; ++i)                      // copy entries and reshape the
      xout[i/components][i%components] = u[i];    // vector blocking from scalar to
  }                                               // whatever is used externally

  template <int m, class Scalar_, class CoarseScalar>
  typename Subdomain<m,Scalar_,CoarseScalar>::XVector Subdomain<m,Scalar_,CoarseScalar>::getSolution() const
  {
    XVector xout(nx);
    getSolution(xout);
    return xout;
  }         

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getCorrection(XVector& xout) const
  {
    for (int i=0; i<nx; ++i)                       // copy entries and reshape the
      xout[i/components][i%components] = du[i];    // vector blocking from scalar to
  }                                                // whatever is used externally

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getRawCorrection(XVector& xout) const
  {
    for (int i=0; i<nx; ++i)                       // copy entries and reshape the
      xout[i/components][i%components] = rawCorrection[i];    // vector blocking from scalar to
  }                                                // whatever is used externally

  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::getResidual(XVector& xout) const
  {
    for (int i=0; i<nx; ++i)                             // copy entries and reshape the
      xout[i/components][i%components] = residual[i];    // vector blocking from scalar to
  }                                                      // whatever is used externally


  template <int m, class Scalar_, class CoarseScalar>
  void Subdomain<m,Scalar_,CoarseScalar>::solveDirichlet(Vector& x, Vector const& b) const
  {
    Vector bi(interiorDofs.size());
    for (int i=0; i<bi.N(); ++i)
      bi[i] = b[interiorDofs[i]];

    Vector xi(bi.N());
    AiiInv.apply(xi,bi);
    for (int i=0; i<bi.N(); ++i)
      x[interiorDofs[i]] = xi[i];
  }

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------
  
  template <int m, class Scalar>
  SharedMemoryDomain<m,Scalar>::SharedMemoryDomain(std::vector<Subdom>& subdomains_)
  : subdomains(&subdomains_)
  {
    for (int i=0; i<subdomains->size(); ++i)
      for (auto const& nid: (*subdomains)[i].neighbors())
        neighbors.push_back({i,nid});
  }
  
  template <int m, class Scalar>
  void SharedMemoryDomain<m,Scalar>::restrict()
  {
    parallelFor(0,subdomains->size(),[&](int i)
    {
      (*subdomains)[i].preRestrict();
    });

    // This is parallel, but essentially memory-bound. Anyways, the restriction is
    // not very time-consuming.
    parallelFor(0,subdomains->size(),[&](int to)
    {
      typename Subdom::Vector c;
      for (int from: (*subdomains)[to].neighbors())
      {
        (*subdomains)[from].getInterfaceResidual(to,c);
        (*subdomains)[to].setInterfaceResidual(from,c);
      }
    });

    // This is in parallel, but essentially memory-bound. 
    parallelFor(0,size(),[&](int i)
    {
      (*subdomains)[i].restrict();
    });
  }
  
  template <int m, class Scalar>
  Dune::FieldVector<double,2> SharedMemoryDomain<m,Scalar>::prolongate()
  {
    typename Subdom::Vector c;
    Dune::FieldVector<double,2> p;
    for (int from=0; from<subdomains->size(); ++from)
      for (int to: (*subdomains)[from].neighbors())
      {
        (*subdomains)[from].getInterfaceValues(to,c);
        (*subdomains)[to].setInterfaceValues(from,c);
      }

    std::vector<Dune::FieldVector<double,2>> ps(size());
    parallelFor(0,size(),[&](int i)
    {
      ps[i] = (*subdomains)[i].prolongate();
    });
    return std::accumulate(ps.begin(),ps.end(),Dune::FieldVector<double,2>());

    return p;
  }
  
  template <int m, class Scalar>  
  std::vector<DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>> SharedMemoryDomain<m,Scalar>::getS() const
  {
    std::vector<DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>> S(size());
    parallelFor(0,size(),[&](int i)
    {
      (*subdomains)[i].getS(S[i]);
    });
    return S;
  }

  template <int m, class Scalar>  
  auto SharedMemoryDomain<m,Scalar>::getSchurResiduals() const -> std::vector<Vector>
  {
    std::vector<Vector> rs(size());
    parallelFor(0,size(),[&](int i)
    {
      (*subdomains)[i].getSchurResidual(rs[i]);
    });

    return rs;
  }
  
  template <int m, class Scalar>  
  void SharedMemoryDomain<m,Scalar>::setCorrections(std::vector<Vector> const& dcs)
  {
    parallelFor(0,size(),[&](int i)
    {
      (*subdomains)[i].setCorrection(dcs[i]);
    });
  }
  
  template <int m, class Scalar>
  void SharedMemoryDomain<m,Scalar>::updateSolution(Scalar alpha)
  {
    parallelFor(0,size(),[&](int i)
    {
      (*subdomains)[i].updateSolution(alpha);
    });
  }
  
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <class Subdomain>
  BDDCSolver<Subdomain>::BDDCSolver(std::vector<Subdomain>& subdomains_,
                                    CoarseConstraints const& constraints_)
  : subdomains(&subdomains_)
  , domain(subdomains_)
  , constraints(constraints_)
  , coarseSolver(domain,constraints)
  {  }

  // ----------------------------------------------------------------------------------------------


  template <class Subdomain>
  double BDDCSolver<Subdomain>::solve()
  {
    if (domain.size()==1)          // With only one subdomain and thus no constraints, the
      return 0;                    // exact solution has already been computed on construction.

    // restriction: communicate residual between neighbors
    Timings::instance().start("restriction");
    domain.restrict();
    Timings::instance().stop("restriction");


    Timings::instance().start("coarse solve");
    coarseSolver.solve();
    Timings::instance().stop("coarse solve");


    // prolongation: communicate correction between neighbors
    Timings::instance().start("prolongation");
    auto p = domain.prolongate();
    Timings::instance().stop("prolongation");

    Timings::instance().start("solution update");
    double alpha = p[0] / p[1];
    std::cout << "choosing step length alpha = " << alpha << ", |dx|_A = " << std::sqrt(p[1]) << "\n";
    domain.updateSolution(alpha);
    Timings::instance().stop("solution update");

    // Return the energy norm of the gradient
    return std::sqrt(p[1]);
  }

}

#endif
