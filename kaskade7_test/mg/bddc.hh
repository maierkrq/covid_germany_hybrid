 
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

#ifndef BDDC_HH
#define BDDC_HH

#include "fem/firstless.hh"
#include "linalg/direct.hh"
#include "linalg/localMatrices.hh"
#include "utilities/enums.hh"

#include "dune/istl/bcrsmatrix.hh"

/**
 * \file
 * \brief  Basic components for (potentially distributed)
 *         Balanced Domain Decomposition with Constraints preconditioners.
 */

namespace Kaskade::BDDC
{
  /**
   * \brief A data structure describing the interfaces of a subdomain.
   *
   * It consists of a sequence of pairs \f$ (i,idx) \f$ of adjacent subdomain
   * number \f$ i \f$ and the (local to the current subdomain) degrees of
   * freedom (sorted in increasing order) that are shared with subdomain \f$ i \f$.
   * The sequence shall be sorted by increasing subdomain number.
   */
  using Interfaces = std::vector<std::pair<int,std::vector<int>>>;

  /**
   * \brief A data structure describing coarse constraints.
   *
   * Each constraint is described by a pair \f$ (s,m) \f$, where \f$ s \f$
   * is a list of subdomains affected by the constraint, and \f$ m \f$ is
   * the dimension of the constraint.
   */
  using CoarseConstraints = std::vector<std::pair<std::vector<int>,int>>;

  /**
   * \brief A pair \f$ (s,i) \f$ of subdomain id \f$ s \f$ and local dof index
   *        \f$ i \f$ denotes a subdomain-local reference to a global dof.
   *
   * \note The data type int (usually 32 bit wide) is thought to be big enough
   *       for both subdomain number s and local index i, since both the
   *       number of subdomains and their size is expected to be at most moderate.
   *       Moreover, efficiency considerations for two-level BDDC preconditioners
   *       suggest that s and i should be of comparable size, such that together
   *       they cover huge problems with up to 10^15 dofs.
   */
  struct LocalDof
  {
    int s;
    int i;
  };

  /**
   * \ingroup multigrid
   * \brief Discrimination of different interface types.
   */
  enum InterfaceType { CORNER = 1, EDGE = 2, FACE = 4 };

  /**
   * \ingroup multigrid
   * \brief A class for computing and providing BDDC coarse space constraints (interface averages).
   * 
   * \tparam m number of components of the variable in the system (usually 1, 2, or 3).
   * \tparam Index integral number type used for indices
   */
  template <int m=1, class Index=size_t>
  class InterfaceAverages
  {
  public:
    constexpr static int components = m;

    /**
     * \brief Constructor
     * \param sharedDofs a set of (global) dofs shared by at least two subdomains. For each such shared dof,
     *                   the vector entry contains a list of LocalDof pairs, stating which local index in
     *                   which subdomain corresponds to the global dof.
     * \param subdomainSize a list of the total number of subdomain-local dofs
     * \param wireframe if true, include averages over all sets of dofs shared by
     *                  the same group of subdomains, usually corners and edges in 3D
     */
    InterfaceAverages(std::vector<std::vector<LocalDof>> const& sharedDofs,
                      std::vector<Index> const& subdomainSize,
                      int types = CORNER | FACE);

    /**
     * \brief For the given subdomain, tells with which neighboring domain which dofs are shared.
     *
     * \param id the subdomain number
     */
    Interfaces const& interfaces(int id) const;

    /**
     * \brief Creates the subdomain-local constraints defining the coarse space.
     *
     * This creates a constraint representing equal averages of shared nodes between
     * adjacent subdomains.
     *
     * \param id the number of the current subdomain
     */
    template <class Scalar=double>
    NumaBCRSMatrix<Dune::FieldMatrix<Scalar,m,m>> coarseConstraint(int id) const;

    /**
     * \brief Provides information about subdomains sharing coarse constraints.
     *
     * \return a sequence of triples \f$ (s_1, s_2,n) \f$ of subdomains \f$ s_1 < s_2 \f$
     *         sharing \f$ n \f$ coarse constraints.
     */
    CoarseConstraints coarseConstraints() const;

  private:
    std::vector<Interfaces> ifs;
    CoarseConstraints ccn;
    std::vector<Index> subdomainSize;

    std::vector<std::vector<std::vector<int>>> cifs; // coarse interfaces
  };




  /**
   * \ingroup multigrid
   * \brief A class representing small KKT systems on a single subdomain.
   * 
   * This is intended to work together with BDDCSolver and KKTSolver, and represents 
   * subdomains, providing subdomain-local storage and processing capabilities. This is 
   * the basis for distributed domain decomposition solvers.
   *
   * A Subdomain object represents the system
   * \f[
   *   B = \begin{bmatrix} A & C^T \\ C \end{bmatrix},
   * \f]
   * and provides in particular access to the "Schur complement"
   * \f[
   * S =   \begin{bmatrix} 0 & I \end{bmatrix}
   *       \begin{bmatrix} A & C^T \\ C \end{bmatrix}^{-1}
   *       \begin{bmatrix} 0 \\ I \end{bmatrix}.
   * \f]
   * Since \f$ A \f$ is only positive semidefinite in general, the true Schur complement
   * \f$ - C A^{-1} C^T \f$ cannot be formed, but is equal in case of invertible \f$ A \f$.
   *
   * The subdomain holds \f$ A_i \f$, \f$ f_i \f$, and the constraints
   * block \f$ C \f$. It also holds a local solution vector \f$ u_i \f$.
   * For the overall picture, we refer to KKTSolver.
   *
   * \tparam m the number of vectorial components in \f$ x \f$ (usually 1
   *           for scalar systems, and 2 or 3 for vectorial ones)
   * \tparam Scalar the field type \f$ K \f$, usually \c double
   * \tparam CoarseScalar the field type \f$ K_c \f$ to be used for representing the coarse solver
   *                      factorization, usually \c double, but could also be \c float
   * 
   * \todo UMFPACK does not provide \c float factorization. Extract L and U factors and convert to float?
   *
   * This is a reference implementation.
   */
  template <int m, class Scalar_=double, class CoarseScalar=Scalar_>
  class Subdomain
  {
  public:
    /**
     * \brief Number of (vectorial) components of \f$ u \f$.
     */
    static const int components = m;

    /**
     * \brief The scalar field type used in coefficients.
     */
    using Scalar = Scalar_;

    using XEntry = Dune::FieldVector<Scalar,components>;
    using XVector = Dune::BlockVector<XEntry>;
    using Vector = Dune::BlockVector<Dune::FieldVector<Scalar,1>>;


    using AMatrix = NumaBCRSMatrix<Dune::FieldMatrix<Scalar,components,components>>;

    /// A matrix type with scalar entries, used for the Schur complement matrix.
    using SMatrix = DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>;

    /**
     * \brief Constructor.
     *
     * This initializes the solution approximation to zero, which implies an admissible
     * state over all subdomains.
     *
     * \param id the subdomain number
     * \param A the subdomain-local stiffness matrix
     * \param f the subdomain-local right hand side
     * \param interfaces a list of pairs (s,ids) of neighboring subdomain number and
     *                   the degrees of freedom shared with this subdomain (in the given order)
     * \param C
     */
    template <int cComponents>
    Subdomain(int id, AMatrix const& A, XVector const& f,
              Interfaces const& interfaces,
              NumaBCRSMatrix<Dune::FieldMatrix<Scalar,cComponents,components>> const& C);

    template <class Interface>
    Subdomain(int id, AMatrix const& A, XVector const& f, Interface const& ifs);

    /**
     * \brief Reports a list of its neighbor subdomains.
     */
    std::vector<int> neighbors() const;
    
    /**
     * \name Restriction support.
     * \brief These methods implement the "restriction" of the residual to the "coarse" space.
     *
     * The restriction is defined as the transpose of the prolongation, i.e.
     * \f$ \Pi^T = \pi^T E^T \f$ with averaging \f$ \pi \f$ and extension
     * \f$ E \f$.
     *
     * The extension operator \f$ E^T \f$ acts locally, realizing the dual harmonic extension
     * \f[
     * \begin{bmatrix} r_i \\ r_D \\ \bar r_D \end{bmatrix}_{\rm new}
     * = \begin{bmatrix} I &  \\
     *                   H^T  & \\
     *                   -H^T & I \end{bmatrix}
     *   \begin{bmatrix} r_i \\ r_D  \end{bmatrix}
     * \f]
     * working on interior residual \f$ r_i \f$,  interface residual \f$ r_D \f$, and averaged
     * interface residual \f$ \bar r_D \f$. Here, \f$ H = A_{ii}^{-1} A_{iD} \f$ denotes the
     * harmonic extension of interface values into the interior. This local operation is
     * performed by the preRestrict() method.
     *
     * The averaging (for a single degree of freedom shared by \f$ k \f$ subdomains)
     * is \f[ \pi^T = \frac{1}{k} \begin{bmatrix} 1 & \dots & 1 \\ \vdots & \ddots & \vdots \\
     *                                          1 & \dots & 1 \end{bmatrix}, \f]
     * therefore symmetric, and acts on the interface residual values \f$ r_D \f$ only. It involves communication
     * between subdomains, exchanging their interface residual \f$ r_D \f$, which is realized by
     * the getInterfaceResidual() and setInterfaceResidual() methods.
     * Finally, restrict() must be called as a sentinel to denote the end of communication, and
     * to perform the actual averaging.
     *
     * @{
     */
    
    /**
     * \brief First step of the restriction: Dirichlet solve
     *
     * This computes the transpose of the prolongation, i.e., with a suitable sorting of
     * interior and bounday degrees of freedom, the
     */
    void preRestrict();

    /**
     * \brief Gets interface residual that shall be shared with neighboring subdomain.
     * \param s the number of the neighboring subdomain
     *
     */
    void getInterfaceResidual(int s, XVector& res) const;

    /**
     * \brief Sets interface residual of the correction from neighboring subdomains.
     * \param s the number of the neighboring subdomain
     *
     */
    void setInterfaceResidual(int s, XVector const& res);

    /**
     * \brief Performs the restriction of the residual vector.
     * \param s the number of the neighboring subdomain
     */
    void restrict();
    
    /// @}

    /**
     * \name Coarse solve support.
     * \brief These methods implement the subdomain-local parts of the coarse space solver.
     * @{
     */

    /**
     * \brief Compute \f$ S \f$.
     *
     * This computes
     * \f[
     * S = \begin{bmatrix} 0 & I \end{bmatrix}
     *     \begin{bmatrix} A & C^T \\ C \end{bmatrix}^{-1}
     *     \begin{bmatrix} 0 \\ I \end{bmatrix}
     * \f]
     */
    void getS(SMatrix& S) const;

    /**
     * \brief Provides the Schur residual.
     *
     * This solves the system
     * \f[ \begin{bmatrix} A & C^T \\ C \end{bmatrix}
     *     \begin{bmatrix} x+\delta x \\ \lambda +\delta\lambda \end{bmatrix} =
     *     \begin{bmatrix} f\\ 0 \end{bmatrix}
     * \f]
     * and returns the solution component \f$ \delta\lambda \f$ in the provided vector \c lambda.
     * The vectors \f$ x \f$ and \f$ \lambda \f$ are maintained internally as solution approximation.
     *
     * \param[out] lambda
     */
    void getSchurResidual(Vector& lambda) const;

    /**
     * \brief Computes the correction vector \f$ [\delta x, \delta\lambda] \f$.
     */
    void setCorrection(Vector const& c);
    
    /// @}

    /**
     * \name Prolongation support.
     * \brief These methods implement the "prolongation" of "coarse" space corrections to the 
     *        "fine" space.
     *
     * The prolongation is defined as \f$ \Pi = E\pi \f$ with averaging \f$ \pi \f$ and extension
     * \f$ E \f$. The averaging (for a single degree of freedom shared by \f$ k \f$ subdomains)
     * is \f[ \pi = \frac{1}{k} \begin{bmatrix} 1 & \dots & 1 \\ \vdots & \ddots & \vdots \\
     *                                          1 & \dots & 1 \end{bmatrix} \f]
     * and acts on the interface degrees of freedom \f$ x_D \f$ only. It involves communication
     * between subdomains, exchanging their interface values \f$ x_D \f$, which is realized by
     * the getInterfaceValues() and setInterfaceValues() methods.
     *
     * The extension operator \f$ E \f$ acts locally, realizing the harmonic extension
     * \f[
     * \begin{bmatrix} x_i \\ x_D \end{bmatrix}_{\rm new}
     * = \begin{bmatrix} I & H & - H \\
     *                     &   &   I \end{bmatrix}
     *   \begin{bmatrix} x_i \\ x_D \\ \bar x_D \end{bmatrix}
     * \f]
     * working on interior dofs \f$ x_i \f$,  interface dofs \f$ x_D \f$, and averaged
     * interface values \f$ \bar x_D \f$. Here, \f$ H = A_{ii}^{-1} A_{iD} \f$ denotes the
     * harmonic extension of interface values into the interior. This local operation is
     * performed by the prolongate() method.
     *
     * @{
     */
    
    /**
     * \brief Gets boundary values that shall be shared with neighboring subdomain.
     * \param s the number of the neighboring subdomain
     *
     */
    void getInterfaceValues(int s, XVector& values) const;

    /**
     * \brief Sets boundary values of the correction from neighboring subdomains.
     * \param s the number of the neighboring subdomain
     *
     */
    void setInterfaceValues(int s, XVector const& values);

    /**
     * \brief Perform averaging on the shared interfaces, projecting the solution to
     *        the subspace of functions continuous across the interfaces.
     *
     * \return a pair \f$ (\delta x ^T r, \delta x^T A \delta x) \f$ for the correction
     *         \f$ \delta x \f$ and the residual \f$ r \f$.
     */
    Dune::FieldVector<double,2> prolongate();
    
    /// @}

    /**
     * \brief Updates the solution approximation.
     */
    void updateSolution(Scalar alpha);

    /**
     * \brief Provides the current solution vector \f$ x \f$.
     *
     * \param[out] x on exit, contains the current solution approximation
     */
    void getSolution(XVector& x) const;
    
    /**
     * \brief Provides the current solution vector \f$ x \f$.
     *
     * \return the current solution approximation
     * 
     * Note that the result is returned by value, which may be less efficient than providing 
     * a result vector to be filled.
     */
    XVector getSolution() const;

    /**
     * \brief Provides the current correction vector \f$ dx \f$.
     *
     * \param[out] dx on exit, contains the current correction
     */
    void getCorrection(XVector& dx) const;
    void getRawCorrection(XVector& dx) const;
    void getResidual(XVector& rr) const;
    Vector const& getRestrictedResidual() const { return restrictedResidual; }

  private:
    using BMatrix = NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>>;

    // These variables describe problem size and connectivity
    int id;                               // my own subdomain number
    int nx, nc;                           // number of dofs and of constraints
    std::vector<std::pair<int,                            // list of pairs (s,ids) of adjacent 
                          std::vector<int>>> interfaces;  // subdomain s and shared dofs
    
    Vector f;                            // the right hand side
    Vector u;                            // approximate solution [x,lambda]
    Vector r;                            // current residual [f-Au,0], restricted after restriction

    Vector du;                           // correction [dx,0]
    
    BMatrix A;
    BMatrix C;
    DirectSolver<Vector,Vector> Binv;    // [ A C^T; C 0 ]^{-1}
    
  
    Vector avg;                          // temporary vector for accumulating interface values
    std::vector<int> averagingFactor;    // for each dof the number of incident subdomains
    std::vector<int> interiorDofs;       // the interior dof indices
    std::vector<int> interfaceDofs;      // the boundary/interface dof indices
    
    // for debugging purposes
    Vector rawCorrection, residual, restrictedResidual;

BMatrix Aii;
DirectSolver<Vector,Vector> AiiInv;
std::vector<int> cIndices;
mutable DynamicMatrix<double> X; // const hack, move computation to constructor

    // Solves the homogeneous Dirichlet problem A_ii x_i = r_i for the interior nodes only.
    // Non-zero boundary values x_D can be imposed by exploiting the splitting
    // [Aii AiD] [x_i x_D]^T = r_i <=> Aii x_i = r_i - AiD xD, i.e. modifying the right hand side.
    // The interface nodes of x are not touched, nor is the constraint/multiplier part.
    void solveDirichlet(Vector& x, Vector const& rhs) const;
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------
  
  template <int m, class Scalar_=double>
  class SharedMemoryDomain
  {
  public:
    using Scalar = Scalar_;
    using Subdom = Subdomain<m,Scalar>;
    using Vector = typename Subdom::Vector;
    
    SharedMemoryDomain(std::vector<Subdom>& subdomains);
    
    void restrict();
    Dune::FieldVector<double,2> prolongate();
    int size() const
    {
      return subdomains->size();
    }
    
    std::vector<DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>> getS() const;
    std::vector<Vector> getSchurResiduals() const;
    void setCorrections(std::vector<Vector> const& dcs);

    void updateSolution(Scalar alpha);
    
  private:
    std::vector<Subdom>* subdomains;
    std::vector<std::pair<int,int>> neighbors;    
  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup multigrid
   * \brief A class for solving a block-diagonal-structured KKT system
   *
   * We treat convex quadratic-linear optimization  of the form
   * \f[
   *  \min_{x_i \in \mathbb{K}^{n_i}} \frac{1}{2}x_i^T A_i x_i - f_i^T x_i \quad
   *  \text{s.t.}\quad C_{ij} x_i = C_{ji} x_j \quad \text{for }(i,j)\in\Gamma,
   * \f]
   * where \f$ A_i \in \mathbb{K}^{n_i\times n_i}\f$ is positive semidefinite,
   * \f$ C_{ij} \in \mathbb{K}^{m_{ij} \times n_i}\f$, \f$ C_{ji}\in\mathbb{K}^{m_{ij} \times n_j}\f$,
   * and \f$ \Gamma \subset \{ (i,j) \in \mathbb{N}^2 \mid i < j \} \f$.
   * For the algorithmic approach and data structures, see below.
   *
   * \tparam components the number of vectorial components in \f$ x \f$.
   *                    1 for scalar problems, 2 or 3 for vectorial problems.
   * \tparam Scalar     the field type \f$ K \f$, usually \c double
   *
   * ## Algorithmic approach: Balanced Domain Decomposition by Constraints
   *
   * We would like to do block elimination of the block-diagonal matrix \f$ A \f$ in the KKT system
   * \f[ \begin{bmatrix} A & C^T \\ C & \end{bmatrix} \begin{bmatrix}u \\ \lambda\end{bmatrix}
   *     = \begin{bmatrix} f \\ 0 \end{bmatrix}
   * \f]
   * with block-diagonal structure
   * \f[ A = \begin{bmatrix} A_1 & & & \\ & A_2 & & \\ & & \ddots & \\ &&& A_n \end{bmatrix} , \quad
   *     f = \begin{bmatrix} f_1 \\ f_2 \\ \vdots \\ f_n \end{bmatrix},
   * \f]
   * but since \f$ A_i \f$ is in general only positive semidefinite (though \f$ A \f$
   * is definite on the nullspace of \f$ C \f$), this fails.
   *
   * ### Reformulation
   *
   * We therefore introduce slacks
   * \f$ s_{ij} \in \mathbb{K}^{m_{ij}} \f$ for \f$ (i,j)\in\Gamma \f$ as additional primal variables,
   * and reformulate equivalently
   * \f[
   *  \min_{x_i \in \mathbb{K}^{n_i}, s_{ij}\in\mathbb{K}^{m_{ij}}} \frac{1}{2} x_i^T A_i x_i - f_i^T x_i \quad
   *  \text{s.t.}\quad \begin{aligned} C_{ij} x_i &= s_{ij} \\
   *                                   C_{ji} x_j &= s_{ij}\end{aligned} \quad \text{for }(i,j)\in\Gamma,
   * \f]
   * The corresponding KKT system (suggestively reordered such that the slacks appear last) is then
   * \f[
   *  \begin{bmatrix} A_1    & C_{1*}^T & \\
   *                  C_{1*} &          &        &        &          & -I_1 \\
   *                         &          & \ddots \\
   *                         &          &        & A_k    & C_{k*}^T \\
   *                         &          &        & C_{k*} &          & -I_k \\
   *                         & -I_1^T   &        &        & -I_k^T
   *  \end{bmatrix} \begin{bmatrix} x_1 \\ \lambda_1 \\ \vdots \\ x_k \\ \lambda_k \\ s \end{bmatrix}
   *  = \begin{bmatrix} f_1 \\ 0 \\ \vdots \\ f_k \\ 0 \\ 0 \end{bmatrix}.
   * \f]
   * Here, the constraints are grouped, such that
   * \f[ C_{i*} = \begin{bmatrix} C_{ij_1} \\ \vdots \\ C_{ij_{m_i}} \end{bmatrix} \quad \text{for }
   *   (ij_k) \in \Gamma \text{ or } (j_k i)\in\Gamma.
   * \f]
   * Similarly, \f$ I_i \f$ does not denote the identity, but a block matrix composed of
   * zero and identity blocks, with exactly one identity block in every row and at most
   * one identity block in every column. Note that now the diagonal blocks
   * \f[ B_i := \begin{bmatrix} A_i    & C_{i*}^T \\ C_{i*} \end{bmatrix} \f]
   * are invertible small KKT systems, such that a block elimination can be done, leading to a Schur
   * complement system for the primal variables \f$ s_{ij} \f$.
   *
   * ### Distributed solve
   *
   * The KKTSolver itself does not know anything about the matrices \f$ A_i, C_{i*}\f$, or the right
   * hand sides \f$ f_i \f$. It only knows about \f$ \Gamma \f$, i.e. which subdomains share a
   * constraint, and about the dimension \f$ m_{ij} \f$ of the constraints. This is enough to
   * (i) form a slack vector \f$ s \f$ and (ii) to decompose it into appropriate right hand sides for
   * the subdomain solves, i.e. defining \f$ I_i \f$.
   *
   * Using the local "Schur complements"
   * \f[
   *    S_i = \begin{bmatrix} 0 & I \end{bmatrix} B_i^{-1} \begin{bmatrix} 0 \\ I \end{bmatrix}
   * \f]
   * (which for invertible \f$ A_i \f$ equal \f$ - C_{i*}A_i^{-1}C_{i*}^T \f$ and are negative definite),
   * the Schur complement system for \f$ s \f$ is
   * \f[
   *  Ss := \left(-\sum_i I_i^T S_i^{-1} I_i \right) s
   *      = \sum_i \begin{bmatrix} 0 & I_i^T \end{bmatrix} B_i^{-1} \begin{bmatrix} f_i \\ 0 \end{bmatrix}.
   * \f]
   * Note that the Schur complement \f$ S \f$ is again *positive definite*.
   * We provide \f$ c_{i*} = I_i s \f$ to the subdomain systems (usually this is some kind of Dirichlet
   * boundary data), and receive the \f$ \lambda_i \f$ solution component from them (which usually is
   * some kind of Neumann boundary data), which is then scattered into a global slack vector using
   * \f$ I_i^T \f$.
   * This defines the action of the Schur complement \f$ B \f$ (a kind of Dirichlet-to-Neumann operator)
   * and allows for a Krylov or gradient algorithm to solve for \f$ s \f$.
   *
   * \related Subdomain
   */
  template <class Domain>
  class KKTSolver
  {
  public:
    using Scalar = typename Domain::Scalar;
    using Constraint = std::tuple<int,int,int>;

  private:
    using Vector = Dune::BlockVector<Dune::FieldVector<Scalar,1>>;

  public:
    /**
     * \brief Constructor
     *
     * \param domain a collection of subdomain objects representing the systems \f$ B [b \lambda]^T = \tilde f\f$
     * \param constr  a sequence of \f$ (i,j,m_{ij}) \f$ with \f$ (i,j)\in\Gamma\f$ and \f$ m_{ij} \f$ the
     *                    number of constraints shared by subdomains \f$ i < j \f$
     */
    KKTSolver(Domain& domain_,
              CoarseConstraints const& constr)
    : domain(&domain_)
    , constraints(domain->size())
    , localConstraintsSize(domain->size())
    , sOffsets(constr.size()+1)
    {
      // Compute the offsets of individual slack vectors s_ij into the global slack vector
      // and the constraints in which each subdomain is involved.
      sOffsets[0] = 0;
      for (int k=0; k<constr.size(); ++k)
      {
        auto const& [s,m] = constr[k];
        for (auto i: s)
        {
          assert(0<=i && i<domain->size());
          constraints[i].push_back(k);
          localConstraintsSize[i] += m;
        }
        sOffsets[k+1] = m;
      }
      std::partial_sum(begin(sOffsets),end(sOffsets),begin(sOffsets));

      // Create Schur complement matrix B. This is sparse, since each constraint, i.e.
      // slack variable, is only affected by two subdomains.
      int const ns = sOffsets.back();                                 // Total number of slacks.
      if (ns == 0)                                                    // If there are no coarse constraints,
        return;                                                       // we're done - there's nothing to do.
      NumaCRSPatternCreator<> sCreator(ns,ns,false,2);
      for (auto const& ss: constraints)                               // For each subdoman find all pairs
        for (int sRow: ss)                                            // of incident slacks and add the
          for (int sCol: ss)                                          // corresponding dense block into
            sCreator.addDenseBlock(sOffsets[sRow],sOffsets[sRow+1],   // the sparsity pattern of S.
                                   sOffsets[sCol],sOffsets[sCol+1]);

      using Entry = Dune::FieldMatrix<Scalar,1,1>;
      NumaBCRSMatrix<Entry> S(sCreator);                              // And then create the matrix.

      // Assemble the Schur complement matrix.
      using SortedIdx = std::vector<std::pair<int,int>>;
      LocalMatrices<Entry,false,SortedIdx,SortedIdx> lm(S);
      auto localSchurComplements = domain->getS();
      for (int i=0; i<domain->size(); ++i)
      {
        SortedIdx rows;
        int localIndex = 0 ;
        for (int k: constraints[i])
          for (int l=sOffsets[k]; l<sOffsets[k+1]; ++l)
            rows.push_back({l,localIndex++});
        assert(std::is_sorted(rows.begin(),rows.end(),FirstLess()));
        lm.push_back(rows,rows);
        lm.back() = localSchurComplements[i];                         // add local stiffness matrix Si
      }
      lm.scatter();
      S *= -1.0;                                                      // take care of correct Schur complement sign

      Sinv = DirectSolver<Vector,Vector>(S,DirectType::UMFPACK,
                                         MatrixProperties::SYMMETRICSTRUCTURE);
    }

    /**
     * \brief solves the system.
     *
     * After the method returns, the subdomains contain their individual solution vectors
     * which can be retieved using Subdomain::getSolution().
     */
    void solve()
    {
      int const ns = sOffsets.back();                 // Total number of slacks.

      if (ns==0)
        return;

      // Obtain right hand side of the Schur complement system. This is the sum of right
      // hand sides from the subdomains, scattered into appropriate positions in the
      // right hand side vector f.
      Timings& timer = Timings::instance();
      timer.start("get Schur rhs");
      Vector f(ns);
      getRightHandSide(f);
      timer.stop("get Schur rhs");

      timer.start("solve Schur system");
      Vector s(ns);                                   // Create zero initial slack vector
      Sinv.apply(s,f);                                // and solve the global Schur system.
      timer.stop("solve Schur system");

      timer.start("set corrections");
      std::vector<Vector> dcs;
      for (int i=0; i<domain->size(); ++i)
      {
        Vector si(localConstraintsSize[i]);           // gather the slacks corresponding to the
        int pos = 0;                                  // constraints affecting subdomain i
        for (auto const k: constraints[i])            // by concatenating all the slacks
        {                                             // corresponding to constraints in which
          int const off = sOffsets[k];                // subdomain i is involved.
          int const nsk = sOffsets[k+1]-off;
          for (int l=0; l<nsk; ++l)
            si[pos++] = s[l+off];
        }
        dcs.push_back(si);
      }
      domain->setCorrections(dcs);
      timer.stop("set corrections");
    }


  private:

    void getRightHandSide(Vector& f) const
    {
      f = 0;
      auto floc = domain->getSchurResiduals();
      
      for (int i=0; i<domain->size(); ++i)
      {
        int pos = 0;                                  // scatter the lambda solution component
        for (auto const k: constraints[i])            // into the right hand side vector
        {
          int const off = sOffsets[k];
          int const nsk = sOffsets[k+1]-off;
          for (int l=0; l<nsk; ++l)
            f[l+off] += floc[i][pos++];
        }
      }
    }


    Domain* domain;
    std::vector<std::vector<int>> constraints;        // involved constraints by subdomain id
    std::vector<int> localConstraintsSize;            // total number of constraints for subdomain id
    std::vector<int> sOffsets;                        // offset of slack vector in total slacks
    DirectSolver<Vector,Vector> Sinv;                 // inverse of Schur complement
  };


  /**
   * \ingroup multigrid
   * \brief Balanced Domain Decomposition with Constraints preconditioner
   *
   * ## Auxiliary space approach to BDDC
   * We solve symmetric positive definite linear equation systems \f$ A x = b \f$ in \f$ \mathbb{R}^n\f$
   * by Balanced Domain Decomposition with Constraints (BDDC). There, the system is
   * rewritten equivalently by doubling certain degrees of freedom and requiring their equality
   * as a constrained block-diagonal system
   * \f[
   *  \begin{bmatrix} A_0 & C^T \\
   *                  C_0 & 0
   *  \end{bmatrix} \begin{bmatrix} x_0 \\ \lambda \end{bmatrix}
   *  = \begin{bmatrix} f_0  \\ 0 \end{bmatrix}, \quad
   *  A_0 = \begin{bmatrix} A_1 & \\ & A_2 \\ & & \ddots \end{bmatrix} \in\mathbb{R}^N.
   * \f]
   * The constraint  \f$ C \f$ contains in each row exactly one entry 1, and one entry -1,
   * making clear that some degrees of freedom are simply split in two (or even more, if
   * several blocks overlap). With that,
   * we can identify \f$ \mathbb{R}^n \f$ with \f$\mathrm{ker} C\f$ and assume that
   * \f$ x^T A_0 x = x^T A x \f$ for \f$ x\in \mathrm{ker} C\f$.
   *
   * The driving use case is solving elliptic PDEs discretized by finite elements by
   * domain decomposition without overlap. The ansatz space is the space of piecewise-continuous
   * functions, whereas the nullspace of the constraints is the usual space of continuous
   * finite element functions. In such a setting, \f$ A_0 \f$ is relatively easy to invert,
   * due to its block-diagonal structure, whereas solving globally coupled systems with \f$ A \f$
   * is expensive.
   *
   * In addition to the constraint \f$ C \f$,
   * we define a smaller set of constraints given as \f$ \tilde C \f$, such that
   * \f$ V:=\mathrm{ker} C \subset \mathrm{ker}\tilde C = \tilde V \f$. We define a projector
   * \f$ \Pi : V_0 \to V \f$ (known as averaging in BDDC literature and "prolongation" in subspace
   * correction literature) and its adjoint \f$ \Pi^T : V^* \to V_0^* \f$ known as "restriction".
   *
   * We then define the "additive" auxiliary space preconditioner \f[ P = \Pi A_0^{-1} \Pi^T, \f]
   * see J. Xu: The auxiliary space method and optimal multigrid preconditioning techniques for
   * unstructured grids. Computing 56:215-235, 1996.
   *
   * ### Prolongation definition
   * For defining the prolongation \f$ \Pi \f$, we take a closer look at the structure of
   * \f$ A_0 \f$ and \f$ C \f$. We assume \f$ A_0 \f$ to be block-diagonal, consisting of
   * square blocks \f$ A_1,\dots,A_m \f$. Accordingly, the degrees of freedom \f$ \{1,\dots,N\}\f$
   * are partitioned into pairwise disjoint sets \f$ \Xi_i \f$. We say degrees of freedom
   * \f$ k\in \Xi_i, l\in\Xi_j\f$ are shared between blocks \f$ i, j \f$ if
   * \f$ x \in V\f$ enforces \f$ x_k = x_l \f$, and call \f$ l \f$ and \f$ k \f$ coupled.
   * Coupled is an equivalence relation, leading to equivalence classes \f$ [k]\f$ of coupled degrees
   * of freedom.
   *
   * Given \f$ x\in \tilde V \f$, a simple prolongation is defined by
   * \f[ (\Pi x)_k = \frac{1}{|[k]|} \sum_{l\in[k]} x_l. \f]
   * In matrix notation, and assuming an ordering of degrees of freedom according to coupling,
   * the projector is block diagonal with symmetric blocks
   * \f[ \pi = \frac{1}{|[k]|} \begin{bmatrix} 1 & \dots & 1 \\ \vdots & \ddots & \vdots \\
   *                            1 & \dots & 1 \end{bmatrix},
   * \f]
   * such that \f$ \Pi^T = \Pi \f$. Unfortunately, since the condition number of the preconditioned
   * system is essentially \f$ \| \Pi \|_A^2 \f$, this simple approach, which just modifies the
   * interface nodes, is mesh-dependent for elliptic problems (notably it depends on \f$ h/H \f$
   * in domain decomposition notation).
   *
   * An optimal subdomain-local prolongation postprocessing is not to add just the difference
   * to the averaged interface nodes, but to extend these harmonically to the interior of the
   * subdomain. This yields an update with minimal energy norm. (Subject to the simple averaging
   * of interface nodes, of course. If they are modified as well, better prolongations such as
   * deluxe averaging can be created - but this is no longer subdomain-local.)
   *
   *
   *
   *
   *
   *
   * If we split the degrees of freedom in a subdomain into interior nodes \f$ x_i \f$ and
   * interface or Dirichlet nodes \f$ x_D \f$, and denote by \f$ \bar x_D \f$ the averaged
   * interface values, then the harmonic extension can be written in assignment form as
   * \f[
   *  \begin{aligned}
   *    q &\gets \bar x_D - x_D \\
   *    x_i & \gets x_i - H q \\
   *    x_D & \gets x_D + q
   *  \end{aligned}
   * \f]
   * with the harmonic extension operator \f$ H = A_{ii}^{-1} A_{iD} \f$,
   * or in matrix form as
   * \f[
   * \begin{bmatrix} x_i \\ x_D \end{bmatrix}_{\rm new}
   * = \begin{bmatrix} I & H & - H \\
   *                     &   &   I \end{bmatrix}
   *   \begin{bmatrix} x_i \\ x_D \\ \bar x_D \end{bmatrix}.
   * \f]
   * For the restriction, its transpose
   * \f[
   *  \begin{bmatrix} r_i \\ r_D \\ \bar r_D \end{bmatrix}_{\rm new}
   *  = \begin{bmatrix} I \\
   *                    H^T \\
   *                    -H^T & I \end{bmatrix}
   *    \begin{bmatrix} r_i \\ r_D \end{bmatrix}
   * \f]
   * is needed, and defines the interface residual \f$ \bar r_D \f$ to be averaged between
   * neighboring subdomains. In assignment form, it reads
   * \f[
   *  \begin{aligned}
   *    q &\gets H^T r_i \\
   *    \bar r_D &\gets r_D -q \\
   *    r_D &\gets q
   *  \end{aligned}
   * \f]
   * with \f$ r_i \f$ untouched.
   *
   * For a problem with just two subdomains partitioned as \f$ [x_i,x_D,y_i,y_D]\f$,
   * the prolongation matrix reads
   * \f[
   *  \Pi = \begin{bmatrix} I &  H & -H  \\
   *                          &    &    & I \\
   *                          &    &    &    & I & H & -H \\
   *                          &    &    &    &   &   & I \end{bmatrix}
   *        \begin{bmatrix} I \\
   *                          & I \\
   *                          & I/2 &    &I/2 \\
   *                          &     &  I \\
   *                          &     &    & I \\
   *                          & I/2 &    & I/2  \end{bmatrix}
   *      = \begin{bmatrix} I & H/2  &   & -H/2 \\
   *                          & I/2  &   &  I/2 \\
   *                          & -H/2 & I & H/2 \\
   *                          &  I/2 &   & I/2 \end{bmatrix}.
   * \f]
   *
   * ### Coarse space correction
   *
   * On the (actually larger) "coarse" space, the less constrained system
   * \f[
   * \begin{bmatrix} A_0 & \tilde C^T \\ \tilde C \end{bmatrix}
   * \begin{bmatrix} \delta x \\ \delta \lambda \end{bmatrix}
   * = \begin{bmatrix} r \\ c \end{bmatrix}
   * \f]
   */
  template <class Subdomain>
  class BDDCSolver // TODO: implement in terms of SymmetricPreconditioner
  {
  public:
    using Constraint = std::tuple<int,int,int>;

    /**
     * \brief Constructor.
     * \param subdomains a list of subdomains, initialized with an admissible solution approximation
     *                   and the corresponding residual
     * \param constraints a list of adjacent subdomain pairs with number of joint coarse constraints
     */
    BDDCSolver(std::vector<Subdomain>& subdomains,
               CoarseConstraints const& constraints);

    /**
     * \brief Performs a gradient step.
     *
     * \return the energy norm of the residual \f$ \|b-Ax\|_{A^{-1}} \f$.
     */
    double solve();

  private:
    using Domain = SharedMemoryDomain<Subdomain::components,typename Subdomain::Scalar>;
    std::vector<Subdomain>* subdomains;
    Domain domain;
    CoarseConstraints constraints;
    KKTSolver<Domain> coarseSolver;
  };

}

#endif

