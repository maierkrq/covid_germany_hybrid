/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FEM_DIFFOPS_MATERIALLAWS_HH
#define FEM_DIFFOPS_MATERIALLAWS_HH

#include <cassert>

#include "dune/common/fmatrix.hh"

#include "fem/diffops/elasto.hh"
#include "linalg/determinant.hh"

namespace Kaskade
{
  namespace Elastomechanics
  {
    /**
     * \ingroup diffopsElasto
     * \brief A namespace containing various material laws.
     * 
     * Hyperelastic material laws define a stored energy density \f$ W(E) \f$ in terms of the
     * Green-Lagrange strain tensor \f$ E \f$.
     */
    namespace MaterialLaws 
    {
      /**
       * \ingroup stationaryElasticity
       * \brief Base class for hyperelastic material laws, providing default implementations of 
       * the stress and the tangent stiffness tensor \f$ C \f$.
       * 
       * The stiffness tensor computation is injected by CRTP:
       * \code
       * class MyMaterial: public MaterialLawBase<MyMaterial> {...};
       * \endcode
       * 
       * \tparam dim the spatial dimension
       * \tparam Scalar the scalar field type (usually double)
       * \tparam MaterialLaw the actual material law implementation (derived class)
       */
      template <int dim_, class Scalar_, class MaterialLaw>
      class MaterialLawBase
      {
      public:
        using Scalar = Scalar_;
        static int const dim = dim_;
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>; 
        using VoigtTensor = Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2>;
        
        /**
         * \brief Returns the second Piola-Kirchhoff stress tensor \f$ S = W'(E) \f$ corresponding to the current strain \f$ E \f$.
         *
         * Note that \f$ S \f$ is always symmetric.
         */
        Tensor stress() const
        {
          Tensor sigma;
          
          for (int i=0; i<dim; ++i)
            for (int j=0; j<=i; ++j)
            {
              Tensor E(0); 
              E[i][j] = 1;
              sigma[i][j] = static_cast<MaterialLaw const*>(this)->d1(E);
            }
            
          for (int i=0; i<dim-1; ++i)
            for (int j=i+1; j<dim; ++j)
              sigma[i][j] = sigma[j][i];
            
          return sigma;
        }
        
        /**
         * \brief Returns the tangent stiffness tensor \f$ C(E) = W''(E)\f$ mapping strain tensor variations to stress tensor (2nd Piola-Kirchhoff) variations.
         * 
         * In this context, the symmetric matrices \f$ \sigma \f$ and \f$ \epsilon \f$ are interpreted as 
         * \f$d(d+1)/2\f$-vectors in Voigt notation as
         * returned by pack. Note that strain tensor off-diagonal entries appear with factor two in the vector.
         * 
         * Now, a matrix \f$ C\in\mathbb{R}^{(d+1)d/2\times(d+1)d/2} \f$ is returned such that for infinitesimal changes 
         * \f$\delta\epsilon\f$ the corresponding change in 
         * \f$ \sigma(\epsilon+\delta\epsilon) = \sigma(\epsilon)+\delta\sigma \f$ is given as 
         * \f$ \delta\sigma = C \delta\epsilon \f$. 
         * Note that due to the symmetry of the stored energy Hessian \f$ W''(E) \f$, the Voigt representation 
         * of the stiffness tensor is a symmetric matrix.
         * 
         * Note that is generic implementation usually exhibits a suboptimal performance. Specialized implementations 
         * can provide a significantly higher performance. For linear material laws, caching of 
         * the (constant) result might be an option.
         */
        VoigtTensor strainToStressMatrix() const
        {
          int const n = dim*(dim+1)/2;
          Dune::FieldMatrix<Scalar,n,n> C;
          
          // For the Voigt notation it holds that x^T C y = W''[unpack(x),unpack(y)]. 
          // Due to symmetry of C, we only compute the lower half directly...
          for (int i=0; i<n; ++i)
          {
            Dune::FieldVector<Scalar,n> x;
            x[i] = 1;
            auto X = unpack(x);
            
            for (int j=0; j<=i; ++j)
            {
              Dune::FieldVector<Scalar,n> y;
              y[j] = 1;
              auto Y = unpack(y);
              
              C[i][j] = static_cast<MaterialLaw const*>(this)->d2(X,Y);
            }
          }
          
          // ... and copy it to the upper half.
          for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
              C[i][j] = C[j][i];
            
          return C;
        }
      };
      
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup stationaryElasticity
       * \brief Adaptor for hyperelastic material laws, providing an easy way to formulate
       *        incompressible material laws in terms of the invariants of the Cauchy-Green
       *        strain tensor \f$ C = I+2E\f$.
       * 
       * This class essentially implements the composition \f$ W(E) = \hat W(I_1,I_2,I_3) \f$
       * with \f$ I_i = I_i(I+2E) \f$, computing derivatives by the chain rule. Note that most
       * hyperelastic material laws given in terms of the invariants are actually designed
       * for incompressible material, i.e. \f$ I_d = 1 \f$, and should be used directly only
       * in incompressible computations. For compressible hyperelasticity consider using
       * \ref CompressibleInvariantsMaterialLaw.
       * 
       * Model of HyperelasticMaterialLaw.
       * 
       * \tparam dim the spatial dimension
       * \tparam Scalar the scalar field type (usually double)
       * \tparam MaterialLaw the actual material law class, a model of InvariantsMaterialConcept
       */
      template <class Material>
      class InvariantsMaterialLaw: public MaterialLawBase<Material::dim,typename Material::Scalar,InvariantsMaterialLaw<Material>>
      {
      public:
        using Scalar = typename Material::Scalar;
        static int const dim = Material::dim;
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        /**
         * \brief Constructor.
         */
        template <class ... Args>
        InvariantsMaterialLaw(Args... args)
        : invariants(Tensor(0.0))
        , material(args...)
        {}
        
        /**
         * \name Hyperelastic material law interface
         * @{
         */
        
        /**
         * \brief Defines a new evaluation/linearization point.
         * \param e the strain tensor
         */
        void setLinearizationPoint(Tensor const& e)
        {
          invariants = ShiftedInvariants<dim,Scalar>(2.0*e);
          material.setLinearizationPoint(invariants.d0());
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(E) \f$.
         */
        Scalar d0() const
        {
          return material.d0();
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(E)E_1 \f$.
         */
        Scalar d1(Tensor const& e1) const
        {
          return 2*material.d1(invariants.d1(e1));
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(E)E_1 E_2 \f$.
         */
        Scalar d2(Tensor const& e1, Tensor const& e2) const
        {
          return 4*(material.d2(invariants.d1(e1),invariants.d1(e2)) + material.d1(invariants.d2(e1,e2)));
        }
        
        /**
         * @}
         */
      private:
        ShiftedInvariants<dim,Scalar> invariants;
        Material                      material;
      };
      
      /**
       * \ingroup stationaryElasticity
       * \brief Adaptor for hyperelastic material laws, providing an easy way to formulate compressible material laws in terms of 
       * the invariants of the isochoric part \f$ \bar C = C/(det C)^(1/d) \f$ of the Cauchy-Green strain tensor and a penalization
       * of the deviatoric part \f$ I_d \f$.
       * 
       * This class essentially implements the law \f$ W(E) = \bar W(\bar E) + D(I_d)\f$, where \f$ \bar E = \frac{E}{z} + \frac{1-z}{2z} I \f$
       * is the isochoric part of the strain tensor and \f$ z = I_d^{1/d} \f$, computing derivatives by the chain rule.
       * 
       * Model of HyperelasticMaterialLaw.
       *
       * \tparam dim the spatial dimension
       * \tparam Scalar the scalar field type (usually double)
       * \tparam MaterialLaw the actual material law class, a model of InvariantsMaterialConcept
       * \tparam DeviatoricPenalty
       */
      template <class Material, class DeviatoricPenalty>
      class CompressibleInvariantsMaterialLaw
      : public MaterialLawBase<Material::dim,typename Material::Scalar,CompressibleInvariantsMaterialLaw<Material,DeviatoricPenalty>>
      {
      public:
        using Scalar = typename Material::Scalar;
        static int const dim = Material::dim;
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        template <class ... Args>
        CompressibleInvariantsMaterialLaw(Tensor const& E, Args... args): bE(E), material(bE.d0(),args...), deviator(bE.determinant().d0()) {}
        
        /**
         * \brief Defines a new evaluation/linearization point.
         * \param e the strain tensor
         */
        void setLinearizationPoint(Tensor const& e)
        {
          bE = IsochoricGreenLagrangeTensor<dim,Scalar>(e);
          deviator = DeviatoricPenalty(bE.determinant().d0());
          material.setLinearizationPoint(bE.d0());
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(\epsilon) \f$.
         */
        Scalar d0() const
        {
          return material.d0() + deviator.d0();
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(\epsilon)\epsilon_1 \f$.
         */
        Scalar d1(Tensor const& e1) const
        {
          return material.d1(bE.d1(e1)) + deviator.d1(bE.determinant().d1(e1));
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(\epsilon)\epsilon_1\epsilon_2 \f$.
         */
        Scalar d2(Tensor const& e1, Tensor const& e2) const
        {
          return material.d2(bE.d1(e1),bE.d1(e2)) + material.d1(bE.d2(e1,e2))
                 + deviator.d2(bE.determinant().d1(e1),bE.determinant().d1(e2)) + deviator.d1(bE.determinant().d2(e1,e2));
        }
        
      private:
        IsochoricGreenLagrangeTensor<dim,Scalar> bE;
        InvariantsMaterialLaw<Material>          material;
        DeviatoricPenalty                        deviator;
      };
      
      // ---------------------------------------------------------------------------------------------------------
      // ---------------------------------------------------------------------------------------------------------
      
      
      /**
       * \ingroup stationaryElasticity
       * \brief The St. Venant-Kirchhoff material, foundation of linear elastomechanics.
       * 
       * The stored energy density is defined as
       * \f[ W(E) = \frac{\lambda}{2} (\mathop{\mathrm{tr}} E)^2 + \mu E:E +
       *            \alpha ( \det(I+2E)^{-1} +2 \mathop{\mathrm{tr}}E - 1), \f]
       * with a classical default value \f$ \alpha = 0 \f$.
       * 
       * The St. Venant-Kirchhoff material, characterized by any two of the five elastic moduli, is most
       * useful for linear elasticity. For finite strain elasticity, its use is not recommended, as it
       * is not polyconvex and, in its original form, does not preserve orientation of the deformation.
       *
       * For \f$ \alpha > 0 \f$, the material is orientation-preserving. For
       * \f$ \alpha > 0.0027(\lambda/2 + \mu) \f$, the material is Drucker-stable under uniaxial
       * compression.
       * 
       * This is a model of the \ref HyperelasticMaterialLaw concept.
       * 
       * \tparam dim the spatial dimension
       * \tparam Scalar the type of real numbers to use
       */
      template <int dim, class Scalar=double>
      class StVenantKirchhoff: public MaterialLawBase<dim,Scalar,StVenantKirchhoff<dim,Scalar>>
      {
        using Det = Determinant<dim,Scalar>;

      public:
        
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        /**
         * \brief Pretty useless default constructor.
         */
        StVenantKirchhoff() = default;
        
        /**
         * \brief Constructor.
         * \param moduli the material parameters
         * \param alpha the compression penalty. Should be larger than 0.0027(lambda/2+mu) for Drucker stability.
         */
        StVenantKirchhoff(ElasticModulus const& moduli, Scalar alpha_=0)
        : lambda(moduli.lame()), mu(moduli.shear()), alpha(alpha_)
        , detC(unitMatrix<Scalar,dim>())
        {
        }
        
        /**
         * \brief Constructor.
         * \param moduli the material parameters
         * \param e the initial Green-Lagrange strain tensor
         * 
         * This is equivalent to:
         * \code
         * StVenantKirchhoff material(moduli);
         * setLinearizationPoint(e);
         * \endcode
         */
        StVenantKirchhoff(ElasticModulus const& moduli, Tensor const& e_)
        : lambda(moduli.lame()), mu(moduli.shear()), e(e_), tre(trace(e))
        , detC(unitMatrix<Scalar,dim>()+2*e)
        {
        }
        
        /**
         * \name Hyperelastic material law interface
         * @{
         */
        
        /**
         * \brief Defines the linearization/evaluation point for subsequent calls to d0, d1, d2.
         */
        void setLinearizationPoint(Tensor const& e_)
        {
          e = e_;
          tre = trace(e);
          detC = Det(unitMatrix<Scalar,dim>()+2*e);
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(\epsilon) \f$.
         */
        Scalar d0() const
        {
          Scalar w = lambda*tre*tre/2 + mu*contraction(e,e);
          if (alpha > 0)
          {
            auto detc = detC.d0();
            if (detc <= 0)
              w = std::numeric_limits<Scalar>::infinity();
            else
              w += alpha*(1/detc-1 + 2*tre);
          }
          return w;
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(\epsilon)\epsilon_1 \f$.
         */
        Scalar d1(Tensor const& e1) const
        {
          Scalar dw = lambda*tre*trace(e1) + 2*mu*contraction(e,e1);
          if (alpha > 0)
            dw += alpha*(-std::pow(detC.d0(),-2)*detC.d1(e1)*2 + 2*trace(e1));
          return dw;
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(\epsilon)\epsilon_1\epsilon_2 \f$.
         */
        Scalar d2(Tensor const& e1, Tensor const& e2) const
        {
          Scalar ddw = lambda*trace(e1)*trace(e2) + 2*mu*contraction(e1,e2);
          if (alpha > 0)
            ddw += alpha*(  2*std::pow(detC.d0(),-3)*detC.d1(e1)*2*detC.d1(e2)*2
                          - std::pow(detC.d0(),-2)*detC.d2(e1,e2)*4 );
          return ddw;
        }

        /**
         * \brief Evaluates the third directional derivative \f$ W'''(\epsilon)[\epsilon_1,\epsilon_2,\epsilon_3] \f$.
         */
        Scalar d3(Tensor const& e1, Tensor const& e2, Tensor const& e3) const
        {
          if (alpha > 0)
            return alpha * ( -6*std::pow(detC.d0(),-4)*detC.d1(e1)*2*detC.d1(e2)*2*detC.d1(e3)*2
                             + 2*std::pow(detC.d0(),-3)*detC.d2(e1,e3)*4*detC.d1(e2)*2
                             + 2*std::pow(detC.d0(),-3)*detC.d1(e1)*2*detC.d2(e2,e3)*4
                             + 2*std::pow(detC.d0(),-3)*detC.d1(e3)*detC.d2(e1,e2)*4
                             - std::pow(detC.d0(),-2)*detC.d3(e1,e2,e3)*8 );
          else
            return 0;
        }

        /**
         * \brief Evaluates the fourth directional derivative \f$ W'''(\epsilon)[\epsilon_1,\epsilon_2,\epsilon_3,\epsilon_4] \f$.
         */
        Scalar d4(Tensor const& , Tensor const& , Tensor const&, Tensor const& ) const
        {
          // That's even easier for a quadratic energy...
          if (alpha > 0)
            abort(); // not yet implemented
          return 0;
        }

        /**
         * @}
         */
        
        private:
          Scalar lambda = 1;            // Lame parameter
          Scalar mu = 1;                // Lame parameter
          Scalar alpha = 0;             // compression penalty factor
          Tensor e = 0;                 // the strain tensor
          Scalar tre = 0;               // trace of strain tensor
          Determinant<dim,Scalar> detC; // determinant of Green-Lagrange tensor
      };
      
      // ---------------------------------------------------------------------------------------------------------
      
       /**
     * \ingroup stationaryElasticity
     * \brief Orthotropic linear material law.
     * 
     * This class defines the material law for linear orthotropic materials 
     * for vector-valued displacements \f$ u \f$, i.e. the variational functional 
     * \f$ u \mapsto \frac{1}{2} \epsilon : \mathcal{C}\epsilon \f$ with \f$ \epsilon = \frac{1}{2}(u_x^Tu_x + u_x + u_x^T) \f$.
     * In an appropriately rotated coordinate system, the stiffness matrix, representing the tensor \f$\mathcal{C}\f$ by interpreting 
     * \f$ \epsilon, \sigma \f$ as vectors, has the form 
     * \f[ C = \begin{bmatrix} * & * & \\ * & * & \\ & & * \end{bmatrix}, \quad C = \begin{bmatrix} * & * & * & & & \\ * & * & * & & & \\ * & * & * & & & \\
     *     & & & * & & \\ & & & & * & \\ & & & & & * \end{bmatrix} \f]
     * depending on the spatial dimension (see <a href="https://en.wikipedia.org/wiki/Orthotropic_material">Wikipedia</a>).
     */
    template <int dim, class Scalar=double>
    class OrthotropicLinearMaterial : public MaterialLawBase<dim,Scalar,OrthotropicLinearMaterial<dim,Scalar>>
    {
      
      using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
    
    public:
      
      /**
       * \brief Constructor.
       * \param orth The material coordinate system matrix. Has to be orthonormal.
       * \param mat  The material parameters.
       * 
       * The columns of the \a orth matrix are the material coordinate directions, i.e. a multiplication by this matrix transforms
       * from material coordinates to world coordinates.
       * 
       * \a mat is a matrix containing the material parameters in the following way:
       * \f[ \begin{bmatrix} E_1 & G_{12}  \\ \nu_{12} & E_2 \end{bmatrix}, \quad 
       *     \begin{bmatrix} E_1 & G_{12} & G_{13} \\ \nu_{12} & E_2 & G_{23} \\ \nu_{13} & \nu_{23} & E_3 \end{bmatrix}, \f]
       * where \f$ E_i \f$ is the elastic modulus in material coordinate direction \f$ i \f$, \f$ \nu_{ij} \f$ is the Poissons ratio
       * for contraction in \f$ j \f$ direction due to strain in \f$ i \f$ direction, and \f$ G_{ij} \f$ is the shear modulus in
       * \f$ ij \f$ plane.
       */
      OrthotropicLinearMaterial(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, Dune::FieldMatrix<Scalar,dim,dim> const& mat_)
      : orth(orth_), orthT(transpose(orth_)),  e(0)
      {
        Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> compliance(0);
        if (dim==3)
        {
          // see http://en.wikipedia.org/wiki/Orthotropic_material
          compliance[0][0] = 1 / mat_[0][0];             // 1/E_11
          compliance[0][1] = - mat_[1][0] / mat_[0][0];   // - nu_12 / E_11 = - nu_21 / E_22
          compliance[0][2] = - mat_[2][0] / mat_[0][0];   // - nu_13 / E_11 = - nu_31 / E_33
          compliance[1][0] = compliance[0][1];          // symmetry
          compliance[1][1] = 1 / mat_[1][1];
          compliance[1][2] = - mat_[2][1] / mat_[1][1];   // - nu_23 / E_22 = - nu_32 / E_33
          compliance[2][0] = compliance[0][2];          // symmetry
          compliance[2][1] = compliance[1][2];          // symmetry
          compliance[2][2] = 1 / mat_[2][2];
          
          compliance[3][3] = 1 / mat_[1][2];             // 1 / G_23
          compliance[4][4] = 1 / mat_[0][2];             // 1 / G_13
          compliance[5][5] = 1 / mat_[0][1];             // 1 / G_12
        } 
        else if (dim==2)
        {
          compliance[0][0] = 1 / mat_[0][0];             // 1/E_11;
          compliance[0][1] = - mat_[1][0] / mat_[0][0];   // - nu_12 / E_11 = - nu_21 / E_22
          compliance[1][0] = compliance[0][1];          // symmetry
          compliance[1][1] = 1 / mat_[1][1];             // 1/E_22
          
          compliance[2][2] = 1/mat_[0][1];               // 1 / G_12
        }
        
        // stiffness tensor is inverse of compliance tensor
        stiffness = compliance;
        stiffness.invert();
      }
      
      /**
       * \brief Constructor.
       * \param orth The material coordinate system matrix. Has to be orthonormal.
       * \param mat1   material parameters
       * \param mat2   material parameters
       * 
       * The columns of the \a orth matrix are the material coordinate directions, i.e. a multiplication by this matrix transforms
       * from material coordinates to world coordinates.
       * 
       * \a mat1 is a \f$ d \times d \f$ matrix containing the top left part of the stiffness matrix \f$ C \f$, and \a mat2
       * is a \f$ d(d-1)/2 \f$ vector containing the diagonal of the bottom right part.
       * \f[ C = \begin{bmatrix} \mathrm{mat1} &  \\  & \mathrm{diag}(\mathrm{mat2}) \end{bmatrix}, \f]
       */
      OrthotropicLinearMaterial(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, 
                            Dune::FieldMatrix<Scalar,dim,dim> const& mat1, Dune::FieldVector<Scalar,dim*(dim-1)/2> const& mat2)
      : orth(orth_), orthT(transpose(orth_)), stiffness(0), e(0)
      {
        for (int i=0; i<dim; ++i)
          for (int j=0; j<dim; ++j)
            stiffness[i][j] = mat1[i][j];
        for (int i=0; i<dim*(dim-1)/2; ++i)
          stiffness[dim+i][dim+i] = mat2[i];
      }
      
      /**
       * \brief Default constructor.
       * 
       * This initializes the material to an isotropic material with given material parameters (2D only).
       */
      OrthotropicLinearMaterial(ElasticModulus const& p)
      : orth(unitMatrix<Scalar,dim>()), orthT(unitMatrix<Scalar,dim>()), stiffness(0), e(0) 
      {
        Scalar lambda = p.lame(), mu = p.shear();
        if (dim==2)
        {
          stiffness[0][0] = stiffness[1][1] = 2*mu + lambda;
          stiffness[1][0] = stiffness[0][1] = lambda;
          stiffness[2][2] = mu;
        }
        if (dim==3)
          abort();
      }
      
      /**
       * \brief Default constructor.
       * 
       * This initializes the material to an isotropic material with Lame parameters set to 1.
       */
      OrthotropicLinearMaterial()
      : OrthotropicLinearMaterial(ElasticModulus()) {}
      
      /**
       * \brief Returns the stiffness tensor.
       *
       * The stiffness tensor \f$ C \f$ is represented in matrix form for Voigt notation (i.e. the mapping from Voigt strain vector to Voigt stress vector). 
       * In 3D this reads
       * \f[ \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} \end{bmatrix} = C
       *     \begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ 2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12} \end{bmatrix} \f]
       */
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> const& stiffnessMatrix() const { return stiffness; }
      
      /**
       * \brief Defines the GreenLagrangeTensor around which to linearize.
       */
      void setLinearizationPoint(Tensor const& e_)
      {
        e = e_;
      }
      
      /**
       * \brief Computes the elastic energy \f$ \frac{1}{2} \epsilon(u_x) : \mathcal{C}\epsilon(u_x) \f$.
       */
      Scalar d0() const
      {
        auto eps = pack(orthT*e*orth);
        auto sigma = stiffness * eps;
        return 0.5 * (sigma*eps); 
      }
      
      /**
       * \brief Computes the first derivative of the elastic energy in the direction \f$ e1\f$.
       * 
       */
      Scalar d1(Tensor const& e1) const
      {
        auto eps  = pack(orthT*e*orth);
        auto eps1 = pack(orthT*e1*orth);
        auto sigma = stiffness * eps1;
        return sigma*eps; 
      }
      
      /**
       * \brief Computes the second derivative of the elastic energy in directions
       * \f$ e1, e2\f$.
       */
      Scalar d2(Tensor const& e1, Tensor const &e2) const
      {
        auto eps1 = pack(orthT*e1*orth);
        auto eps2 = pack(orthT*e2*orth);
        auto sigma = stiffness * eps2;
        return sigma*eps1; 
      }
      
      Scalar d3(Tensor const& e1, Tensor const& e2, Tensor const &e3) const
      {
        return 0.0;
      }
      
      Scalar d4(Tensor const& e1, Tensor const& e2, Tensor const &e3, Tensor const &e4) const
      {
        return 0.0;
      }
      
    private:
      Dune::FieldMatrix<Scalar,dim,dim> orth; // the orthonormal matrix mapping material to global coordinates
      Dune::FieldMatrix<Scalar,dim,dim> orthT; // the transposed orthonormal matrix mapping material to global coordinates (for convenience)
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> stiffness;
      Dune::FieldMatrix<Scalar,dim,dim> e; // the linearization point 
    };
      // ---------------------------------------------------------------------------------------------------------
      
      template <int dim, class Scalar=double>
    class OrthotropicNonLinearMaterial : public MaterialLawBase<dim,Scalar,OrthotropicLinearMaterial<dim,Scalar>>
    {
      
      using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
    
    public:
      
      /**
       * \brief Constructor.
       * \param orth The material coordinate system matrix. Has to be orthonormal.
       * \param mat  The material parameters.
       * 
       * The columns of the \a orth matrix are the material coordinate directions, i.e. a multiplication by this matrix transforms
       * from material coordinates to world coordinates.
       * 
       * \a mat is a matrix containing the material parameters in the following way:
       * \f[ \begin{bmatrix} E_1 & G_{12}  \\ \nu_{12} & E_2 \end{bmatrix}, \quad 
       *     \begin{bmatrix} E_1 & G_{12} & G_{13} \\ \nu_{12} & E_2 & G_{23} \\ \nu_{13} & \nu_{23} & E_3 \end{bmatrix}, \f]
       * where \f$ E_i \f$ is the elastic modulus in material coordinate direction \f$ i \f$, \f$ \nu_{ij} \f$ is the Poissons ratio
       * for contraction in \f$ j \f$ direction due to strain in \f$ i \f$ direction, and \f$ G_{ij} \f$ is the shear modulus in
       * \f$ ij \f$ plane.
       */
      OrthotropicNonLinearMaterial(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, Dune::FieldMatrix<Scalar,dim,dim> const& mat_)
      : orth(orth_), orthT(transpose(orth_)),  e(0), a(0), G(0)
      {
		 


		Scalar D = 0.75 * (2/(mat_[0][0]*mat_[1][1]) + 2/(mat_[1][1]*mat_[2][2]) + 2/(mat_[2][2]*mat_[0][0]) 
		               - 1/(mat_[0][0]*mat_[0][0]) - 1/(mat_[1][1]*mat_[1][1]) - 1/(mat_[2][2]*mat_[2][2]));

				//Materialparameters aij
a[0][0] = 1/D * (2/(3*mat_[1][1]) + 2/(3*mat_[2][2]) - 1/(3*mat_[0][0]));    
a[0][1] = 1/D * (1/(6*mat_[0][0]) + 1/(6*mat_[1][1]) - 5/(6*mat_[2][2])); 	
a[0][2] = 1/D * (1/(6*mat_[0][0]) + 1/(6*mat_[2][2]) - 5/(6*mat_[1][1]));
		
a[1][0] = 1/D * (1/(6*mat_[1][1]) + 1/(6*mat_[0][0]) - 5/(6*mat_[2][2]));  	
a[1][1] = 1/D * (2/(3*mat_[2][2]) + 2/(3*mat_[0][0]) - 1/(3*mat_[1][1])); 	
a[1][2] = 1/D * (1/(6*mat_[1][1]) + 1/(6*mat_[2][2]) - 5/(6*mat_[0][0]));
		
a[2][0] = 1/D * (1/(6*mat_[2][2]) + 1/(6*mat_[0][0]) - 5/(6*mat_[1][1]));   
a[2][1] = 1/D * (1/(6*mat_[2][2]) + 1/(6*mat_[1][1]) - 5/(6*mat_[0][0]));	
a[2][2] = 1/D * (2/(3*mat_[1][1]) + 2/(3*mat_[0][0]) - 1/(3*mat_[2][2]));
        
        //std::cout << a << "\n";
        
		//Materialparameters Gij
								G[0][1] = mat_[0][1]; 	G[0][2] = mat_[0][2]; 
		G[1][0] = mat_[0][1]; 							G[1][2] = mat_[1][2]; 
		G[2][0] = mat_[0][2]; 	G[2][1] = mat_[1][2]; 
		     
        // stiffness tensor see http://en.wikipedia.org/wiki/Orthotropic_material
        stiffness[0][0] = a[0][0] ; stiffness[0][1] = a[0][1] ; stiffness[0][2] = a[0][2] ;
        stiffness[1][0] = a[1][0] ; stiffness[1][1] = a[1][1] ; stiffness[1][2] = a[1][2] ;
        stiffness[2][0] = a[2][0] ; stiffness[2][1] = a[2][1] ; stiffness[2][2] = a[2][2] ;
        
        stiffness[3][3] = mat_[1][2];
        stiffness[4][4] = mat_[0][2];
        stiffness[5][5] = mat_[0][1];         
          
        //std::cout << stiffness << "\n";  
          
        // Structural tensor 
        x[0] =1;x[1] =0;x[2] =0;  
        y[0] =0;y[1] =1;y[2] =0;  
        z[0] =0;z[1] =0;z[2] =1;      
                        
        L[0] = outerProduct(orth*x,orth*x);
        L[1] = outerProduct(orth*y,orth*y);
        L[2] = outerProduct(orth*z,orth*z);
        
        /*
        std::cout << "L0"<< "\n" <<L[0] << "\n" << std::endl;
        std::cout << "L1"<< "\n" <<L[1] << "\n" << std::endl;
        std::cout << "L2"<< "\n" <<L[2] << "\n" << std::endl;
		*/
        
       }
      
      /**
       * \brief Constructor.
       * \param orth The material coordinate system matrix. Has to be orthonormal.
       * \param mat1   material parameters
       * \param mat2   material parameters
       * 
       * The columns of the \a orth matrix are the material coordinate directions, i.e. a multiplication by this matrix transforms
       * from material coordinates to world coordinates.
       * 
       * \a mat1 is a \f$ d \times d \f$ matrix containing the top left part of the stiffness matrix \f$ C \f$, and \a mat2
       * is a \f$ d(d-1)/2 \f$ vector containing the diagonal of the bottom right part.
       * \f[ C = \begin{bmatrix} \mathrm{mat1} &  \\  & \mathrm{diag}(\mathrm{mat2}) \end{bmatrix}, \f]
       */
      OrthotropicNonLinearMaterial(Dune::FieldMatrix<Scalar,dim,dim> const& orth_, 
                            Dune::FieldMatrix<Scalar,dim,dim> const& mat1, Dune::FieldVector<Scalar,dim*(dim-1)/2> const& mat2)
      : orth(orth_), orthT(transpose(orth_)), stiffness(0), e(0)
      {
        for (int i=0; i<dim; ++i)
          for (int j=0; j<dim; ++j)
            stiffness[i][j] = mat1[i][j];
        for (int i=0; i<dim*(dim-1)/2; ++i)
          stiffness[dim+i][dim+i] = mat2[i];
      }
      
      /**
       * \brief Default constructor.
       * 
       * This initializes the material to an isotropic material with given material parameters (2D only).
       */
      OrthotropicNonLinearMaterial(ElasticModulus const& p)
      : orth(unitMatrix<Scalar,dim>()), orthT(unitMatrix<Scalar,dim>()), stiffness(0), e(0) 
      {
        Scalar lambda = p.lame(), mu = p.shear();
        if (dim==2)
        {
          stiffness[0][0] = stiffness[1][1] = 2*mu + lambda;
          stiffness[1][0] = stiffness[0][1] = lambda;
          stiffness[2][2] = mu;
        }
        if (dim==3)
          abort();
      }
      
      /**
       * \brief Default constructor.
       * 
       * This initializes the material to an isotropic material with Lame parameters set to 1.
       */
      OrthotropicNonLinearMaterial()
      : OrthotropicNonLinearMaterial(ElasticModulus()) {}
      
      /**
       * \brief Returns the stiffness tensor.
       *
       * The stiffness tensor \f$ C \f$ is represented in matrix form for Voigt notation (i.e. the mapping from Voigt strain vector to Voigt stress vector). 
       * In 3D this reads
       * \f[ \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} \end{bmatrix} = C
       *     \begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ 2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12} \end{bmatrix} \f]
       */
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> const& stiffnessMatrix() const { return stiffness; }
      
      /**
       * \brief Defines the GreenLagrangeTensor around which to linearize.
       */
      void setLinearizationPoint(Tensor const& e_)
      {
        e = e_;
      }
      
      /**
       * \brief Computes the elastic energy \f$ \frac{1}{2} \epsilon(u_x) : \mathcal{C}\epsilon(u_x) \f$.
       */
      Scalar d0() const
      {
			
	    Scalar shear = 0 , eps = 0;
	    Dune::FieldMatrix<Scalar, 3 ,3> A(0),B(0);
	
		for(unsigned int i = 0 ; i < dim ; ++i)
		{
			for(unsigned int j = 0 ; j < dim ; ++j)
			{
				//eps += a[i][j]*contraction(transpose(e),L[i])*contraction(transpose(e),L[j]);
				eps += a[i][j]*trace(e*(L[i]))*trace(e*(L[j]));
			}
		}		
				
		for(unsigned int i = 0 ; i < dim ; ++i)
		{
			auto A = e*L[i];
			for(unsigned int j = 0 ; j < dim ; ++j)
			{
				if(i==j)
				{continue;}		
				
				auto B = e.rightmultiplyany(L[j]);
				
				shear += G[i][j]*trace(A.rightmultiplyany(B));
			}
		}			
		

        return 0.5 * eps + shear;
      }
      
      /**
       * \brief Computes the first derivative of the elastic energy in the direction \f$ e1\f$.
       * 
       */
      Scalar d1(Tensor const& e1) const
      {
	    Scalar shear = 0 , eps = 0;
	    Dune::FieldMatrix<Scalar, 3 ,3> A(0), B(0), C(0), D(0);
		  
		for(unsigned int i = 0 ; i < dim ; ++i)
		{
			for(unsigned int j = 0 ; j < dim ; ++j)
			{
				eps += a[i][j]*(trace(e*(L[i]))*trace(e1*(L[j])) 
							   +trace(e1*(L[i]))*trace(e*(L[j])));
			}
		}	
		
		for(unsigned int i = 0 ; i < dim ; ++i)
		{
			auto A = e.rightmultiplyany(L[i]);
			auto C = e1.rightmultiplyany(L[i]);
			
			for(unsigned int j = 0 ; j < dim ; ++j)
			{
				if(i==j)
				{continue;}
								
				auto B = e1.rightmultiplyany(L[j]);
				auto D = e.rightmultiplyany(L[j]);
				
				shear += G[i][j]*(trace(A.rightmultiplyany(B))+trace(C.rightmultiplyany(D)));
			}
		}

		
		

        return 0.5 * eps + shear;
      }
      
      /**
       * \brief Computes the second derivative of the elastic energy in directions
       * \f$ e1, e2\f$.
       */
      Scalar d2(Tensor const& e1, Tensor const &e2) const
      {
 	    
 	    Scalar shear = 0 , eps = 0;
 	    Dune::FieldMatrix<Scalar, 3 ,3> A(0),B(0),C(0),D(0);
		  
		for(unsigned int i = 0 ; i < dim ; ++i)
		{
			for(unsigned int j = 0 ; j < dim ; ++j)
			{
				eps += a[i][j]*(trace(e2*L[i])*trace(e1*L[j]) 
							   +trace(e1*L[i])*trace(e2 *L[j]));
			}
		}	
		
		for(unsigned int i = 0 ; i < dim ; ++i)
		{
			auto A = e2*L[i];
			auto C = e1*L[i];
			
			for(unsigned int j = 0 ; j < dim ; ++j)
			{
				if(i==j)
				{continue;}				
				
				auto B = e1*L[j];
				auto D = e2*L[j];
				
				shear += G[i][j]*(trace(A.rightmultiplyany(B))+trace(C.rightmultiplyany(D)));
			}
		}	

        return 0.5 * eps + shear;
      }
      
      Scalar d3(Tensor const& e1, Tensor const& e2, Tensor const &e3) const
      {
        return 0.0;
      }
      
      Scalar d4(Tensor const& e1, Tensor const& e2, Tensor const &e3, Tensor const &e4) const
      {
        return 0.0;
      }
      
    private:
      Dune::FieldMatrix<Scalar,dim,dim> orth; // the orthonormal matrix mapping material to global coordinates
      Dune::FieldMatrix<Scalar,dim,dim> orthT; // the transposed orthonormal matrix mapping material to global coordinates (for convenience)
      Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> stiffness;
      Dune::FieldMatrix<Scalar,dim,dim> e; // the linearization point 
      Dune::FieldMatrix<Scalar,dim,dim> a; // Materialparameters
	  Dune::FieldMatrix<Scalar,dim,dim> G; // Materialparameters
	  Dune::FieldVector< Dune::FieldMatrix<Scalar,dim,dim> , 3 > L;//Structural tensor
	  Dune::FieldVector< Scalar , dim > x,y,z;

    };
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup stationaryElasticity
       * \brief Mooney-Rivlin material law formulated in terms of the (shifted) invariants \f$ i_1, i_2, i_3 \f$ of the 
       * doubled Green-Lagrange strain tensor \f$ 2E \f$.
       * 
       * The Mooney-Rivlin hyperelastic material energy for incompressible materials is defined in terms of the 
       * Cauchy-Green strain tensor \f$ C = I+2E \f$ as
       * \f[ W = C_1(I_1-3) + C_2(I_2-3) = C_1 i_1 + C_2 i_2. \f]
       * The material parameters \f$ 2(C_1+C_2) \f$ correspond to the shear modulus \f$ \mu \f$ of linear elasticity. 
       * 
       * Here we extend the material law directly to compressible materials by adding a volumetric part:
       * \f[ W = C_1(I_1-3) + C_2(I_2-3) + C_3(I_d+I_d^{-1}-2) = C_1 i_1 + C_2 i_2 + C_3\frac{i_d^2}{1+i_d}. \f]
       * 
       * For 2D problems, \f$ C_2=0 \f$ is used independent of what is specified in the constructor, such that for 2D
       * the Mooney-Rivlin and the Neo-Hookean model coincide.
       * 
       * The implementation is numerically stable in the vicinity of the reference configuration.
       * 
       * Use this with the \ref IncompressibleInvariantsMaterialLaw adaptor for creating a hyperelastic material law.
       * If \f$ C_3\ne 0 \f$ this will also work for compressible situations. Alternatively,  
       * the \ref CompressibleInvariantsMaterialLaw adaptor can be used for compressible materials.
       * 
       * Model of InvariantsMaterialConcept.
       *
       * \see Kaskade::InvariantsMaterialLaw
       */
      template <int dimension>
      class MooneyRivlin
      {
      public:
        /**
         * \brief The scalar field type (usually double).
         */
        using Scalar = double;
        
        /**
         * \brief The spatial dimension (2 or 3).
         */
        static int const dim = dimension;
        
        /**
         * \brief A d-dimensional vector type.
         */
        using Invariants = Dune::FieldVector<Scalar,dim>;
        
        /**
         * \brief Constructor.
         * \param c1 material parameter (factor for I1)
         * \param c2 material parameter (factor for I2), ignored for 2D problems
         * \param c3 material parameter (factor for volumetric penalty)
         */
        MooneyRivlin(double c1_, double c2_, double c3_=0.0)
        : c1(c1_), c2(dim==3?c2_:0), c3(c3_)
        {}
        
        /**
         * \brief Sets a new linearization point.
         */
        void setLinearizationPoint(Invariants const& i)
        {
          i1 = i[0];
          i2 = i[1];
          id = i[dim-1];
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(I) \f$.
         */
        Scalar d0() const
        {
          return c1*i1 + c2*i2 + c3*id*id/(1+id);
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(I)d_1I \f$.
         */
        Scalar d1(Invariants const& di1) const
        {
          return c1*di1[0] + c2*di1[1] + c3*id*(2+id)/square(1+id)*di1[dim-1];
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(I)\d_1I d_2I \f$.
         */
        Scalar d2(Invariants const& di1, Invariants const& di2) const
        {
          return 2 / power(1+id,3) *di1[dim-1]*di2[dim-1];
        }
        
      private:
        double i1, i2, id, c1, c2, c3;
      };
      
      // ------------------------------------------------------------------------------------------

      /**
       * \ingroup stationaryElasticity
       * \brief Neo-Hookean material law formulated in terms of the shifted invariants \f$ i_1, i_2, i_3 \f$ 
       * of the doubled Green-Lagrange strain tensor \f$ 2E \f$.
       * 
       * The Neo-Hookean hyperelastic material energy for incompressible materials is defined 
       * in terms of the Cauchy-Green strain tensor \f$ C = I+2E\f$ as
       * \f[ W = \frac{G}{2} (I_1-d) = \frac{G}{2} i_1. \f]
       * The material parameter \f$ G \f$ corresponds to the shear modulus of linear elasticity. 
       * The material law is a special case of the more general MooneyRivlin material law.
       * 
       * The implementation is numerically stable in the vicinity of the reference configuration.
       * 
       * Use the \ref InvariantsMaterialLaw adaptor for creating a hyperelastic material law, 
       * or the \ref CompressibleInvariantsMaterialLaw adaptor in case of compressible materials.
       * 
       * Model of InvariantsMaterialConcept.
       * 
       * \see \ref InvariantsMaterialLaw
       */
      template <int dimension>
      class NeoHookean: public MooneyRivlin<dimension>
      {
      public:
        using typename MooneyRivlin<dimension>::Invariants;

        /**
         * \brief Constructor.
         * \param moduli the elastic moduli of the material 
         */
        NeoHookean(ElasticModulus const& moduli): MooneyRivlin<dimension>(moduli.shear()/2,0)
        {}
      };

      // ------------------------------------------------------------------------------------------
      
      /**
       * \ingroup stationaryElasticity
       * \brief Blatz-Ko material law for rubber foams in terms of the shifted invariants \f$ i_1, i_2, i_3 \f$ of 
       * the doubled Green-Lagrange strain tensor \f$ 2E \f$.
       * 
       * In contrast to many other material laws given in terms of the invariants, the Blatz-Ko material model 
       * for foams is explicitly formulated for compressible materials. The stored energy density is 
       * \f[ W = \frac{\mu}{2} \left(i_1 + \frac{2}{a} ( (1+i_3)^{-a/2}-1 ) \right), \f] with \f$ \mu \f$ the shear modulus
       * and \f$ a = 2\nu / (1-2\nu) \f$ related to Poisson's ratio \f$ \nu > 0 \f$.
       * 
       * Use the \ref InvariantsMaterialLaw adaptor for creating a hyperelastic material law.
       * 
       * Model of InvariantsMaterialConcept.
       * 
       * \see \ref InvariantsMaterialLaw
       */
      template <int dimension>
      class BlatzKo
      {
      public:
        /**
         * \brief The scalar field type (usually double).
         */
        using Scalar = double;
        
        /**
         * \brief The spatial dimension (2 or 3).
         */
        static int const dim = dimension;
        
        /**
         * \brief A d-dimensional vector type.
         */
        using Invariants = Dune::FieldVector<Scalar,dim>;
        
        /**
         * \brief Constructor.
         * \param moduli The elastic moduli. The poisson ratio has to be positive (\f$ \nu > 0 \f$).
         */
        BlatzKo(ElasticModulus const& moduli)
        : mu(moduli.shear()), a(2*moduli.poisson()/(1-2*moduli.poisson())), pstable(1,0)
        {
          assert(a>0);
        }
        
        /**
         * \brief Sets a new linearization point.
         */
        void setLinearizationPoint(Invariants const& i)
        {
          i1 = i[0];
          i3 = i[dimension-1];
          pstable = Pstable(-a/2,i3);
        }
        
        /**
         * \brief Evaluates the stored energy density \f$ W(I) \f$.
         */
        Scalar d0() const
        {
          return mu*(i1 + 2*pstable.d0()/a)/2;
        }
        
        /**
         * \brief Evaluates the first directional derivative \f$ W'(I)d_1I \f$.
         */
        Scalar d1(Invariants const& di1) const
        {
          return mu*(di1[0] + 2*pstable.d1(di1[dimension-1])/a)/2;
        }
        
        /**
         * \brief Evaluates the second directional derivative \f$ W''(I)\d_1I d_2I \f$.
         */
        Scalar d2(Invariants const& di1, Invariants const& di2) const
        {
          return mu*pstable.d2(di1[dimension-1],di2[dimension-1])/a;
        }
        
      private:
        double i1, i3, mu, a;
        Pstable pstable;
      };
      
      
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup viscoPlasticity
       * \brief An adaptor for using hyperelastic stored energies for viscoplasticity.
       * 
       * The stored energy density for viscoplasticity is defined as \f$ W^{vp}(E) = W(E-E^{vp}). \f$ 
       * The viscoplastic strain \f$ E^{vp} \f$ acts as an offset, such that the material is stress-free
       * in a deformed configuration. In time-dependent viscoplastic models, the viscoplastic strain 
       * usually evolves according to some law depending on the current stress.
       * 
       * Model of HyperelasticMaterialLaw.
       * 
       * \tparam HyperelasticEnergy the material law of elastic response
       */
      template <class HyperelasticEnergy>
      class ViscoPlasticEnergy: public HyperelasticEnergy
      {
      public:
        using Scalar = typename HyperelasticEnergy::Scalar;
        static int const dim = HyperelasticEnergy::dim;
        
        using Tensor = Dune::FieldMatrix<Scalar,dim,dim>;
        
        /**
         * \brief Default constructor.
         */
        ViscoPlasticEnergy() = default;
        
        /**
         * \brief Constructor.
         * \param evp The viscoplastic stain.
         * The remaining parameters are forwarded to the hyperelastic energy density constructor.
         */
        template<typename... Args>
        ViscoPlasticEnergy(Tensor const& evp_, Args&&... args): HyperelasticEnergy(std::forward<Args>(args)...), evp(evp_)
        {}
        
        /**
         * \brief Sets a new linearization point.
         */
        void setLinearizationPoint(Tensor const& e)
        {
          HyperelasticEnergy::setLinearizationPoint(e-evp);
        }
                
      private:
        Tensor             evp;
      };
      
      // ---------------------------------------------------------------------------------------------------------
      
      /**
       * \ingroup viscoPlasticity
       * \brief Computes the von Mises equivalent stress \f$ \sigma_v \f$.
       * 
       * The equivalent von Mises stress is defined as 
       * \f[ 2\sigma_v^2 = (\sigma_{11}-\sigma_{22})^2 + (\sigma_{22}-\sigma_{33})^2 + (\sigma_{33}-\sigma_{11})^2 + 6 (\sigma_{23}^2+\sigma_{31}^2 + \sigma_{12}^2) \f]
       * in terms of the Cauchy stress \f$ \sigma \f$.
       * 
       * \tparam dim the spatial dimension. For d=2, plane stress is assumed.
       * \tparam Scalar the scalar field type of the stress tensor, usually double
       * 
       * Assuming there are finite element functions lambda and mu giving the material parameters, and u the displacement, a scalar finite element function sv of
       * von Mises equivalent stress values can be obtained as follows:
       * \code
       * interpolateGlobally<PlainAverage>(sv,makeFunctionView(u.space(), [&] (auto const& evaluator)
       * {
       *   ElasticModulus modulus(lambda.value(evaluator.cell(),evaluator.xloc()),
       *                          mu.value(evaluator.cell(),evaluator.xloc()));
       *   HyperelasticVariationalFunctional<StVenantKirchhoff<dim>,GreenLagrangeTensor<double,dim>>  material(modulus);
       *   material.setLinearizationPoint(u.derivative(evaluator));
       *   return Dune::FieldVector<double,1>(vonMisesStress(material.cauchyStress()));
       * }));
       * \endcode
       */
      template <int dim, class Scalar>
      Scalar vonMisesStress(Dune::FieldMatrix<Scalar,dim,dim> const& stress);

      // ---------------------------------------------------------------------------------------------------------

      /**
       * \ingroup viscoPlasticity
       * \brief A simple Duvaut-Lions flow rule with \f$ J_2 \f$ (von Mises) yield surface.
       * 
       * The Duvaut-Lions model of viscoplasticity defines the flow rate, i.e. the time derivative of the viscoplastic strain,
       * as \f[ \dot \epsilon^{\rm vp} = \tau^{-1} C^{-1} (\sigma - P\sigma), \f]
       * where \f$ P \f$ is the closest point projector onto the admissible stress state. In \f$ J_2 \f$ viscoplasticity, the 
       * admissible stresses are characterized by \f$ \|\sigma\| \le \sqrt{2/3}\,\sigma_Y \f$, and the boundary of that set 
       * is known as von Mises yield surface.
       * 
       * \param cauchyStress the true stress tensor
       * \param C the strain to stress matrix in Voigt notation
       * \param tau the relaxation time (>0)
       * \param sigmaY the yield strength (>=0)
       * \return the flow rate \f$ \dot \epsilon^{\rm vp} \f$ in Voigt notation
       */
      template <int dim, class Scalar=double>
      Dune::FieldVector<Scalar,dim*(dim+1)/2> 
      duvautLionsJ2Flow(Dune::FieldMatrix<Scalar,dim,dim> const& cauchyStress, Dune::FieldMatrix<Scalar,dim*(dim+1)/2,dim*(dim+1)/2> const& C, double tau, double sigmaY);
      
      // ---------------------------------------------------------------------------------------------------------
      
    }
  }
}

#endif
