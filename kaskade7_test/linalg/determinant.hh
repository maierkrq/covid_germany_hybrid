/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/en/numerik/software/kaskade-7.html               */
/*                                                                           */
/*  Copyright (C) 2015-2015 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DETERMINANT_HH
#define DETERMINANT_HH

namespace Kaskade
{
  /**
   * \cond internals
   */
  namespace DeterminantDetail
  {
    template <int dim, class Scalar>
    double determinantD1(Dune::FieldMatrix<Scalar,dim,dim> const& A, Dune::FieldMatrix<Scalar,dim,dim> const& dA)
    {
      if (dim==1)
        return A[0][0];

      if (dim==2)
        return A[0][0]*dA[1][1] + dA[0][0]*A[1][1] - A[0][1]*dA[1][0] - dA[0][1]*A[1][0];

      if (dim==3)
      {
        double ret = 0;

        // This is the rule of Sarrus, coded here for easy derivation of the derivative code.
        // for (int i=0; i<3; ++i)
        //   ret += A[0][i]*A[1][(i+1)%3]*A[2][(i+2)%3] - A[2][i]*A[1][(i+1)%3]*A[0][(i+2)%3];

        for (int i=0; i<3; ++i)
          ret += A[0][i]*A[1][(i+1)%3]*dA[2][(i+2)%3] + A[0][i]*dA[1][(i+1)%3]*A[2][(i+2)%3] + dA[0][i]*A[1][(i+1)%3]*A[2][(i+2)%3]
               - A[2][i]*A[1][(i+1)%3]*dA[0][(i+2)%3] - A[2][i]*dA[1][(i+1)%3]*A[0][(i+2)%3] - dA[2][i]*A[1][(i+1)%3]*A[0][(i+2)%3];

        return ret;
      }
    }

    template <int dim, class Scalar>
    double determinantD2(Dune::FieldMatrix<Scalar,dim,dim> const& A, Dune::FieldMatrix<Scalar,dim,dim> const& dA1, Dune::FieldMatrix<Scalar,dim,dim> const& dA2)
    {
      if (dim==1)
        return 0;

      if (dim==2)
        return dA2[0][0]*dA1[1][1] + dA1[0][0]*dA2[1][1] - dA2[0][1]*dA1[1][0] - dA1[0][1]*dA2[1][0];

      if (dim==3)
      {
        double ret = 0;
        for (int i=0; i<3; ++i)
          ret += A[0][i]*dA2[1][(i+1)%3]*dA1[2][(i+2)%3] + A[0][i]*dA1[1][(i+1)%3]*dA2[2][(i+2)%3] + dA1[0][i]*A[1][(i+1)%3]*dA2[2][(i+2)%3]
               + dA2[0][i]*A[1][(i+1)%3]*dA1[2][(i+2)%3] + dA2[0][i]*dA1[1][(i+1)%3]*A[2][(i+2)%3] + dA1[0][i]*dA2[1][(i+1)%3]*A[2][(i+2)%3]
               - A[2][i]*dA2[1][(i+1)%3]*dA1[0][(i+2)%3] - A[2][i]*dA1[1][(i+1)%3]*dA2[0][(i+2)%3] - dA1[2][i]*A[1][(i+1)%3]*dA2[0][(i+2)%3]
               - dA2[2][i]*A[1][(i+1)%3]*dA1[0][(i+2)%3] - dA2[2][i]*dA1[1][(i+1)%3]*A[0][(i+2)%3] - dA1[2][i]*dA2[1][(i+1)%3]*A[0][(i+2)%3];
        // TODO: does it pay off to factorize out common terms?
        return ret;
      }
    }

    template <int dim, class Scalar>
    double determinantD3(Dune::FieldMatrix<Scalar,dim,dim> const& A, Dune::FieldMatrix<Scalar,dim,dim> const& dA1,
                         Dune::FieldMatrix<Scalar,dim,dim> const& dA2, Dune::FieldMatrix<Scalar,dim,dim> const& dA3)
    {
      if (dim<=2)
        return 0;

      if (dim==3)
      {
        double ret = 0;
        for (int i=0; i<3; ++i)
          ret += dA3[0][i]*dA2[1][(i+1)%3]*dA1[2][(i+2)%3] + dA3[0][i]*dA1[1][(i+1)%3]*dA2[2][(i+2)%3] + dA1[0][i]*dA3[1][(i+1)%3]*dA2[2][(i+2)%3]
               + dA2[0][i]*dA3[1][(i+1)%3]*dA1[2][(i+2)%3] + dA2[0][i]*dA1[1][(i+1)%3]*dA3[2][(i+2)%3] + dA1[0][i]*dA2[1][(i+1)%3]*dA3[2][(i+2)%3]
               - dA3[2][i]*dA2[1][(i+1)%3]*dA1[0][(i+2)%3] - dA3[2][i]*dA1[1][(i+1)%3]*dA2[0][(i+2)%3] - dA1[2][i]*dA3[1][(i+1)%3]*dA2[0][(i+2)%3]
               - dA2[2][i]*dA3[1][(i+1)%3]*dA1[0][(i+2)%3] - dA2[2][i]*dA1[1][(i+1)%3]*dA3[0][(i+2)%3] - dA1[2][i]*dA2[1][(i+1)%3]*dA3[0][(i+2)%3];

        return ret;
      }
    }
  }
  /**
   * \endcond
   */

  /**
   * \brief A class for computing determinants and their derivatives.
   * \tparam dim the size of the matrix
   * \tparam Scalar the matrix entry type
   * \tparam byValue if true, the matrix is stored inside the determinant object, otherwise only a reference is stored
   */
  template <int dim, class Scalar=double, bool byValue=true>
  class Determinant
  {
    static_assert(0<dim && dim<=3,"determinant implemented only for dimensions 1-3");

  public:
    using Matrix = Dune::FieldMatrix<Scalar,dim,dim>;

    explicit Determinant(Matrix const& A_) : A(A_) {}

    Determinant(Determinant const&) = default;
    Determinant& operator=(Determinant const&) = default;

    double d0() const
    {
      return A.determinant();
    }

    double d1(Matrix const& dA) const
    {
      return DeterminantDetail::determinantD1(A,dA);
    }

    double d2(Matrix const& dA1, Matrix const& dA2) const
    {
      return DeterminantDetail::determinantD2(A,dA1,dA2);
    }

    double d3(Matrix const& dA1, Matrix const& dA2, Matrix const& dA3) const
    {
      return DeterminantDetail::determinantD3(A,dA1,dA2,dA3);
    }

  private:
    std::conditional_t<byValue,Matrix,Matrix const&> A;
  };

}

#endif // DETERMINANT_HH
