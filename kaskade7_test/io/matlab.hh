/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MATLAB_HH
#define MATLAB_HH

#include <fstream>
#include <string>

#include "linalg/triplet.hh"

namespace Kaskade
{
  
  /**
   * \ingroup IO
   * \brief Writes the given matrix to an executable matlab file.
   */
  template <class Entry, class Index>
  void writeToMatlab(NumaBCRSMatrix<Entry,Index> const& A, std::string const& basename, std::string const& outPath,int precision=16)
  {
    std::string fname = outPath + basename + ".m";
    std::ofstream f(fname.c_str());
    f.precision(precision);
    
    // Write header.
    f << "function A = " << basename << '\n';
    
    // write matrix in triplet format
    f << "data = [\n";
    int n = Entry::rows;
    int m = Entry::cols;
    
    for (auto row = A.begin(); row!=A.end(); ++row)
    {
      for (int i=0; i<n; ++i)
      {
        auto ridx = row.index()*n + i;
        for (auto col=row.begin(); col!=row.end(); ++col)
        {
          for (int j=0; j<m; ++j)
          {
            auto cidx = col.index()*m + j;
            f << ridx+1 << ' ' << cidx+1 << ' ' << (*col)[i][j] << std::endl;
          }
        }
      }
    }
    f << "];\nA = sparse(data(:,1),data(:,2),data(:,3)," << n*A.N() << "," << m*A.M() << ");\n";
  }
  
  /**
   * \ingroup IO
   * \brief Writes the given matrix and vector to an executable matlab file.
   */
  template <class Entry, class Index, class VEntry>
  void writeToMatlab(NumaBCRSMatrix<Entry,Index> const& A, Dune::BlockVector<VEntry> const& b, std::string const& basename, int precision=16)
  {
    std::string fname = basename + ".m";
    std::ofstream f(fname.c_str());
    f.precision(precision);
    
    // Write vector.
    f << "function [A,b] = " << basename << '\n'
    << " b = [\n";
    for (size_t i=0; i<b.N(); ++i)
      f << b[i] << std::endl;
    f << "];\n";
    
    // write matrix in triplet format
    f << "data = [\n";
    int n = Entry::rows;
    int m = Entry::cols;
    
    for (auto row = A.begin(); row!=A.end(); ++row)
    {
      for (int i=0; i<n; ++i)
      {
        auto ridx = row.index()*n + i;
        for (auto col=row.begin(); col!=row.end(); ++col)
        {
          for (int j=0; j<m; ++j)
          {
            auto cidx = col.index()*m + j;
            f << ridx+1 << ' ' << cidx+1 << ' ' << (*col)[i][j] << std::endl;
          }
        }
      }
    }
    f << "];\nA = sparse(data(:,1),data(:,2),data(:,3)," << n*A.N() << "," << m*A.M() << ");\n";
  }
  
  /**
   * \ingroup IO
   *
   * Writes the assembled matrix and right hand side of a
   * VariationalFunctionalAssembler to the given file. The contents is a
   * Matlab function that can be executed. The file extension .m is
   * appended automatically.
   *
   * \param[in] assembler a VariationalFunctionalAssembler containing the matrix and the right hand side
   * \param[in] basename the file name to be written to. A ".m" suffix is appended automatically.
   * \param[in] precision the number of valid decimal digits written for each scalar
   */
  template <class Assembler>
  void writeToMatlab(Assembler const& assembler, std::string const& basename, int precision=16)
  {
    std::string fname = basename + ".m";
    std::ofstream f(fname.c_str());
    f.precision(precision);

    f << "function [A,b] = " << basename << '\n'
        << " b = [\n";

    assembler.toSequence(0,Assembler::TestVariableSetDescription::noOfVariables,
        std::ostream_iterator<typename Assembler::Scalar>(f,"\n"));
    f << "\n];\ndata = [\n";

    typedef MatrixAsTriplet<typename Assembler::Scalar> Matrix;

    Matrix A = assembler.template get<Matrix>(false);
    
    for (size_t i=0; i<A.ridx.size(); ++i)
      f << A.ridx[i]+1 << ' ' << A.cidx[i]+1 << ' ' << A.data[i] << '\n';
    f << "];\nA = sparse(data(:,1),data(:,2),data(:,3)," << A.N() << "," << A.M() << ");\n";
  }

  template <class Assembler>
  void writeToMatlab(Assembler const& assembler, std::string const& path, std::string const& basename, int precision=16)
  {
    std::string fname = path + ".m";
    std::ofstream f(fname.c_str());
    f.precision(precision);

    f << "function [A,b] = " << basename << '\n'
        << " b = [\n";

    assembler.toSequence(0,Assembler::TestVariableSetDescription::noOfVariables,
        std::ostream_iterator<typename Assembler::Scalar>(f,"\n"));
    f << "\n];\ndata = [\n";

    typedef MatrixAsTriplet<typename Assembler::Scalar> Matrix;

    Matrix A = assembler.template get<Matrix>(false);
    
    for (size_t i=0; i<A.ridx.size(); ++i)
      f << A.ridx[i]+1 << ' ' << A.cidx[i]+1 << ' ' << A.data[i] << '\n';
    f << "];\nA = sparse(data(:,1),data(:,2),data(:,3)," << A.N() << "," << A.M() << ");\n";
  }
}


#endif
