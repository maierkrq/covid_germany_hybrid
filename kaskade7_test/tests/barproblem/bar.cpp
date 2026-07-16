/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2009 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <cmath>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include "dune/grid/geometrygrid.hh"

#include "fem/assemble.hh"
#include "fem/norms.hh"
#include "fem/embedded_errorest.hh"
#include "fem/istlinterface.hh"
#include "fem/functional_aux.hh"
//#include "fem/hierarchicspace.hh"
#include "fem/lagrangespace.hh"
#include "linalg/direct.hh"
#include "io/vtk.hh"

#include "utilities/gridGeneration.hh"
#include "utilities/kaskopt.hh"

/**
  * \defgroup tests Tests
  * \brief Classes and functions that test whether certain parts of Kaskade work as they should.
  */
 /**
 * @file
 * @ingroup tests
 * @brief  Test with 2-dimensionale linear elasticity equation.
 * 
 * Testprogram for Kaskade: tests wether the error between calculated and exact solution behaves as expected. 
 * This test was built by giving an exact solution u on \Omega:= \{ (x_1,x_2) \in R³ | 0<=x<=4, 0<=x_2<=1 \} (Rectangle, a=4, b=1} and constructing a corresponding BVP.
 *  
 * The solution is given by 
 * 
 * u(x)=1/ \mu * g(x_1)*e_2, 
 * 
 * with
 *  
 * g(x)=sin(0.5* \pi*x).
 * 
 * We get the BVP:
 * 
 * - \bigtriangledown \cdot \sigma(u)(x) = f(x) = -g_2''(x_1)*e_2-g_3''(x_1)e_3 on \Omega
 *  \sigma(u)(x)* \eta = (g_2''(x_1)* \eta_2)*e_1 on \Theta_1 := \{ x \in \Omega | x_2=0 or x_2=1 \}
 *  u= u0=0 on \Theta_2 := \{ x \in \Omega | x_1=0 or x_1=4 \}
 * 
 * Boundary conditions and f are adapted in elasto.hh. 
 * 
 * Uses a direct solver to solve linear problem.
 */

using namespace Kaskade;
#include "bar.hh"
struct ExactSolution //exact solution of the  BVP
{
   using Scalar = double;
   static int const components = 2;
   using ValueType = Dune::FieldVector<Scalar,components>;

   template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }

   template <class Cell>
   ValueType value(Cell const& cell,Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
   {
     Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
     double pi=3.141592653589793238462;
     ValueType y; y[0]=0; y[1]=0.125*sin(0.5*pi*x[0]); 
     return y;
   }
};



int main(int argc, char *argv[])
{
  std::cout << "Start Bar program" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int maxRefinements, maxOrder, result, quick, verbosity;
  double length;
  bool valid=true;

  if (getKaskadeOptions(argc,argv,Options
    ("refinements",      maxRefinements,        5,          "maximum number of uniform grid refinements")
    ("order",            maxOrder,              4,          "maximum finite element ansatz order")
    ("result",           result,                0,          "set 1 if testResult.text should be written")
    ("verbosity",        verbosity,             1,          "set 1 if testResult.text should be written")
    ("length",			     length,			          4.0,	     	"length of the cylinder")  //shouldn't be changed for test to work
    ("quick",            quick,                 0,			" set 1 for short test")
    ))
    return 1;


  

 // THREE-dimensional space: dim=3
  constexpr int dim=2;        
  using Grid = Dune::UGGrid<dim>;
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView>>;
  using Vector = Dune::FieldVector<double, dim>;
  
  if (quick)
  {
	maxRefinements=4;
	maxOrder=3;
  }
  else
  {
    std::cout << "Test can take a while, use '--quick 1' for shorter test\n";
  }  
    
  int orders[4]={1,2,3,4}; //max refinements possible for corresponding order
  
  // logs of meshsizes for corresponding refinements (0 to 4)
  double hs[6]={1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125}; 
  for (int m=0; m<6; m++)
  {
    hs[m]=std::log(hs[m]);
  }
  
  double preErrors[4][6]={ {0.029, 0.014, 0.0044, 0.0013, 0.00032, 8.5e-5},
                       {0.0022, 0.00023, 1.9e-5, 1.4e-6, 1.0e-7, 7.4e-9},
                       {0.00015, 7.2e-6, 3.8e-7, 2.1e-8, 1.3e-9, 7.8e-11},
                       {9.7e-6, 3.0e-7, 8.8e-9, 2.7e-10, 8.5e-12, 7.9e-13}  }; //allowed errors, pre-calculated and with added tolerance
  double convOrders[5]={1.6, 3.6, 4.15, 5.0}; //precalculated order of convergence, rounded down for tolerance
  
  std::stringstream message("Test succeeded", std::stringstream::out);
  
  
  for (int numOrder=0; numOrder<std::min(4,maxOrder); numOrder++)
  { 
    int maxRefs=std::min(5, maxRefinements);
    double postErrors[maxRefs+1]; //errors determinated while running test

  
  
    //create a rectancle grid   
    Vector x0(0);
    Vector edgeLengths(0); edgeLengths[0]=length; edgeLengths[1]=1.0;
    GridManager<Grid> gridManager( createRectangle<Grid, double >(x0, edgeLengths, 1) );
    
    for (int refs=0; refs<=maxRefs; refs++)
    {
      
      
      
      boost::timer::cpu_timer gridTimer;
      if (refs>0)
        gridManager.globalRefine(1);
      if (verbosity)
      {
        std::cout << "Grid: " << gridManager.grid().size(0) << " tetrahedra, " << std::endl;
        std::cout << "      " << gridManager.grid().size(dim-1) << " edges, " << std::endl;
        std::cout << "      " << gridManager.grid().size(dim) << " points" << std::endl;
        std::cout << "computing time for refinement of mesh: " << boost::timer::format(gridTimer.elapsed());
      } 
      // construction of finite element space for the scalar solution u
      H1Space h1Space(gridManager,gridManager.grid().leafGridView(),orders[numOrder]);
      
      auto varSetDesc = makeVariableSetDescription(makeSpaceList(&h1Space),
                                                   boost::fusion::make_vector(Variable<SpaceIndex<0>,Components<dim>>("u")));
      using VarSetDesc = decltype(varSetDesc);
    
      using Functional = ElasticityFunctional<VarSetDesc>;
      using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
      using CoefficientVectors = VarSetDesc::CoefficientVectorRepresentation<0,1>::type;
    
      //variational Functional
      Functional F(ElasticModulus(1,8), length);
    
      //construct Galerkin representation
      Assembler assembler(varSetDesc.spaces);
      VarSetDesc::VariableSet u(varSetDesc);
    
    
      boost::timer::cpu_timer assembTimer;
      assembler.assemble(linearization(F,u));
      if (verbosity)
        std::cout << "computing time for assemble: " << boost::timer::format(assembTimer.elapsed());
      
      AssembledGalerkinOperator<Assembler> A(assembler);
      CoefficientVectors solution(VarSetDesc::CoefficientVectorRepresentation<>::init(varSetDesc.spaces));
      CoefficientVectors rhs(assembler.rhs());
    
      boost::timer::cpu_timer directTimer;
      directInverseOperator(A).applyscaleadd(-1.0,rhs,solution);
      if (verbosity)
        std::cout << "computing time for direct solve: " << boost::timer::format(directTimer.elapsed());
      component<0>(u) = component<0>(solution);
      
      VarSetDesc::VariableSet func( varSetDesc ) ;
      
      interpolateGloballyWeak<PlainAverage>(boost::fusion::at_c<0>(func.data),ExactSolution());
      func-=u;
      L2Norm l2Norm;
      postErrors[refs] = l2Norm( boost::fusion::at_c<0>(func.data) ) ;
      if (verbosity)
        std::cout << "L2norm of error = " << postErrors[refs] << " for order " << orders[numOrder] << " and refinements " << refs << " , has to be smaller than " << preErrors[numOrder][refs] << std::endl;

      if (preErrors[numOrder][refs]<postErrors[refs])
      { 
          message << "Test failed: The error after " << refs << 
                 " refinements was too high at the test with ansatz functions of order " << orders[numOrder] << ".\n";
          valid=false;
      }
      
      
    }
    
     // calculate order of convergence by linear regression with \delta_i~log(e_i)-p*log(h_i)+log(c)
    double a11=0, a12=0, a22=5, b1=0, b2=0;
    if (maxRefs<=3 && verbosity)
       std::cout << "not enough points for useful convergence order\n";
    else
    {  
      for(int j=0; j<=4; j++)
      {
        
        postErrors[j] = std::log(postErrors[j]);
        a11 += hs[j]*hs[j];
        a12 += hs[j];
        b1 += hs[j]*postErrors[j];
        b2 += postErrors[j];
      }
      std::cout << a11 << " "<< a12 << " "<< a22 << " "<< b1 << " "<< b2 << " \n";
      double det = 1.0/(a11*a22-a12*a12);
      
      double convOrd = det*(a22*b1-a12*b2);
      double logc = det*(a11*b2-a12*b1);
      if (verbosity)
      {
        std::cout << "error = c*h^p with" << std::endl;
        std::cout << "p = " << convOrd << std::endl;
        std::cout << "log(c) = " << logc << std::endl;
      }
      if (convOrd< convOrders[numOrder])
        {
          message << "Test failed: The convergence order for order " << orders[numOrder] << 
                 " was too öw.\n";
          valid = false;
        }
    }
    
    
  }
  boost::timer::cpu_timer outputTimer;
  if (valid) std::cout <<"********************\nTest succeeded!!! \n********************\n";
  else if (!result) std::cout << "Test failed, try again with --result 1 to get a testResult.txt\n";
  else           std::cout << "Test failed, see testResult.txt\n";
  
  if(result)
  {
    std::string description = "Test with 2D bar elasticity problem \n";
    std::ofstream outfile("testResult.txt");
    outfile << description << std::endl << message.str() << std::endl << std::endl;
    outfile.close();
  }
  
  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed());
  std::cout << "End Bar program" << std::endl;
}
