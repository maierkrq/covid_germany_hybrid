/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2019-2020 Zuse Institute Berlin                            */
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
#include "fem/lagrangespace.hh"
#include "linalg/direct.hh"
#include "linalg/icc0precond.hh"
#include "linalg/hyprecond.hh"       
#include "linalg/iluprecond.hh"
#include "linalg/cg.hh"
#include "io/vtk.hh"
#include "utilities/kaskopt.hh"
/**
  * \defgroup tests Tests
  * \brief Classes and functions that test whether certain parts of Kaskade work as they should.
  */
 /**
 * @file
 * @ingroup tests
 * @brief  Test with 3-dimensionale linear elasticity equation.
 *
 * Testprogram for Kaskade: tests wether the error between calculated and exact solution behaves as expected. 
 * This test was built by giving an exact solution u on
 * \f-4 \Omega:= \{ (x,y,z) \in R³ | 0<=x<=4, y²+z²<=1 \} \f$ (Cylinder, l=4, r=1}
 * and constructing a corresponding BVP.
 *  
 * The solution is given by 
 * 
 * \f[ u(x)=1/ \mu * (g_2(x_1)*e_2+g_3(x_1)*e_3), \f]
 * 
 * with
 *  
 * \f[ g_2(x)=sin(0.5* \pi*x) and g_3(x)=-4/(x+2)-x/3+2. \f]
 * 
 * We get the BVP:
 * 
 * \f[ - \bigtriangledown \cdot \sigma(u)(x) = f(x) = -g_2''(x_1)*e_2-g_3''(x_1)e_3 on \Omega
 *  \sigma(u)(x)* \eta = (g_2''(x_1)* \eta_2+g_3''(x_1)*\eta_3)*e_1 on \Theta_1 := \{ x \in \Omega | x_2²+x_3²=1 \}
 *  u= u0=0 on \Theta_2 := \{ x \in \Omega | x_1=0 oder x_1=4 \} \f]
 * 
 * Boundary conditions and f are adapted in elasto.hh. 
 * 
 * Uses a direct solver or cg Method to solve linear problem. However cg-Method produces greater errors for higher refinements/orders.
 */
using namespace Kaskade;
#include "elasto.hh"
struct ExactSolution //exact solution of the  BVP
{
   using Scalar = double;
   static int const components = 3;
   using ValueType = Dune::FieldVector<Scalar,components>;

   template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }

   template <class Cell>
   ValueType value(Cell const& cell,Dune::FieldVector<typename Cell::Geometry::ctype,Cell::dimension> const& localCoordinate) const
   {
     Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x = cell.geometry().global(localCoordinate);
    double pi=3.141592653589793238462;
     ValueType y; y[0]=0; y[1]=0.125*sin(0.5*pi*x[0]); y[2]=-0.125*(4.0/(x[0]+2)+1.0/3.0*x[0]-2);
     return y;
   }
};

//transforms cuboid of lenght l to cylinder of length l, radius
class CuboidToCylinder: public Dune::AnalyticalCoordFunction< double,3,3,CuboidToCylinder>
{
  using Base = Dune::AnalyticalCoordFunction<double,3,3,CuboidToCylinder>;

  public:
  using DomainVector = Base::DomainVector;
  using RangeVector = Base::RangeVector;
  
  CuboidToCylinder (double _radius=1.0)
  : radius(_radius)
  {}
  
  void evaluate(const DomainVector &u, RangeVector &y) const
  {
    double enorm = sqrt(u[2]*u[2]+u[1]*u[1]);
    if(enorm > 1e-5)
    {
      double radiusFactor = radius*std::max(std::fabs(u[2]),std::fabs(u[1]));
      double scaling = radiusFactor/sqrt(u[2]*u[2]+u[1]*u[1]);
      y[ 2 ] = u[ 2 ]*scaling;
      y[ 1 ] = u[ 1 ]*scaling;
      y[ 0 ] = u[ 0 ];
    } 
    else
    {
      y[ 2 ] = u[ 2 ];
      y[ 1 ] = u[ 1 ];
      y[ 0 ] = u[ 0 ];
    }
  }
  
  private:
    double radius;
};



int main(int argc, char *argv[])
{
  std::cout << "Start Cylinder program\n" ;

  boost::timer::cpu_timer totalTimer;

  int maxRefinements, maxOrder, writeResult, longer, iterative, verbosity;
  double radius,length;
  bool valid=true;

  if (getKaskadeOptions(argc,argv,Options
    ("refinements",      maxRefinements,    2,      "maximum number of uniform grid refinements")
    ("order",            maxOrder,          2,      "maximum finite element ansatz order")
    ("result",           writeResult,            0,      "set 1 if testResult.text should be written")
    ("longer",           longer,            0,      "set 1 for a longer test, set 2 for an even longer test (not practical)")
    ("iterative",        iterative,         0,      "set 1 for using a cg-method")
    ("verbosity",        verbosity,         1,      "reporting level, 0 or 1")
    ("radius",           radius,  		    	1.0,    "radius of the cylinder") //shouldn't be changed for test to work
    ("length",			     length,		      	4.0,    "length of the cylinder")  //shouldn't be changed for test to work
    ))
    return 1;


  

 // THREE-dimensional space: dim=3
  constexpr int dim=3;        
  using Grid = Dune::UGGrid<dim>;
  using GeoGrid = Dune::GeometryGrid<Grid,CuboidToCylinder>;
  using LeafView = GeoGrid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  
  if (longer)
  {
    maxRefinements=4;
    maxOrder=4;
  }
  
  int orders[4]={1,2,3,4}; //max refinements possible for corresponding order so that direct solver can still solve the problem
  int refinements[4];
  if (longer==2 && iterative)
  {
    refinements[0]=4; refinements[1]=4; refinements[2]=4; refinements[3]=4;
  }
  else
  { refinements[0]=4; refinements[1]=3; refinements[2]=2; refinements[3]=2;}

  // logs of meshsizes for corresponding refinements (0 to 4)
  double hs[5]={2.0,sqrt(2),sqrt(0.5),sqrt(0.125), sqrt(0.03125)}; 
  for (int m=0; m<5; m++)
    hs[m] = std::log(hs[m]);
  
  double preErrors[4][5]={ {0.11, 0.064, 0.025, 0.0071, 0.0018},
                           {0.02, 0.003, 0.00027, 2.2e-5, 2.0e-6},
                           {0.007, 0.0003, 1.6e-5, 8.5e-7, 7.9e-8},
                           {0.00014, 1.8e-5, 5.6e-7, 2.8e-8, -1}  }; //allowed errors, pre-calculated and rounded up for tolerance, depend on order and refinements
  double convOrders[5]={1.4, 3.6, 5, 5}; //convorder depending on order of ansatz func
  
  std::stringstream message("Test succeeded", std::stringstream::out);
  
  
  for (int numOrder=0; numOrder<std::min(4,maxOrder); numOrder++)
  { 
    int maxRefs=std::min(refinements[numOrder], maxRefinements);
    double postErrors[maxRefs+1]; //errors determinated while running test
  
    // Creating a Rectancle out of cubes with edgelength 2
    
    Dune::GridFactory<Grid> factory;
    
    
    Dune::FieldVector<double,dim> v;  // vertex coordinates v[0], v[1], v[2]
    double inb=1.0;
    int cubeNumber=std::max(floor(length/2+0.5), inb);
    v[0]=0; v[1]=-1; v[2]=-1; factory.insertVertex(v);
    v[0]=0; v[1]=-1; v[2]=1; factory.insertVertex(v);
    v[0]=0; v[1]=1; v[2]=1; factory.insertVertex(v);
    v[0]=0; v[1]=1; v[2]=-1; factory.insertVertex(v);
    int k=3; //number of points in grid -1
    
    
    std::vector<unsigned int> vid(4);    // pyramid defined by 4 vertex indices
    Dune::GeometryType gt(Dune::GeometryType::simplex,3);
    for (int i=1; i<=cubeNumber; i++)
    {  
      v[0]=length/cubeNumber*i-length/cubeNumber/2; v[1]=0; v[2]=0; factory.insertVertex(v);
      v[0]=length/cubeNumber*i; v[1]=-1; v[2]=-1; factory.insertVertex(v);
      v[0]=length/cubeNumber*i; v[1]=-1; v[2]=1; factory.insertVertex(v);
      v[0]=length/cubeNumber*i; v[1]=1; v[2]=1; factory.insertVertex(v);
      v[0]=length/cubeNumber*i; v[1]=1; v[2]=-1; factory.insertVertex(v);
      k=k+5;
      
      vid[0]=k-8; vid[1]=k-7; vid[2]=k-5; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-6; vid[1]=k-7; vid[2]=k-5; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-8; vid[1]=k-3; vid[2]=k-5; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k  ; vid[1]=k-3; vid[2]=k-5; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-8; vid[1]=k-7; vid[2]=k-3; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-2; vid[1]=k-7; vid[2]=k-3; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-5; vid[1]=k  ; vid[2]=k-6; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-1; vid[1]=k  ; vid[2]=k-6; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-7; vid[1]=k-2; vid[2]=k-6; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-1; vid[1]=k-2; vid[2]=k-6; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-3; vid[1]=k-2; vid[2]=k  ; vid[3]=k-4; factory.insertElement(gt,vid);
      vid[0]=k-1; vid[1]=k-2; vid[2]=k  ; vid[3]=k-4; factory.insertElement(gt,vid);
    }
  
    Grid* grid = factory.createGrid();
    GridManager<GeoGrid> gridManager(new GeoGrid(grid,new CuboidToCylinder(radius)));    
    
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
      using Operator=AssembledGalerkinOperator<Assembler>;
      //variational Functional
      Functional F(ElasticModulus(1,8), length, radius);
    
      //construct Galerkin representation
      Assembler assembler(varSetDesc.spaces);
      VarSetDesc::VariableSet u(varSetDesc);
    std::cout << "assembling\n";
    
      boost::timer::cpu_timer assembTimer;
      assembler.assemble(linearization(F,u));
      if (verbosity)
        std::cout << "computing time for assemble: " << boost::timer::format(assembTimer.elapsed());
      
      Operator A(assembler);
      CoefficientVectors solution(VarSetDesc::CoefficientVectorRepresentation<>::init(varSetDesc.spaces));
      CoefficientVectors rhs(assembler.rhs());
      
      
      if (iterative) //solving linear problem using cg method with iluk preconditioner.    
      {
        Dune::InverseOperatorResult res;
        int iteSteps = 5000;
        double iteEps = 1.0e-8;
        StrakosTichyPTerminationCriterion<double> termination(iteEps,iteSteps);
        int lookAhead=50;
        termination.setLookAhead(lookAhead);
        
        const DefaultDualPairing<CoefficientVectors,CoefficientVectors> defaultScalarProduct{};
            int fill_lev =  1;
            ILUKPreconditioner<AssembledGalerkinOperator<Assembler,0,1,0,1> > iluk(A,fill_lev,verbosity);
            Dune::BiCGSTABSolver<CoefficientVectors> cg(A,iluk,iteEps,iteSteps,verbosity);
            cg.apply(solution,rhs,res);
        solution *= -1.0;
        u.data = solution.data;
      }
      
      else 
      {
        directInverseOperator(A).applyscaleadd(-1.0,rhs,solution);
        component<0>(u) = component<0>(solution);  
      }
      
      // exact solution in Ansatzspace
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
    double a11=0, a12=0, a22=maxRefs+1, b1=0, b2=0;
 
    for(int j=0; j<=maxRefs; j++)
    {

      postErrors[j] = std::log(postErrors[j]);
      a11 += hs[j]*hs[j];
      a12 += hs[j];
      b1 += hs[j]*postErrors[j];
      b2 += postErrors[j];
    }

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
      message << "Test failed: The convergence order for ansatz functions of order " << orders[numOrder] <<
      " was too low.\n";
      valid = false;
    }
  }

  if (valid)
    std::cout << "********************\nTest succeeded!!! \n********************\n";
  else if (writeResult)
    std::cout << "Test failed, see testResult.txt\n";
  else
    std::cout << "Test failed, try again with --result 1 to get a testResult.txt\n";

  if(writeResult)
  {
    std::string description = "Test with simple 3D elasticity problem and cg-solver\n";
    std::ofstream outfile("testResult.txt");
    outfile << description << std::endl << message.str() << std::endl << std::endl;
    outfile.close();
  }
  
  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed());
  std::cout << "End Cylinder program" << std::endl;
}
