/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/en/numerik/software/kaskade-7.html               */
/*                                                                           */
/*  Copyright (C) 2002-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/norms.hh"
#include "fem/lagrangespace.hh"
//#include "fem/hierarchicspace.hh"   // ContinuousHierarchicMapper
#include "linalg/direct.hh"
#include "linalg/partialDirectPreconditioner.hh"
#include "linalg/iluprecond.hh"      // ILUT, ILUK, ARMS
#include "linalg/icc0precond.hh"
#include "linalg/hyprecond.hh"       // BoomerAMG, Euclid
#include "linalg/jacobiPreconditioner.hh"
#include "linalg/cg.hh"
#include "mg/hb.hh"
#include "utilities/enums.hh"
#include "utilities/gridGeneration.hh" //  createUnitSquare, createUnitCube
#include "io/vtk.hh"
//#include "io/amira.hh"
#include "utilities/kaskopt.hh"

#include "cubus.hh"

using namespace Kaskade;
#include "ht.hh"

#define DIM 3

using namespace boost::fusion;

boost::timer::cpu_times operator / (boost::timer::cpu_times times, double d)
{
  return {(long int) (times.wall/d), (long int) (times.user/d), (long int) (times.system/d)};
}

int main(int argc, char *argv[])
{
  using std::min;

  std::cout << "Start heat transfer tutorial program" << std::endl;

  boost::timer::cpu_timer totalTimer;

  int verbosityOpt = 1;
  bool dump = false; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  int  refinements = getParameter(pt, "refinements", 5),
       order       = getParameter(pt, "order", 2),
       verbosity   = getParameter(pt, "verbosity", 1),
       init        = getParameter(pt, "init", 0);
  std::cout << "original mesh shall be refined : " << refinements << " times" << std::endl;
  std::cout << "discretization order           : " << order << std::endl;
  std::cout << "output level (verbosity)       : " << verbosity << std::endl;

  int nthreads = getParameter(pt,"threads",4);
  NumaThreadPool::instance(nthreads);

  int  direct;
  bool onlyLowerTriangle = false;
    
  DirectType directType;
//   IterateType iterateType = IterateType::CG;
  MatrixProperties property = MatrixProperties::SYMMETRIC;
  PrecondType precondType = PrecondType::NONE;
  std::string empty;

  std::string s("names.type.");
  s += getParameter(pt, "solver.type", empty);
  direct = getParameter(pt, s, 0);
    
  s = "names.direct." + getParameter(pt, "solver.direct", empty);
  directType = static_cast<DirectType>(getParameter(pt, s, 0));

//   s = "names.iterate." + getParameter(pt, "solver.iterate", empty);
//   iterateType = static_cast<IterateType>(getParameter(pt, s, 0));
  s = "names.preconditioner." + getParameter(pt, "solver.preconditioner", empty);
  precondType = static_cast<PrecondType>(getParameter(pt, s, 0));
  
  bool output = getParameter(pt,"output",false),
       solve  = getParameter(pt,"solve",true);
  int blocks = getParameter(pt,"blocks",40);
  double rowBlockFactor = getParameter(pt,"rowBlockFactor",2.0);

  property = MatrixProperties::SYMMETRIC;

  if ( (directType == DirectType::MUMPS)||(directType == DirectType::PARDISO) || ( (precondType == PrecondType::ICC) && !direct ) )
  {
    onlyLowerTriangle = true;
    std::cout << 
      "Note: direct solver MUMPS/PARADISO or ICC preconditioner ===> onlyLowerTriangle is set to true!" 
      << std::endl;
  }

  boost::timer::cpu_timer gridTimer;
#if DIM==2
  //   two-dimensional space: dim=2
  int const dim=2;        
  typedef Dune::UGGrid<dim> Grid;
  // There are alternatives to UGGrid: ALUSimplexGrid (red refinement), ALUConformGrid (bisection)
  //typedef Dune::ALUSimplexGrid<2,2> Grid;
  //typedef Dune::ALUConformGrid<2,2> Grid;
  Grid::setDefaultHeapSize(2048);
  Dune::GridFactory<Grid> factory;

  // vertex coordinates v[0], v[1]
  Dune::FieldVector<double,dim> v;    
  v[0]=0; v[1]=0; factory.insertVertex(v);
  v[0]=1; v[1]=0; factory.insertVertex(v);
  v[0]=1; v[1]=1; factory.insertVertex(v);
  v[0]=0; v[1]=1; factory.insertVertex(v);
  // triangle defined by 3 vertex indices
  std::vector<unsigned int> vid(3);
  Dune::GeometryType gt(Dune::GeometryType::simplex,2);
  vid[0]=0; vid[1]=1; vid[2]=2; factory.insertElement(gt,vid);
  vid[0]=0; vid[1]=2; vid[2]=3; factory.insertElement(gt,vid);
  std::unique_ptr<Grid> grid( factory.createGrid() ) ;
  // the coarse grid will be refined three times
  grid->globalRefine(refinements);
  // some information on the refined mesh
  std::cout << std::endl << "Grid: " << grid->size(0) << " triangles, " << std::endl;
  std::cout << "      " << grid->size(1) << " edges, " << std::endl;
  std::cout << "      " << grid->size(2) << " points" << std::endl;

  // a gridmanager is constructed 
  // as connector between geometric and algebraic information
  GridManager<Grid> gridManager(std::move(grid));
#else
  //  three-dimensional space: dim=3
  //  we offer 2 geometries:
  //  - very simple:  1 tetrahedron defined by 4 vertices
  //  - more complex: 1 cube defined by 48 tetrahedra, provided in cubus.hh
  int const dim=3; 
    
  /*    
  //definition of 1 tetrahedron by 4 vertices
  typedef Dune::UGGrid<dim> Grid;
  Dune::GridFactory<Grid> factory;
  // vertex coordinates v[0], v[1]
  Dune::FieldVector<double,dim> v;    
  v[0]=0; v[1]=0; v[2]=0; factory.insertVertex(v);
  v[0]=1; v[1]=0; v[2]=0; factory.insertVertex(v);
  v[0]=0; v[1]=1; v[2]=0; factory.insertVertex(v);
  v[0]=0; v[1]=0; v[2]=1; factory.insertVertex(v);
  // tetrahedron defined by 4 vertex indices
  std::vector<unsigned int> vid(4);
  Dune::GeometryType gt(Dune::GeometryType::simplex,dim);
  vid[0]=0; vid[1]=1; vid[2]=2; vid[3]=3; factory.insertElement(gt,vid);
  std::unique_ptr<Grid> grid( factory.createGrid() ) ;
  // the coarse grid will be refined three times
  grid->globalRefine(refinements);
  */
    
  // definition of a more complex mesh using cubus.hh
  int heapSize = getParameter(pt,"heapsize",2048);
  typedef Dune::UGGrid<dim> Grid;
  std::unique_ptr<Grid> grid( RefineGrid<Grid>(refinements, heapSize) );

  // some information on the refined mesh
  std::cout << std::endl << "Grid: " << grid->size(0) << " tetrahedra, " << std::endl;
  std::cout << "      " << grid->size(1) << " triangles, " << std::endl;
  std::cout << "      " << grid->size(dim-1) << " edges, " << std::endl;
  std::cout << "      " << grid->size(dim) << " points" << std::endl;
  // a gridmanager is constructed 
  // as connector between geometric and algebraic information
  GridManager<Grid> gridManager(std::move(grid));
#endif
  std::cout << "computing time for generation of initial mesh: " << boost::timer::format(gridTimer.elapsed()) << "\n";
    
  using LeafView = Grid::LeafGridView;
  using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> >;
  // using H1Space = FEFunctionSpace<ContinuousHierarchicMapper<double,LeafView> >;
  using Spaces = boost::fusion::vector<H1Space const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>,Components<1>>>;
  using VariableSet = VariableSetDescription<Spaces,VariableDescriptions>;
  using Functional = HeatFunctional<double,VariableSet>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional> >;
  constexpr int neq = Functional::TestVars::noOfVariables;
  constexpr int nvars = Functional::AnsatzVars::noOfVariables;
  using CoefficientVectors = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
  using LinearSpace = VariableSet::CoefficientVectorRepresentation<0,neq>::type;
    
  // construction of finite element space for the scalar solution T.
  H1Space temperatureSpace(gridManager,gridManager.grid().leafGridView(),order);
    
  Spaces spaces(&temperatureSpace);
    
  // construct variable list.
  // VariableDescription<int spaceId, int components, int Id>
  // spaceId: number of associated FEFunctionSpace
  // components: number of components in this variable
  // Id: number of this variable
        
  std::string varNames[1] = { "u" };
    
  VariableSet variableSet(spaces,varNames);

  // construct variational functional
    
  double kappa = 1.0;
  double q = 1.0;
  Functional F(kappa,q);
  std::cout << std::endl << "no of variables = " << nvars << std::endl;
  std::cout << "no of equations = " << neq   << std::endl;
  size_t dofs = variableSet.degreesOfFreedom(0,nvars);
  std::cout << "number of degrees of freedom = " << dofs   << std::endl;

    
  Assembler assembler(gridManager,spaces);
  VariableSet::VariableSet u(variableSet);
  VariableSet::VariableSet du(variableSet);

  size_t nnz = assembler.nnz(0,1,0,1,onlyLowerTriangle);
  std::cout << "number of nonzero elements in the stiffness matrix: " << nnz << std::endl << std::endl;
  
  CoefficientVectors solution(VariableSet::CoefficientVectorRepresentation<0,1>::init(spaces));
  solution = 0;
  
  std::cout << "Maximum number of hardware threads: " << boost::thread::hardware_concurrency() << std::endl;
  //gridManager.setEnforceConcurrentReads(true);
  gridManager.enforceConcurrentReads(true);
  assembler.setNSimultaneousBlocks(blocks);
  assembler.setRowBlockFactor(rowBlockFactor);
  
  boost::timer::cpu_times assembTimes;
  if ( init == 0 ) 
  {
    assembler.assemble(linearization(F,u),assembler.MATRIX|assembler.RHS|assembler.VALUE,nthreads,verbosity);
    boost::timer::cpu_timer assembTimer;
    assembler.assemble(linearization(F,u),assembler.MATRIX|assembler.RHS|assembler.VALUE,nthreads,verbosity);
    assembTimes = assembTimer.elapsed();
    std::cout << "computing time for assemble: " << boost::timer::format(assembTimes);
  }
  if ( init == 1 ) 
  {
    boost::timer::cpu_timer assembTimer;
    for (int i=1;i<=5;i++)
      assembler.assemble(linearization(F,u),assembler.MATRIX|assembler.RHS|assembler.VALUE,nthreads,verbosity);
    assembTimes = assembTimer.elapsed()/5;
    std::cout << "computing time for assemble: " << boost::timer::format(assembTimes);
    std::ofstream initialout;
    if ( nthreads == 1 ) initialout.open("oneThread3d.time",std::ofstream::out);
    if ( nthreads == 2 ) initialout.open("twoThread3d.time",std::ofstream::out);
    initialout << assembTimes.wall << std::endl;
    initialout.close();
    return 0; 
  }
  
  std::ofstream sumout;
  sumout.open("assembsummary3d.stat",std::ofstream::app);

  double oneThreadTime, twoThreadTime;
  std::ifstream initialin;
  initialin.open("oneThread3d.time",std::ifstream::in);
  initialin >> oneThreadTime;
  initialin.close();
  initialin.open("twoThread3d.time",std::ifstream::in);
  initialin >> twoThreadTime;
  initialin.close();

  sumout << std::setw(3) << nthreads << "  " << std::setw(15) << std::setprecision(5);
  
  if ( nthreads == 1 )
    sumout << std::fixed << "1.00000" << "  " << 2.0*twoThreadTime/(double)(assembTimes.wall);
  
  if ( nthreads == 2 )
    sumout << std::fixed << oneThreadTime/(double)(assembTimes.wall) << "  " << "2.00000";

  if ( nthreads > 2 )
    sumout << std::fixed << oneThreadTime/(double)(assembTimes.wall) << "  " << 2.0*twoThreadTime/(double)(assembTimes.wall);
  
  sumout << "  " << std::setw(11) << std::setprecision(0) << std::fixed << assembTimes.wall << std::endl;

  sumout.close();
  
  if ( solve ) {
    CoefficientVectors rhs(assembler.rhs());
    AssembledGalerkinOperator<Assembler,0,1,0,1> A(assembler, onlyLowerTriangle);
    
    // matrix may be used in triplet format, e.g.,
    // MatrixAsTriplet<double> tri = A.get<MatrixAsTriplet<double> >();
    //     for (k=0; k< nnz; k++)
    //       {
    //         printf("%3d %3d %e\n", tri.ridx[k], tri.cidx[k], tri.data[k]);
    //       }
  
    if (direct)
    {
      boost::timer::cpu_timer directTimer;
      directInverseOperator(A,directType,property).applyscaleadd(-1.0,rhs,solution);
      u.data = solution.data;
      std::cout << "computing time for direct solve: " << boost::timer::format(directTimer.elapsed());
    }
    else
   {
    boost::timer::cpu_timer iteTimer;
    Dune::InverseOperatorResult res;
    const DefaultDualPairing<LinearSpace,LinearSpace> defaultScalarProduct{};
    int iteSteps = getParameter(pt, "solver.iteMax", 2000);
    double iteEps = getParameter(pt, "solver.iteEps", 1.0e-10);
    StrakosTichyPTerminationCriterion<double> termination(iteEps,iteSteps);
    int lookAhead;
    switch (precondType)
    {
      case PrecondType::NONE:
      case PrecondType::HB:   lookAhead=50; break;
      default:                lookAhead=3; break;
    }
    lookAhead = getParameter(pt, "solver.lookAhead", lookAhead);
    termination.setLookAhead(lookAhead);

    switch (precondType)
    {
      case PrecondType::NONE:
      {
        std::cout << "selected preconditioner: NONE" << std::endl;
        IdentityPreconditioner<typename AssembledGalerkinOperator<Assembler,0,neq,0,nvars>::domain_type> trivial;
        CG<LinearSpace,LinearSpace> cg(A,trivial,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ILUT:
      {
        std::cout << "selected preconditioner: ILUT" << std::endl;
        std::cout << "needs matrix.property = GENERAL" << std::endl;
        int lfil = getParameter(pt, "solver.ILUT.lfil", 140);
        double dropTol = getParameter(pt, "solver.ILUT.dropTol", 0.01);
        ILUTPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > ilut(A,lfil,dropTol,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,ilut,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ILUK:
      {
        std::cout << "selected preconditioner: ILUK" << std::endl;
//        std::cout << "Note that this preconditioner combined with the BICGSTAB solver" << std::endl;
        std::cout << "needs matrix.property = GENERAL" << std::endl;
        int fill_lev = getParameter(pt, "solver.ILUK.fill_lev", 3);
        ILUKPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > iluk(A,fill_lev,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,iluk,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::BiCGSTABSolver<LinearSpace> cg(A,iluk,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ARMS:
      {
        int lfil = getParameter(pt, "solver.ARMS.lfil", 140);
        int lev_reord = getParameter(pt, "solver.ARMS.lev_reord", 1);
        double dropTol = getParameter(pt, "solver.ARMS.dropTol", 0.01);
        double tolind = getParameter(pt, "solver.ARMS.tolind", 0.2);
        ARMSPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > iluk(A,lfil,dropTol,lev_reord,tolind,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,iluk,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::ICC0:
      {
        std::cout << "selected preconditioner: ICC0" << std::endl;
        ICC_0Preconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > icc0(A);
        CG<LinearSpace,LinearSpace> cg(A,icc0,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::HB:
      {
        std::cout << "selected preconditioner: HB" << std::endl;
        HierarchicalBasisPreconditioner<Grid,AssembledGalerkinOperator<Assembler,0,neq,0,nvars>::range_type, AssembledGalerkinOperator<Assembler,0,neq,0,nvars>::range_type > hb(gridManager.grid());
        CG<LinearSpace,LinearSpace> cg(A,hb,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::BOOMERAMG:
      {
        int steps = getParameter(pt, "solver.BOOMERAMG.steps", iteSteps);
        int coarsentype = getParameter(pt, "solver.BOOMERAMG.coarsentype", 21);
        int interpoltype = getParameter(pt, "solver.BOOMERAMG.interpoltype", 0);
        int cycleType = getParameter(pt, "solver.BOOMERAMG.cycleType", 1);
        int relaxType = getParameter(pt, "solver.BOOMERAMG.relaxType", 3);
        int variant = getParameter(pt, "solver.BOOMERAMG.variant", 0);
        int overlap = getParameter(pt, "solver.BOOMERAMG.overlap", 1);
        double tol = getParameter(pt, "solver.BOOMERAMG.tol", iteEps);
        double strongThreshold = getParameter(pt, "solver.BOOMERAMG.strongThreshold", (dim==2)?0.25:0.6);
        BoomerAMG<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> >
                  boomerAMGPrecon(A,steps,coarsentype,interpoltype,tol,cycleType,relaxType,
                  strongThreshold,variant,overlap,1,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,boomerAMGPrecon,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
//         Dune::LoopSolver<LinearSpace> cg(A,boomerAMGPrecon,iteEps,iteSteps,verbosity);
//         cg.apply(solution,rhs,res);
      }
      break;
      case PrecondType::EUCLID:
      {
        std::cout << "selected preconditioner: EUCLID" << std::endl;
        int level      = getParameter(pt, "solver.EUCLID.level",1);
        double droptol = getParameter(pt, "solver.EUCLID.droptol",0.01);
        int printlevel = 0;
        if (verbosity>2) printlevel=verbosity-2;
        printlevel = getParameter(pt,"solver.EUCLID.printlevel",printlevel);
        int bj = getParameter(pt, "solver.EUCLID.bj",0);
        Euclid<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > EuclidPrecon(A,level,droptol,printlevel,bj,verbosity);
        CG<LinearSpace,LinearSpace> cg(A,EuclidPrecon,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
     }
      break;
      case PrecondType::JACOBI:
      default:
      {
        std::cout << "selected preconditioner: JACOBI" << std::endl;
        JacobiPreconditioner<AssembledGalerkinOperator<Assembler,0,neq,0,nvars> > jacobi(A,1.0);
        CG<LinearSpace,LinearSpace> cg(A,jacobi,defaultScalarProduct,termination,verbosity);
        cg.apply(solution,rhs,res);
      }
      break;
    }
    solution *= -1.0;
    u.data = solution.data;
    
    std::cout << "iterative solve eps= " << iteEps << ": " 
              << (res.converged?"converged":"failed") << " after "
              << res.iterations << " steps, rate="
              << res.conv_rate << ", computing time=" << (double)(iteTimer.elapsed().user)/1e9 << "s\n";
  }
   
      
    
    // compute L2 norm of the solution
    boost::timer::cpu_timer outputTimer;
    L2Norm l2Norm;
    std::cout << "L2norm(solution) = " << l2Norm(boost::fusion::at_c<0>(u.data)) << std::endl;
    

    // output of solution in VTK format for visualization,
    // the data are written as ascii stream into file temperature.vtu,
    // possible is also binary
    if ( output ) {
      IoOptions options;
      options.outputType = IoOptions::ascii;
      LeafView leafView = gridManager.grid().leafGridView();
      writeVTKFile(leafView,u,"temperature",options,min(order,2));
      std::cout << "graphical output finished, data in VTK format is written into file temperature.vtu \n";
        
      // output of solution for Amira visualization,
      // the data are written in binary format into file temperature.am,
      // possible is also ascii
      // IoOptions options;
      // options.outputType = IoOptions::ascii;
      // LeafView leafView = gridManager.grid().leafView();
      // writeAMIRAFile(leafView,variableSet,u,"temperature",options);
      std::cout << "computing time for output: " << boost::timer::format(outputTimer.elapsed()) << "\n";
    }
  }


  std::cout << "total computing time: " << boost::timer::format(totalTimer.elapsed()) << "\n";
  std::cout << "End heat transfer tutorial program" << std::endl;
}
