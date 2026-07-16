include ../../Makefile.Local
include ../../Makefile.Rules

INCLUDE = $(DUNEINC) $(UGINC) $(BOOSTINC) $(AMIRAINC) $(KASKADEINC) $(ITSOLINC) $(HYPREINC) $(TAUCSINC) $(TRILINOSINC) $(UMFPACKINC) $(MUMPSINC) $(PSURFACEINC) -I.
FLAGS = $(WFLAGS) $(OPTFLAGS)
PLATFORM= $(shell uname)


amiraTest: amiraTest.o Makefile $(KASKADE7)/libs/libkaskade.a 
	$(CXX) $(FLAGS) $< $(KASKADE7)/libs/umfpack_solve.o $(KASKADE7)/libs/mumps_solve.o \
	$(KASKADE7)/libs/superlu_solve.o $(KASKADELIB) \
 $(DUNELIB) \
 $(UGLIB) \
 $(BOOSTLIB) \
 $(UMFPACKLIB) \
 $(MUMPSLIB) \
 $(ITSOLLIB) \
 $(HYPRELIB) \
 $(TAUCSLIB) \
 $(TRILINOSLIB) \
 $(PARDISOLIB) \
 $(AMIRALIB) \
 $(BLASLIB) $(FTNLIB) $(NUMALIB) -o $@
 
 clean:
	rm -f gccerr.txt tmpFile.txt tmp.txt amiraTest amiraTest.o checkDefault* randomCube_out.am .goutputstream-* 
 
 amiraTest.o: amiraTest.cpp \
 $(KASKADE7)/utilities/enums.hh \
 $(KASKADE7)/fem/assemble.hh \
 $(KASKADE7)/fem/firstless.hh \
 $(KASKADE7)/fem/fixfusion.hh \
 $(KASKADE7)/fem/functional_aux.hh \
 $(KASKADE7)/fem/functionspace.hh \
 $(KASKADE7)/fem/gridmanager.hh \
 $(KASKADE7)/utilities/threading.hh \
 $(KASKADE7)/utilities/kalloc.hh \
 $(KASKADE7)/utilities/timing.hh \
 $(KASKADE7)/fem/pshapefunctions.hh \
 $(KASKADE7)/fem/barycentric.hh \
 $(KASKADE7)/fem/fixdune.hh \
 $(KASKADE7)/fem/mllgeometry.hh \
 $(KASKADE7)/linalg/dynamicMatrix.hh \
 $(KASKADE7)/utilities/scalar.hh \
 $(KASKADE7)/fem/shapefunctioncache.hh \
 $(KASKADE7)/fem/fetransfer.hh \
 $(KASKADE7)/linalg/crsutil.hh \
 $(KASKADE7)/utilities/detailed_exception.hh \
 $(KASKADE7)/fem/quadrature.hh \
 $(KASKADE7)/fem/variables.hh \
 $(KASKADE7)/fem/linearspace.hh \
 $(KASKADE7)/io/iobase.hh \
 $(KASKADE7)/linalg/localMatrices.hh \
 $(KASKADE7)/linalg/threadedMatrix.hh \
 $(KASKADE7)/utilities/duneInterface.hh \
 $(KASKADE7)/utilities/is_function.hh \
 $(KASKADE7)/utilities/member_variable_macro.hh \
 $(KASKADE7)/utilities/typeTraits.hh \
 $(KASKADE7)/fem/spaces.hh \
 $(KASKADE7)/fem/hierarchicspace.hh \
 $(KASKADE7)/fem/converter.hh \
 $(KASKADE7)/fem/gridcombinatorics.hh \
 $(KASKADE7)/fem/hierarchicshapefunctions.hh \
 $(KASKADE7)/fem/combinatorics.hh \
 $(KASKADE7)/linalg/simpleLAPmatrix.hh \
 $(KASKADE7)/linalg/linearsystem.hh \
 $(KASKADE7)/linalg/triplet.hh \
 $(KASKADE7)/linalg/umfpack_tools.hh \
 $(KASKADE7)/fem/scalarspace.hh \
 $(KASKADE7)/fem/views.hh \
 $(KASKADE7)/fem/sign.hh \
 $(KASKADE7)/utilities/power.hh \
 $(KASKADE7)/fem/lagrangespace.hh \
 $(KASKADE7)/fem/lagrangeshapefunctions.hh \
 $(KASKADE7)/utilities/deprecation.hh \
 $(KASKADE7)/io/vtk.hh \
 $(KASKADE7)/io/check_endianness.hh \
 $(KASKADE7)/io/vtkwriter.hh \
 $(KASKADE7)/io/matlab.hh \
 $(KASKADE7)/linalg/direct.hh \
 $(KASKADE7)/fem/istlinterface.hh \
 $(KASKADE7)/linalg/factorization.hh \
 $(KASKADE7)/linalg/factory.hh \
 $(KASKADE7)/linalg/symmetricOperators.hh \
 $(KASKADE7)/linalg/apcg.hh \
 $(KASKADE7)/utilities/geometric_sequence.hh \
 $(KASKADE7)/mg/multigrid.hh \
 $(KASKADE7)/mg/additiveMultigrid.hh \
 $(KASKADE7)/linalg/conjugation.hh \
 $(KASKADE7)/linalg/domainDecompositionPreconditioner.hh \
 $(KASKADE7)/fem/supports.hh \
 $(KASKADE7)/fem/scalarspace.hh \
 $(KASKADE7)/linalg/jacobiPreconditioner.hh \
 $(KASKADE7)/linalg/matrixTraits.hh \
 $(KASKADE7)/mg/prolongation.hh \
 $(KASKADE7)/fem/deforminggridmanager.hh \
 $(KASKADE7)/fem/boundaryInterpolation.hh \
 $(KASKADE7)/fem/forEach.hh \
 $(KASKADE7)/mg/multiplicativeMultigrid.hh \
 $(KASKADE7)/utilities/memory.hh \
 $(KASKADE7)/utilities/gridGeneration.hh \
 $(KASKADE7)/utilities/GridGeneration/cube.hh \
 $(KASKADE7)/utilities/vectorTools.hh \
 $(KASKADE7)/utilities/kaskopt.hh \
 $(KASKADE7)/fem/diffops/elastoVariationalFunctionals.hh \
 $(KASKADE7)/fem/diffops/materialLaws.hh \
 $(KASKADE7)/fem/diffops/elasto.hh \
 $(KASKADE7)/fem/integration.hh \
 $(KASKADE7)/linalg/determinant.hh \
 $(KASKADE7)/linalg/adjugate.hh \
 $(KASKADE7)/linalg/cofactor.hh \
 $(KASKADE7)/linalg/wrappedMatrix.hh

