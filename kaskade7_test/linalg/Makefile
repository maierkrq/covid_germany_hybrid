

LIBOBJ += $(patsubst %, linalg/%, domainDecompositionPreconditioner.o dynamicMatrix.o factorization.o umfpack_tools.o simpleLAPmatrix.o apcg.o threadedMatrix.o qp.o pqp.o qpmg.o)


LINATARGETOBJ = umfpack_solve.o mumps_solve.o

ifdef SUPERLUINC
  LINATARGETOBJ += superlu_solve.o
endif

ifdef PARDISOINC
  LINATARGETOBJ += pardiso_solve.o
endif

ifeq ($(USE_MPI),1)
  LINATARGETOBJ += mumpsmpi_solve.o
endif


# The target kasklib (complete Kaskade 7 library) depends not only on libkaskade.a, but also on
# the (optional) target objects.
kasklib: $(patsubst %, libs/%, $(LINATARGETOBJ))

# We need the libs directory created for copying target objects there
$(patsubst %, linalg/%, $(LINATARGETOBJ)):

# A static rule for copying target objects to libs directory. 
$(patsubst %, libs/%, $(LINATARGETOBJ)): libs/%: linalg/%
	mkdir -p libs
	cp $< $@
	
# Target objects do not end up in LIBOBJ, thus we mark the source files here
# for dependency generation.
DEPENDSRC += $(patsubst %.o, linalg/%.cpp, $(LINATARGETOBJ))

# 
# # unit test targets
linalg/ndtest: linalg/domainDecompositionPreconditioner.cpp kasklib
	$(CXX) -g -DUNITTEST $(INCLUDE) $(FLAGS) $< $(KASKADELIB) $(DUNELIB) $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB)  -o $@

linalg/tensor: linalg/tensor.cpp linalg/tensor.hh
	$(CXX) -DUNITTEST $(INCLUDE) $(FLAGS) $< $(KASKADELIB) $(DUNELIB)  $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@



# dynamicMatrix: dynamicMatrix.cpp
# 	$(CXX) -g -DUNITTEST $(INCLUDE) $(FLAGS) $<  $(DUNELIB) $(BOOSTLIB) $(BLASLIB) $(FTNLIB)  -o $@
# 
# numaMatrix: numaMatrix.cpp numaMatrix.hh 
# 	$(CXX) -g -DUNITTEST $(INCLUDE) $(FLAGS) $<  $(KASKADELIB) $(DUNELIB) $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB) -o $@
# 
# mumpsmpi_solve.o: mumpsmpi_solve.cpp mumpsmpi_solve.hh
# 	$(CXX) $(FLAGS) $(MUMPSMPIINC) $(INCLUDE) $(CLANGINC) $< -c -o $@ 2> gccerr.txt 
# 
# qp: qp.cpp qp.hh
# 	$(CXX) -g -DUNITTEST $(INCLUDE) $(FLAGS) $<  $(KASKADELIB) $(DUNELIB) $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB) -o $@
# 
# mumps: mumps_solve.cpp mumps_solve.hh
# 	$(CXX) -g -DUNITTEST $(INCLUDE) $(FLAGS) $<  $(KASKADELIB) $(MUMPSLIB) $(BLASLIB) $(FTNLIB) -o $@
# 
# pardiso: pardiso_solve.cpp pardiso_solve.hh
# 	$(CXX) -g -DUNITTEST $(INCLUDE) $(FLAGS) $<  $(KASKADELIB) $(BLASLIB) $(FTNLIB) -o $@
# 
