
LIBOBJ += $(patsubst %, mg/%, prolongation.o add.o)


# # unit test targets
mg/prolongationtest: mg/prolongation.cpp
	$(CXX) -g -DUNITTEST $(INCLUDE) $(FLAGS) $< $(KASKADELIB) $(DUNELIB) $(UGLIB) $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB)  -o $@
# 
# multigrid.o: multigrid.cpp multigrid.hh prolongation.hh
# 	$(CXX) $(FLAGS) -DUNITTEST $(INCLUDE)  $(WFLAGS) -c $< $
# 
# multigridtest: multigrid.o
# 	$(CXX) $(FLAGS)   $(WFLAGS) $< $(KASKADE7)/libs/umfpack_solve.o $(KASKADE7)/libs/mumps_solve.o $(KASKADELIB) $(DUNELIB) $(UGLIB) $(BOOSTLIB) $(UMFPACKLIB) \
#  $(MUMPSLIB) \
#  $(BLASLIB) $(FTNLIB) $(NUMALIB)  -o $@
# 
# 
