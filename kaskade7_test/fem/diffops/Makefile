
LIBOBJ += $(patsubst %, fem/diffops/%, contactConstraints.o elasto.o membraneModels.o elastoVariationalFunctionals.o \
                                       materialLaws.o shellKinematics.o)



# unit tests

DEPENDSRC += fem/diffops/shelltest.cpp

fem/diffops/shelltest: fem/diffops/shelltest.o kasklib
	$(CXX) -g -DUNITTEST $(FLAGS) $< $(KASKADELIB) $(DUNELIB) $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@

# diffops/elastoTest.o: diffops/elastoVariationalFunctionals.cpp diffops/elastoVariationalFunctionals.hh
# 	$(CXX) $(FLAGS) $(INCLUDE) -DUNITTEST $< -c -o $@ 2> gccerr.txt 
# 
# diffops/elasto: diffops/elastoTest.o
# 	$(CXX) $<  $(KASKADELIB) \
# 	    $(BOOSTLIB) \
# 	    $(DUNELIB) \
# 	    $(BLASLIB) $(FTNLIB) -o $@ 



	
