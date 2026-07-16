
TESTRESULTS += tests/elastomechanics3d/testResult.txt
CLEAN += $(patsubst %, tests/elastomechanics3d/%, elasto elasto.o)
DEPENDSRC += tests/elastomechanics3d/elasto.cpp


$(KASKADE7)/tests/elastomechanics3d/testResult.txt: $(KASKADE7)/tests/elastomechanics3d/elasto
	cd $(KASKADE7)/tests/elastomechanics3d; ./elasto --result 1 --verbosity 0 --iterative 1 

$(KASKADE7)/tests/elastomechanics3d/elasto: $(KASKADE7)/tests/elastomechanics3d/elasto.o kasklib
	$(CXX) $(FLAGS) $< $(KASKADE7)/libs/*.o $(KASKADELIB) \
	$(DUNELIB) \
	$(UGLIB) \
	$(BOOSTLIB) \
	$(DIRECTSOLVERLIB) \
	$(ITSOLLIB) \
	$(HYPRELIB) \
	$(AMIRALIB) \
	$(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@ 
