
TESTRESULTS += tests/barproblem/testResult.txt
CLEAN += $(patsubst %, tests/barproblem/%, bar bar.o)
DEPENDSRC += tests/barproblem/bar.cpp


$(KASKADE7)/tests/barproblem/testResult.txt: $(KASKADE7)/tests/barproblem/bar
	cd $(KASKADE7)/tests/barproblem; ./bar --result 1 --quick 1 --verbosity 0

$(KASKADE7)/tests/barproblem/bar: $(KASKADE7)/tests/barproblem/bar.o kasklib
	$(CXX) $(FLAGS) $< $(KASKADE7)/libs/*.o $(KASKADELIB) \
	$(DUNELIB) \
	$(UGLIB) \
	$(BOOSTLIB) \
	$(DIRECTSOLVERLIB) \
	$(ITSOLLIB) \
	$(HYPRELIB) \
	$(AMIRALIB) \
	$(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@ 
