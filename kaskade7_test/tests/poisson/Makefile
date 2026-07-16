
TESTRESULTS += tests/poisson/testResult.txt
CLEAN += $(patsubst %, tests/poisson/%, poisson poisson.o)
DEPENDSRC += tests/poisson/poisson.cpp


$(KASKADE7)/tests/poisson/testResult.txt: $(KASKADE7)/tests/poisson/poisson
	cd tests/poisson; ./poisson --verbosity 0 --result 1

$(KASKADE7)/tests/poisson/poisson: tests/poisson/poisson.o kasklib
	$(CXX) $< libs/*.o $(KASKADELIB) \
	$(DUNELIB) \
	$(UGLIB) \
	$(BOOSTLIB) \
	$(DIRECTSOLVERLIB) \
	$(ITSOLLIB) \
	$(HYPRELIB) \
	$(AMIRALIB) \
	$(BLASLIB) $(FTNLIB) $(NUMALIB) \
	$(LINKFLAGS) -o $@

