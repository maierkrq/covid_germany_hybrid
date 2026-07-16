

TESTRESULTS += tests/statHeatTransfer/testResult.txt
CLEAN += $(patsubst %, tests/statHeatTransfer/%, statHeatTransfer statHeatTransfer.o)
DEPENDSRC += tests/statHeatTransfer/statHeatTransfer.cpp


$(KASKADE7)/tests/statHeatTransfer/testResult.txt: $(KASKADE7)/tests/statHeatTransfer/statHeatTransfer
	cd tests/statHeatTransfer; ./statHeatTransfer --verbosity 0 --result 1 --solver.type direct

$(KASKADE7)/tests/statHeatTransfer/statHeatTransfer: tests/statHeatTransfer/statHeatTransfer.o kasklib
	$(CXX) $< libs/*_solve.o $(KASKADELIB) \
 $(DUNELIB) \
 $(UGLIB) \
 $(BOOSTLIB) \
 $(DIRECTSOLVERLIB) \
 $(ITSOLLIB) \
 $(HYPRELIB) \
 $(AMIRALIB) \
 $(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@















