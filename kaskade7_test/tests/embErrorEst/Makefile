
TESTRESULTS += tests/embErrorEst/testResult.txt
CLEAN += $(patsubst %, tests/embErrorEst/%, embErrorEst embErrorEst.o)
DEPENDSRC += tests/embErrorEst/embErrorEst.cpp


$(KASKADE7)/tests/embErrorEst/testResult.txt: $(KASKADE7)/tests/embErrorEst/embErrorEst
	cd $(KASKADE7)/tests/embErrorEst; ./embErrorEst --verbosity 0 --result 1 --solver.direct MUMPS

$(KASKADE7)/tests/embErrorEst/embErrorEst: $(KASKADE7)/tests/embErrorEst/embErrorEst.o kasklib
	$(CXX) $< $(KASKADE7)/libs/*.o $(KASKADELIB) \
 $(DUNELIB) \
 $(UGLIB) \
 $(BOOSTLIB) \
 $(DIRECTSOLVERLIB) \
 $(ITSOLLIB) \
 $(HYPRELIB) \
 $(AMIRALIB) \
 $(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@












