
TESTRESULTS += tests/reentrantCorner/testResult.txt
CLEAN += $(patsubst %, tests/reentrantCorner/%, reentrantCorner reentrantCorner.o)
DEPENDSRC += tests/reentrantCorner/reentrantCorner.cpp


$(KASKADE7)/tests/reentrantCorner/testResult.txt: $(KASKADE7)/tests/reentrantCorner/reentrantCorner
	cd $(KASKADE7)/tests/reentrantCorner; ./reentrantCorner --verbosity 0 --result 1

$(KASKADE7)/tests/reentrantCorner/reentrantCorner: $(KASKADE7)/tests/reentrantCorner/reentrantCorner.o kasklib
	$(CXX) $< $(KASKADE7)/libs/*.o $(KASKADELIB) \
 $(DUNELIB) \
 $(UGLIB) \
 $(BOOSTLIB) \
 $(DIRECTSOLVERLIB) \
 $(ITSOLLIB) \
 $(HYPRELIB) \
 $(AMIRALIB) \
 $(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@



