
TESTRESULTS += tests/hbPrec/testResult.txt
CLEAN += $(patsubst %, tests/hbPrec/%, hbPrec hbPrec.o)
DEPENDSRC += tests/hbPrec/hbPrec.cpp


$(KASKADE7)/tests/hbPrec/testResult.txt: $(KASKADE7)/tests/hbPrec/hbPrec
	cd $(KASKADE7)/tests/hbPrec; ./hbPrec --verbosity 0 --result 1

$(KASKADE7)/tests/hbPrec/hbPrec: $(KASKADE7)/tests/hbPrec/hbPrec.o kasklib
	$(CXX) $< $(KASKADE7)/libs/*.o $(KASKADELIB) \
 $(DUNELIB) \
 $(UGLIB) \
 $(BOOSTLIB) \
 $(DIRECTSOLVERLIB) \
 $(ITSOLLIB) \
 $(HYPRELIB) \
 $(AMIRALIB) \
 $(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@

