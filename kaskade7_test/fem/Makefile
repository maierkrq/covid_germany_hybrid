
# include module makefiles to add their targets and dependencies
FEMMODULES ::= fem/diffops

include $(patsubst %, %/Makefile, $(FEMMODULES))


LIBOBJ += $(patsubst %, fem/%, assemble.o boundaryFace.o boundaryLocator.o cellLocator.o boundaryInterpolation.o errorest.o hierarchicshapefunctions.o lagrangeshapefunctions.o nedelecshapefunctions.o )

# unit test targets

fem/linearspacetest: fem/linearspace.cpp
	$(CXX) $(FLAGS) $(INCLUDE) -DUNITTEST $< $(KASKADELIB) \
	$(DUNELIB) \
	$(BOOSTLIB) \
	$(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@


# hdivshapefunctiontest: nedelecshapefunctions.cpp nedelecshapefunctions.hh
# 	$(CXX) -g -DUNITTEST $(INCLUDE)  $(WFLAGS) $(DEBUGFLAGS) $< $(KASKADELIB) $(DUNELIB)  $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB)  -o $@
# 	
# advectiontest: advect.cpp advect.hh
# 	$(CXX) -g -DUNITTEST $(INCLUDE)  $(WFLAGS) $(DEBUGFLAGS) $< $(KASKADELIB) $(DUNELIB) $(UGLIB)  $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB)  -o $@
# 	
# intersectiontest: boundaryLocator.cpp boundaryLocator.hh
# 	$(CXX) -g -DUNITTEST $(INCLUDE)  $(WFLAGS) $(DEBUGFLAGS) $< $(KASKADELIB) $(DUNELIB) $(UGLIB)  $(BOOSTLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB)  -o $@


