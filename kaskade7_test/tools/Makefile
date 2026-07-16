TOOLSRC ::= filterOutput.cpp sumDiffFile.cpp remove_install_deps.cpp
DEPENDSRC += $(patsubst %, tools/%, $(TOOLSRC))
CLEAN += $(patsubst %.cpp, tools/%.o, $(TOOLSRC) Makefile.Sentinel)

$(KASKADE7)/tools/remove_install_deps: tools/remove_install_deps.cpp
	$(CXX) $< -o $@


# We need to make sure that remove_install_deps is built first, as it is
# needed for the dependency file generation. Here's the trick: Inclusion of
# the dummy file Makefile.Sentinel triggers its building, and as it depends on
# remove_install_deps *explicitly* and not by a pattern rule, it will be built
# - and its dependency remove_install_deps as well.
include $(KASKADE7)/tools/Makefile.Sentinel

$(KASKADE7)/tools/Makefile.Sentinel: $(KASKADE7)/tools/remove_install_deps
	touch $@


