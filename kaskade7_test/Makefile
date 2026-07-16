###################################################################################################
#
# Targets
#
# help
# Provides a list of targets
#
# kasklib (default)
# Builds the Kaskade 7 library.
#
# doc
# Creates the doxygen documentation for Kaskade 7
#
# clean
# Removes intermediate files.
#
# test
# Builds and executes the tests.
#
# tutorial
# Builds the tutorial examples.
#
# run-tutorial
# Builds and executs the tutorial examples.
#
#
# Variables
#
# Defined for use in Kaskade 7 programs, tests, and tutorials:
# ============================================================
#
# CXX
# The C++ compiler to use.
#
# INCLUDE
# Compiler flags specifying include directories for Kaskade and third party libraries
#
# WFLAGS
# Compiler flags for warnings
#
# OPTFLAGS, DEBUGFLAGS, OPTLOWFLAGS, QUICKCOMPILEFLAGS
# Compiler flags defining optimized or debug builds. Use only one of these.
#
# LINKFLAGS
# Compiler flags for linking. They specify the paths to search at run time for shared libraries.
#
#
# To be extended in Kaskade 7 submodules, tests, and tutorial makefiles:
# ======================================================================
#
# LIBOBJ
# Object files that shall end up in libkaskade. These are file paths relative to the Kaskade
# root directory.
#
# DEPENDSRC
# Source files for which dependencies shall be created. The source files for LIBOBJ
# do not need to be listed here, they are added automatically.
#
# CLEAN
# Files to remove on "make clean". LIBOBJ and TESTRESULTS files do not need to be listed here,
# they are added automatically. Files listed here shall be given relative to the KASKADE7 root directory.
#
# TESTRESULTS
# Output files of tests. Those are the dependencies of "make test" - make sure that their
# dependencies and rules are correct. The result files shall be given relative to the KASKADE7
# root directory.
#
# TUTORIALTARGETS
# Targets for tutorials. These are typically phony targets that build the tutorial programs.
#
# TUTORIALRUNTARGETS
# Targets for running the tutorials. These are typically phony targets that execute the tutorial programs.
#
###################################################################################################

.PHONY: default
default: kasklib

include Makefile.Local
include $(KASKADE7)/Makefile.Rules
include $(KASKADE7)/tools/Makefile



.PHONY: all
all: kasklib doc tests


###################################################################################################
#
# Provide help text
#
help:
	@echo Supported make targets:
	@echo "help:		        this help text"
	@echo "kasklib (default):	build the Kaskade 7 library"
	@echo "doc:			create the Doxygen documentation"
	@echo "manual:			build the LaTeX programmer\'s manual"
	@echo "test:			build and execute self test programs"
	@echo "tutorial:		build tutorial examples"
	@echo "run-tutorial:		run tutorial examples"
	@echo "clean:			remove intermediate files"


###################################################################################################
#
# Build the libkaskade library.
#

# The library consists of several "modules". Each module resides in a subdirectory and provides 
# its own list of object files and submodules. The modules are listed here.
LIBMODULES = algorithm fem io linalg mg timestepping utilities
LIBTARGET = $(KASKADE7)/libs/libkaskade.a
	
# Module makefiles add their objects to LIBOBJ if they end up in libkaskade. Files to be cleaned can be added to CLEAN.
# Additional source files for which to create dependency information can be added to DEPENDSRC.
include $(patsubst %, $(KASKADE7)/%/Makefile, $(LIBMODULES))

# Target for building the object files ending up in the Kaskade library.
# DEPRECATED - WILL BE REMOVED IN FUTURE VERSIONS
.PHONY: objects
objects: $(patsubst %, $(KASKADE7)/%, $(LIBOBJ))

# Target for building the Kaskade 7 library. This includes libkaskade.a and the target objects.
# The intermediate object files in the module directories are added as prerequisites in order to 
# prevent automatic deletion of intermediate files.
.PHONY: kasklib
kasklib: $(LIBTARGET)

$(LIBTARGET): $(patsubst %, $(KASKADE7)/%, $(LIBOBJ))
	mkdir -p $(KASKADE7)/libs
	$(AR)cr $(LIBTARGET) $(patsubst %, $(KASKADE7)/%, $(LIBOBJ))


###################################################################################################
#
# Build and execute the tests
#
# Test modules shall add to the TESTRESULTS variable the output files of their tests,
# and define rules for creating the test results. They shall also add source files to the
# DEPENDSRC variable and intermediate files to the CLEAN variable, if necessary.

			
TESTMODULES = barproblem elastomechanics3d embErrorEst hbPrec poisson reentrantCorner statHeatTransfer

-include $(patsubst %, $(KASKADE7)/tests/%/Makefile, $(TESTMODULES))

CLEAN += tests/testResult.txt $(TESTRESULTS)

.PHONY: test
test: $(KASKADE7)/Makefile.Local $(patsubst %, $(KASKADE7)/%, $(TESTRESULTS))
	cat $(patsubst %, $(KASKADE7)/%, $(TESTRESULTS)) > $(KASKADE7)/tests/testResult.txt


###################################################################################################
#
# Build and execute the tutorials
#
# Tutorial modules shall add to the TUTORIALTARGETS variable. They shall also add source files to the
# DEPENDSRC variable and intermediate files to the CLEAN variable, if necessary.

TUTORIALMODULES = covid_germany_abm_pde_ode #paper_covid_abm_pde_Berlin_Brandenburg #covid_germany_abm #
#covid_test_abm_pde_Berlin_Brandenburg #covid_test_abm_Berlin_Brandenburg #covid_abm_pde_Berlin_Brandenburg #covid_abm_Berlin_Brandenburg #covid_pde_ode_faster #covid_pde_separate_Berlin_Brandenburg #
#covid_frac_pde_separate_Berlin_Brand #
include $(patsubst %, $(KASKADE7)/tutorial/%/Makefile.gen, $(TUTORIALMODULES))


tutorial: Makefile.Local $(TUTORIALTARGETS)

.PHONY: run-tutorial
run-tutorial: Makefile.Local $(TUTORIALRUNTARGETS)

###################################################################################################
#
# Build the Doxygen documentation.
#

.PHONY: doc manual
doc: $(KASKADE7)/doc/html/index.html

$(KASKADE7)/doc/html/index.html: $(shell find $(KASKADE7)/ -name "*.hh") $(KASKADE7)/doc/macros.js \
				 $(KASKADE7)/doc/Doxyfile $(KASKADE7)/doc/kaskade7_test.css
	cd $(KASKADE7)/doc; sh makeDocu

manual: $(KASKADE7)/tutorial/ProgrammersManual/manual.pdf

$(KASKADE7)/tutorial/ProgrammersManual/manual.pdf: $(KASKADE7)/tutorial/ProgrammersManual/manual.tex
	cd $(KASKADE7)/tutorial/ProgrammersManual; pdflatex manual.tex; pdflatex manual.tex



###################################################################################################
#
# Build dependency information
#
# For each source file *.cpp we maintain a dependency file *.d. Source files are listed in
# DEPENDSRC (library sources end up there as defined here, but modules can define additional
# source files). All dependency files are included, which means they are built by make if
# they do not yet exist. They are updated together with their source files due to a dependency
# in the .cpp => .o compilation rule.


#
# WARNING: DO NOT INCLUDE MODULE MAKEFILES *BELOW* THIS LINE! (ONLY ABOVE)
#          THE REASON IS, THAT MODULE MAKEFILES TYPICALLY ADD TO THE DEPENDSRC VARIABLE
#          SUCH THAT MISSING DEPENDENCY FILES ARE BUILT AUTOMATICALLY. THIS HAPPENS HERE,
#          AND THE CONTENTS OF DEPENDSRC IS USED AT THIS POSITION. IF A MODULE MAKEFILE
#          IS INCLUDED LATER, THE CHANGES TO DEPENDSRC WILL HAVE NO EFFECT, AND DEPENDENCY
#          FILES REMAIN MISSING. THIS, IN TURN, PREVENTS MAKE FROM FINDING THE APPROPRIATE
#          RULE FOR BUILDING OBJECT FILES...
#
DEPENDSRC += $(LIBOBJ:.o=.cpp)
include $(DEPENDSRC:.cpp=.d)



CLEAN += gccerr.txt

.PHONY: clean
clean:
	rm -f $(LIBOBJ) $(patsubst %, $(KASKADE7)/%, $(CLEAN)) $(DEPENDSRC:.cpp=.d)
	rm -f remove_install_deps *.bak $(KASKADE7)/libs/*.o $(KASKADE7)/libs/*.a


KASKADE7NAME = $(shell basename $(KASKADE7))

