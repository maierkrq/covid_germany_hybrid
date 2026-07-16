
LIBOBJ += $(patsubst %, utilities/%, kaskopt.o threading.o enums.o kalloc.o date_and_time.o save_file.o timing.o rotate.o deprecation.o gaussNodes.o)



# include ../Makefile.Local
# include ../Makefile.Rules
# 
# INCLUDE = $(BOOSTINC) $(KASKADEINC) $(DUNEINC)
# FLAGS = $(WFLAGS) $(OPTFLAGS)  -DKASKADE=\"$(KASKADE7)\" -DVERSION=\"0000\"
# PLATFORM= $(shell uname)
# 
# OBJ = kaskopt.o threading.o enums.o kalloc.o date_and_time.o save_file.o timing.o rotate.o deprecation.o gaussNodes.o
# SRC = $(OBJ:.o=.cpp)
# 
# 
# # default target ist to build (update) libkaskade
# TARGETLIB = ../libs/libkaskade.a
# default: $(TARGETLIB)
# 
# # the libkaskade library depends on its member object files (in particular the ones built here)
# $(TARGETLIB): $(TARGETLIB)($(OBJ)) Makefile
# 
# depend:
# 	cp Makefile.gen Makefile
# 	../tools/gccmakedep -f Makefile $(INCLUDE) $(MAKEDEPENDFLAGS) $(SRC)
# 	../tools/remove_install_deps Makefile $(KASKADE7)
# 
# clean:
# 	rm -f gccerr.txt ptree $(OBJ) *.bak run*
# 
# # unit test targets
# ptree: ptree.o kaskopt.o Makefile
# 	$(CXX) -g $<  \
#  $(KASKADELIB) \
#  $(BOOSTLIB) \
#  -o $@
# 
# 
# 
