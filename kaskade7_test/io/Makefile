

LIBOBJ += $(patsubst %, io/%, readPoly.o iobase.o ansysmeshreader.o base64.o vtkwriter.o checkPoint.o vtkreader.o check_endianness.o)



# include ../Makefile.Local
# include ../Makefile.Rules
# 
# INCLUDE = $(DUNEINC) $(UGINC) $(FUSIONINC) $(BOOSTINC) $(UMFPACKINC) $(KASKADEINC) 
# FLAGS = $(WFLAGS) $(OPTFLAGS)
# FLAGS2 = $(WFLAGS) -O1 -std=c++0x -DBOOST_NO_VARIADIC_TEMPLATES -DZIBHTC
# FLAGS3 = $(WFLAGS) $(DEBUGFLAGS)
# 
# 
# OBJ = readPoly.o iobase.o ansysmeshreader.o base64.o vtkwriter.o checkPoint.o vtkreader.o
# SRC = $(OBJ:%.o=%.cpp)
# OBJ2 = check_endianness.o
# SRC2 = check_endianness.cpp
# 
# TARGETLIB = ../libs/libkaskade.a
# 
# # default target ist to build (update) libkaskade, and the check for endianness
# default: $(TARGETLIB) $(OBJ2)
# 
# # the libkaskade library depends on its member object files (in particular the ones built here)
# $(TARGETLIB): $(TARGETLIB)($(OBJ)) $(TARGETLIB)($(OBJ2)) Makefile
# 
# $(OBJ2): $(SRC2)
# 	$(CXX) $(FLAGS2) -c $(SRC2)
# 
# 
# depend:
# 	cp Makefile.gen Makefile
# 	../tools/gccmakedep -f Makefile $(INCLUDE) $(MAKEDEPENDFLAGS) $(SRC) $(SRC2)
# 	../tools/remove_install_deps Makefile $(KASKADE7)
# 
# clean:
# 	rm -f gccerr.txt $(OBJ) $(OBJ2) *.bak *.tmp
# 
# 	
# vtktest: vtkwriter2.cpp vtkwriter2.hh
# 	$(CXX) $(FLAGS3) $(INCLUDE) -DUNITTEST $< $(KASKADELIB) $(BOOSTLIB) $(DUNELIB) $(UGLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB) -o $@
# 
# checkPointTest: checkPoint.cpp checkPoint.hh
# 	$(CXX) $(FLAGS3) $(INCLUDE) -DUNITTEST $< $(KASKADELIB) $(BOOSTLIB) -lboost_filesystem $(DUNELIB) $(UGLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB) -o $@
# 
# vtkreaderTest: vtkreader.cpp vtkreader.hh
# 	$(CXX) $(FLAGS) $(INCLUDE) -DUNITTEST $< $(KASKADELIB) $(BOOSTLIB) $(DUNELIB) $(UGLIB) $(BLASLIB) $(FTNLIB) $(NUMALIB) -o $@
# 
# 
# 	
# # DO NOT DELETE
# 
