#! /bin/sh
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $SCRATCH

downloadAndUnpack () {
  FROM=$1
  cd $DOWNLOAD
  wget -c $FROM
  cd $SCRATCH
  tar xzf $DOWNLOAD/`basename $FROM`
}

# get the DUNE core modules
downloadAndUnpack $DUNE_SRC_PATH/dune-common-$DUNE_VERSION.tar.gz
downloadAndUnpack $DUNE_SRC_PATH/dune-geometry-$DUNE_VERSION.tar.gz
downloadAndUnpack $DUNE_SRC_PATH/dune-grid-$DUNE_VERSION.tar.gz
downloadAndUnpack $DUNE_SRC_PATH/dune-istl-$DUNE_VERSION.tar.gz

# get further modules
git clone https://gitlab.dune-project.org/staging/dune-uggrid.git --branch releases/$DUNE_BRANCH
git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git --branch releases/$DUNE_BRANCH


# building DUNE
# QUESTION: HOW TO CONVINCE CMAKE TO FIND MKL BLAS/LAPACK LIBRARIES INSTALLED IN A NONSTANDARD DIRECTORY?
#           FOR NOW WE LET DUNE USE THE GENERIC BLAS/LAPACK THAT CMAKE FINDS, AND PROVIDE OUR OWN
#           WHEN LINKING KASKADE EXECUTABLES
echo -n "CMAKE_FLAGS=\"" > opts
echo -n "-DCMAKE_CXX_COMPILER=$CXX -DCMAKE_INSTALL_PREFIX=$TARGET -DUMFPACK_ROOT=$TARGET "  >> opts
echo \" >> opts

# DUNE uggrid builds parallel (MPI) version of UG if cmake finds an MPI library. Here we 
# try to convince UG to build a sequential version even if an MPI library is available.
sed -e "s/find_package(MPI)/set(MPI_C_FOUND 0)/g" < dune-uggrid/CMakeLists.txt > out
mv out dune-uggrid/CMakeLists.txt

./dune-common-$DUNE_VERSION/bin/dunecontrol --opts=$PWD/opts all
./dune-common-$DUNE_VERSION/bin/dunecontrol --opts=$PWD/opts make install

# The config.h include files are not installed by the DUNE installer script, but we need them in Kaskade.
cp dune-common-$DUNE_VERSION/build-cmake/config.h $TARGET/include/dune/common/
cp dune-grid-$DUNE_VERSION/build-cmake/{config,FC}.h $TARGET/include/dune/grid/

# There is an old bug in multilevellineargeometry.hh that we need to patch out.
patch $TARGET/include/dune/geometry/multilineargeometry.hh $INSTALLHOME/patches/multilineargeometry.hh.patch


#------------------------------------------------------------------------------
# Enter paths into the Makefile. Note that DUNELIB must come before UGLIB, but UG needs dunegeometry behind it.
echo DUNEINC = -I$TARGET/include -I$TARGET/include/dune       >> $TARGET/Makefile.Local
echo DUNELIB = -L$TARGET/lib -ldunegeometry -ldunegrid -ldunecommon >> $TARGET/Makefile.Local
echo ALUGRIDINC = -I$TARGET/include/dune/alugrid              >> $TARGET/Makefile.Local
echo ALUGRIDLIB = -ldunealugrid                               >> $TARGET/Makefile.Local
echo UGINC = -I$TARGET/include/ug                             >> $TARGET/Makefile.Local
echo UGLIB = -lugS2 -lugS3 -lugL -ldunegeometry -lmpi         >> $TARGET/Makefile.Local
echo HAVE_UG = -DHAVE_UG=1 -DENABLE_UG=1 -DUG_USE_NEW_DIMENSION_DEFINES -DUG_USE_SYSTEM_HEAP=1 -DUGLIB  >> $TARGET/Makefile.Local

#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


