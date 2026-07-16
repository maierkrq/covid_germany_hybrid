#! /bin/bash
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $DOWNLOAD
wget -c $SUITESPARSE_SRC_PATH

cd $SCRATCH
TARFILE=`basename $SUITESPARSE_SRC_PATH`
tar xzf $DOWNLOAD/$TARFILE

# Build instructions according to SuiteSparse 5.3.0 README
cd SuiteSparse
make BLAS="$BLASLIB" LAPACK="$LAPACKLIB" library
make BLAS="$BLASLIB" LAPACK="$LAPACKLIB" INSTALL=$TARGET install

echo UMFPACKINC = -I$TARGET/include/             >> $TARGET/Makefile.Local
echo UMFPACKLIB = -L$TARGET/lib -lumfpack -lamd  >> $TARGET/Makefile.Local


#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


