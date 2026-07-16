#! /bin/sh
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $DOWNLOAD
wget -c $HYPRE_SRC_PATH

cd $SCRATCH
TARFILE=`basename $HYPRE_SRC_PATH`
tar xzf $DOWNLOAD/$TARFILE

# Build instructions according to Hypre 2.11.2
DIR=`basename -s .tar.gz $TARFILE`
DIR=`basename -s .tgz $DIR`
cd $DIR/src 

(export CC CXX FC; ./configure --prefix=$TARGET --without-MPI --with-blas-lib="$BLASLIB" --with-lapack-lib="$LAPACKLIB")
make install

echo HYPREINC = -I$TARGET/include/     >> $TARGET/Makefile.Local
echo HYPRELIB = -L$TARGET/lib -lHYPRE  >> $TARGET/Makefile.Local



#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


