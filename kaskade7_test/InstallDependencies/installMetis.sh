#! /bin/sh
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $DOWNLOAD
wget -c $METIS_SRC_PATH

cd $SCRATCH
TARFILE=`basename $METIS_SRC_PATH`
tar xzf $DOWNLOAD/$TARFILE

# Build instructions according to ITSOL 2
DIR=`basename -s .tar.gz $TARFILE`
DIR=`basename -s .tgz $DIR`
cd $DIR

# Configuration according to Metis 5.1.0 BUILD.TXT
make config prefix=$TARGET cc=$CC
make
make install


#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


