#! /bin/sh
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $DOWNLOAD
wget -c $ITSOL_SRC_PATH

cd $SCRATCH
TARFILE=`basename $ITSOL_SRC_PATH`
tar xzf $DOWNLOAD/$TARFILE

# Build instructions according to ITSOL 2
DIR=`basename -s .tar.gz $TARFILE`
DIR=`basename -s .tgz $DIR`
cd $DIR

# Edit the makefile. Probably all we need to do is to specify the desired
# C and Fortran compilers and add symbol table generation to ar.
sed -e "s|AR[[:blank:]]*=[[:blank:]]*ar|AR = ar -s|g" \
    -e "s|CC[[:blank:]]*=[[:blank:]]*gcc|CC = $CC|g" \
    -e "s|FC[[:blank:]]*=[[:blank:]]*gfortran|FC = $FC|g" makefile > makefile.sed
mv makefile.sed makefile

make lib
cp INC/*.h $TARGET/include/
cp LIB/lib*.a $TARGET/lib/

# ITSOL's ios.h defines HB as 1, which is a pretty bad idea colliding with Kaskade's HB
# preconditioner type.
sed -e "s/#define HB   1/#define ITSOL_HB 1/g" INC/ios.h > $TARGET/include/ios.h

echo ITSOLINC = -I$TARGET/include/     >> $TARGET/Makefile.Local
echo ITSOLLIB = -L$TARGET/lib -litsol  >> $TARGET/Makefile.Local




#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


