#! /bin/sh
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $DOWNLOAD
TARFILE=`basename $LIBAMIRAMESH_SRC_PATH`
curl -C -  -o $TARFILE $LIBAMIRAMESH_SRC_PATH

cd $SCRATCH
tar xzf $DOWNLOAD/$TARFILE

# The provided Makefile uses $(_path_gcc)/bin/$(_ccomp) as C compiler.
# We extract the path components from the CXX path.
export _path_gcc=`echo $CXX | sed -e "s?/bin/.*??g"`
export _ccomp=`echo $CC | sed -e "s?.*/bin/??g"`
export _cxxcomp=`echo $CXX | sed -e "s?.*/bin/??g"`
(cd libamiramesh/src; make)
cp -r libamiramesh/include/amiramesh $TARGET/include/
cp libamiramesh/lib/libamiramesh.a $TARGET/lib/

echo AMIRAINC = -I$TARGET/include/         >> $TARGET/Makefile.Local
echo AMIRALIB = -L$TARGET/lib -lamiramesh  >> $TARGET/Makefile.Local



#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


