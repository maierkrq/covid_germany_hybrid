#! /bin/bash
# Copyright 2019 - 2019 Zuse Institute Berlin

cd $DOWNLOAD
echo Downloading boost >> $LOGFILE
wget -c $BOOST_SRC_PATH

cd $SCRATCH
TARFILE=`basename $BOOST_SRC_PATH`
tar xjf $DOWNLOAD/$TARFILE

# Build instructions according to Boost getting started on Unix
cd `basename -s .tar.bz2 $BOOST_SRC_PATH`
echo Building boost >> $LOGFILE
bash bootstrap.sh --prefix=$TARGET
./b2 install


echo BOOSTINC = -I$TARGET/include/             >> $TARGET/Makefile.Local
echo BOOSTLIB = -L$TARGET/lib                  >> $TARGET/Makefile.Local


#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


