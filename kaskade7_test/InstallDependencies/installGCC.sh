#! /bin/sh
# Copyright 2018 - 2020 Zuse Institute Berlin




GCC_CONFIGURE=

#------------------------------------------------------------------------------
# Install GMP if required. This library is a prerequisite of MPC, but
# if there is a system library, available, we can use that.
if [ "$GMP_SRC_PATH" != "" ]; then
  cd $DOWNLOAD
  echo downloading GMP library from $GMP_SRC_PATH >> $LOGFILE
  wget -c $GMP_SRC_PATH

  cd $SCRATCH
  TARFILE=$DOWNLOAD/`basename $GMP_SRC_PATH`
  tar xf $TARFILE

  SRCDIR=`basename -s .gz $TARFILE`
  SRCDIR=`basename -s .xz $SRCDIR`
  SRCDIR=`basename -s .tar $SRCDIR`
  SRCDIR=$SCRATCH/`basename -s .tgz $SRCDIR`

  cd $SRCDIR
  echo configuring GMP >> $LOGFILE
  ./configure --prefix=$TARGET/GCC

  echo building GMP >> $LOGFILE
  make
  if [ $? -ne 0 ]; then
    echo BUILDING GMP FAILED | tee -a $LOGFILE
    exit
  fi
  make install
fi

#------------------------------------------------------------------------------
# Install MPFR if required. This library is a prerequisite of MPC, but
# if there is a system library, available, we can use that.
if [ "$MPFR_SRC_PATH" != "" ]; then
  cd $DOWNLOAD
  echo downloading MPFR library from $MPFR_SRC_PATH >> $LOGFILE
  wget -c $MPFR_SRC_PATH

  cd $SCRATCH
  TARFILE=$DOWNLOAD/`basename $MPFR_SRC_PATH`
  tar xf $TARFILE

  SRCDIR=`basename -s .gz $TARFILE`
  SRCDIR=`basename -s .xz $SRCDIR`
  SRCDIR=`basename -s .tar $SRCDIR`
  SRCDIR=$SCRATCH/`basename -s .tgz $SRCDIR`

  cd $SRCDIR
  echo configuring MPFR >> $LOGFILE
  ./configure --prefix=$TARGET/GCC --with-gmp=$TARGET/GCC

  echo building MPFR >> $LOGFILE
  make
  if [ $? -ne 0 ]; then
    echo BUILDING MPFR FAILED | tee -a $LOGFILE
    exit
  fi
  make install
fi

#------------------------------------------------------------------------------
# Install MPC if required. This library is a prerequisite of GCC, but
# if there is a system library, available, we can use that.
if [ "$MPC_SRC_PATH" != "" ]; then
  cd $DOWNLOAD
  echo downloading MPC library from $MPC_SRC_PATH >> $LOGFILE
  wget -c $MPC_SRC_PATH

  cd $SCRATCH
  TARFILE=$DOWNLOAD/`basename $MPC_SRC_PATH`
  tar xf $TARFILE

  SRCDIR=`basename -s .gz $TARFILE`
  SRCDIR=`basename -s .tar $SRCDIR`
  SRCDIR=$SCRATCH/`basename -s .tgz $SRCDIR`

  cd $SRCDIR
  echo configuring MPC
  ./configure --prefix=$TARGET/GCC --with-gmp=$TARGET/GCC

  echo building MPC >> $LOGFILE
  make
  if [ $? -ne 0 ]; then
    echo BUILDING MPC FAILED | tee -a $LOGFILE
    exit
  fi
  make install
  GCC_CONFIGURE="$GCC_CONFIGURE --with-mpc=$TARGET/GCC"
fi

#------------------------------------------------------------------------------
# Install GCC itself.

cd $DOWNLOAD
echo downloading GCC from $GCC_SRC_PATH >> $LOGFILE
wget -c $GCC_SRC_PATH

cd $SCRATCH
TARFILE=$DOWNLOAD/`basename $GCC_SRC_PATH`
tar xf $TARFILE

# remove popular archive compression suffixes
SRCDIR=`basename -s .xz $TARFILE`
SRCDIR=`basename -s .gz $SRCDIR`
SRCDIR=`basename -s .tar $SRCDIR`
SRCDIR=$SCRATCH/`basename -s .tgz $SRCDIR`

OBJDIR=$SCRATCH/objdir
mkdir $OBJDIR

# GCC configuration (consider --diable-bootstrap for faster installation)
cd $OBJDIR
echo configuring GCC >> $LOGFILE
$SRCDIR/configure --prefix=$TARGET --disable-multilib $GCC_CONFIGURE

# GCC building (consider executing the self test as well)
echo building GCC >> $LOGFILE
make -j 4
if [ $? -ne 0 ]; then
  echo BUILDING GCC FAILED: status $? | tee -a $LOGFILE
  exit
fi

# GCC installation
echo installing GCC >> $LOGFILE
make install

# Report C++ compiler to use
CXX=$TARGET/bin/g++
CC=$TARGET/bin/gcc
FC=$TARGET/bin/gfortran
FTNLIB="-lgfortran -lpthread"


#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi



