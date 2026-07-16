#! /bin/bash
# Copyright 2018 - 2019 Zuse Institute Berlin

# -----------------------------------------------------------------------------
# Configuration section
# You will probably need to set paths and package sources here.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Installation target
# This specifies the path where to install the third party libraries. You need to
# have write access to that directory. ON INSTALLATION, PREVIOUSLY EXISTING
# CONTENTS IN THAT DIRECTORY WILL BE DELETED.
# TARGET=$(KASKADE7)/work/input/installed
TARGET=/data/numerik/software/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed

# -----------------------------------------------------------------------------
# Download directory
# If the files to download are already available there, they'll not be fetched again.
DOWNLOAD=/data/numerik/software/KaskadeDependencies/Kaskade7.5Dependencies-10.2/download

# -----------------------------------------------------------------------------
# Scratch space
# In this directory, temporary builds may be stored. AFTER INSTALLATION, THIS WILL BE DELETED.
SCRATCH=/tmp/KaskadeDependenciesScratch


# -----------------------------------------------------------------------------
# GCC
# Kaskade 7 requires a decent C++ compiler. GCC 8.x is known to work,
# but newer versions or other compilers may work as well. This installs
# GCC from source as specified in the following lines. It may be easier and
# quicker to rely on the compiler that comes with your distribution.

# The location from where to fetch the GCC source. If empty, GCC is not installed
# and the standard C++ compiler "c++" is used. The downloaded file shall be
# a gzip-ed tar archive.
# GCC_SRC_PATH=
GCC_SRC_PATH=ftp://ftp.fu-berlin.de/unix/languages/gcc/releases/gcc-10.2.0/gcc-10.2.0.tar.xz

# The location of the MPC library required by GCC. If your distribution has MPC already
# installed, you may set this as an empty string. The downloaded file shall be
# a gzip-ed tar archive.
MPC_SRC_PATH=https://ftp.gnu.org/gnu/mpc/mpc-1.1.0.tar.gz

# The location of the GMP library required by GCC. If your distribution has GMP already
# installed, you may set this as an empty string.
GMP_SRC_PATH=https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz

# The location of the MPFR library required by GCC. If your distribution has MPFR already
# installed, you may set this as an empty string.
MPFR_SRC_PATH=https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.xz

# -----------------------------------------------------------------------------
# boost
BOOST_SRC_PATH=https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2

# -----------------------------------------------------------------------------
# BLAS & LAPACK
# Right now, installation of BLAS and LAPACK is not supported here. Viable options are
# Atlas or MKL. Make the libraries to be used known here, e.g., for SuiteSparse.
# Use BLASLIB="-L/path/to/library/directory/ -lblas" and LAPACKLIB="-L/path/to/lib/ -llapack".

# Example: Intel MKL
# Leave MKL_ROOT empty if you don't want to use MKL
# MKL_ROOT=/home/numerik/bzfweise/numerik/htc/libraries/MKL/compilers_and_libraries_2019.1.144/linux/mkl/lib/intel64
MKL_ROOT=$TARGET/../../MKL/mkl/lib/intel64

# These are required for SuiteSparse make finding MKL
BLASLIBDIR=$MKL_ROOT
BLASLIB="-L$BLASLIBDIR -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
LAPACKLIB=$BLASLIB



# Example: Atlas
# Leave ATLAS_SRC_PATH empty if you don't want to use Atlas
ATLAS_SRC_PATH=

# -----------------------------------------------------------------------------
# SuiteSparse with UMFPACK direct solver and including Metis graph partitioning
SUITESPARSE_SRC_PATH=http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.3.0.tar.gz

# -----------------------------------------------------------------------------
# MUMPS - direct multifrontal solver 
# This uses Metis and is therefore installed *after* SuiteSparse
MUMPS_SRC_PATH=http://mumps.enseeiht.fr/MUMPS_5.3.3.tar.gz

# -----------------------------------------------------------------------------
# ITSOL - iterative solvers 
ITSOL_SRC_PATH=http://www-users.cs.umn.edu/~saad/software/ITSOL/ITSOL_2.tar.gz

# -----------------------------------------------------------------------------
# Hypre - preconditioners
HYPRE_SRC_PATH=https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.11.2.tar.gz

# -----------------------------------------------------------------------------
# libamiramesh
# This is a library that is not publically available.
LIBAMIRAMESH_SRC_PATH=file://$TARGET/../libamiramesh.tar.gz

# -----------------------------------------------------------------------------
# DUNE
# Kaskade 7 uses Dune, currently in version 2.7.

# The version of DUNE packages
DUNE_VERSION=2.7.0

# The locations from which to download the required DUNE modules.
DUNE_SRC_PATH=https://dune-project.org/download/$DUNE_VERSION

# Some grid modules, in particular alugrid and uggrid, do no longer offer packages for
# download. Thus, we check out the release branches from the git repository.
DUNE_BRANCH=2.7

# -----------------------------------------------------------------------------
# Florian's FP library
# Give the path to the git repository
FP_SRC_PATH=https://github.com/flwende/fp.git







# -----------------------------------------------------------------------------
# Execution section
# Conceptually, you should not need to modify anything below. Reality may be
# different, though, in particular in different environments or if package
# versions change, you may need to touch the installation logic below.
# -----------------------------------------------------------------------------

LOGFILE=$TARGET/install.log
rm -f $LOGFILE

if [ -n "$CXX" -o -n "$CC" ]; then
  echo Environment variables CC or CXX set:
  echo CC = $CC
  echo CXX = $CXX
  echo Wrongly set compilers can lead to all kinds of strange errors when installing GCC.
  echo Check the paths to make sure they are as intended.
  sleep 20
fi

# current path - here we expect to find the installation scripts
INSTALLHOME=$PWD

if [ ! -f "$INSTALLHOME/installGCC.sh" ]; then
  echo "Error: installGCC.sh not found in current working directory."
  echo "Please execute the install script from the directory where it is located."
  exit 1
fi

# Create empty target directory
echo clearing target directory
mkdir -p $TARGET
if [ -n "$TARGET" ]; then
  rm -rf $TARGET/*
fi
PATH=$TARGET/bin:$PATH
export LD_RUN_PATH=$LD_RUN_PATH:$TARGET/lib

# do not delete download directory contents, as this is intended to prevent
# unnecessary downloads
mkdir -p $DOWNLOAD

# create scratch directory
echo clearing scratch directory
mkdir -p $SCRATCH
if [ -n "$SCRATCH" ]; then
  rm -rf $SCRATCH/*
fi


# install GCC, if requested
if [ -n "$GCC_SRC_PATH" ]; then
  . $INSTALLHOME/installGCC.sh
  echo GCC installation status: $? | tee -a $LOGFILE
  sleep 20
else
  CXX=`which g++`
  CC=`which gcc`
  FC=`which gfortran`
  FTNLIB="-lgfortran -lpthread"
fi

export CXX
export CC
echo CXX = $CXX       >> $TARGET/Makefile.Local
echo FTNLIB = $FTNLIB >> $TARGET/Makefile.Local

# make BLAS & LAPACK available
echo BLASLIB = $LAPACKLIB $BLASLIB >> $TARGET/Makefile.Local

# install required packages
# not called in subshells such that packages can set environment variables
echo Installing libraries >> $LOGFILE
. $INSTALLHOME/installBoost.sh
. $INSTALLHOME/installSuiteSparse.sh
. $INSTALLHOME/installMUMPS.sh
. $INSTALLHOME/installITSOL.sh
. $INSTALLHOME/installLibamiramesh.sh
. $INSTALLHOME/installHypre.sh
. $INSTALLHOME/installDUNE.sh
. $INSTALLHOME/installFP.sh

echo LD_RUN_PATH = $TARGET/lib $TARGET/lib64 $BLASLIBDIR >> $TARGET/Makefile.Local


echo ==========================================================================
echo For executing programs you\'ll probably want to set the LD_LIBRARY_PATH:
echo LD_LIBRARY_PATH=$TARGET/lib:$TARGET/lib64:$BLASLIBDIR


# -----------------------------------------------------------------------------
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi




