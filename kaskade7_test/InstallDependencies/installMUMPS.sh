#! /bin/sh
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $DOWNLOAD
TARFILE=`basename $MUMPS_SRC_PATH`
curl -C -  -o $TARFILE $MUMPS_SRC_PATH

cd $SCRATCH
tar xzf $DOWNLOAD/$TARFILE

SRCDIR=`basename -s .gz $TARFILE`
SRCDIR=`basename -s .tar $SRCDIR`
SRCDIR=$SCRATCH/`basename -s .tgz $SRCDIR`

cd $SRCDIR
# Create a MUMPS configuration. Templates for different architectures are provided
# in Make.inc/. Right now we don't use Scotch (since not yet installed),
# but this could be improved for performance reasons.
# We also use the sequential version of MUMPS and therefore install the pseudo MPI
# library below. Consider to use the parallel version later.
#
# WARNING: As of MUMPS 5.3.3, the gfortran of GCC 10.2 will report an error, where GCC 8.x
# versions of gfortran emit a warning.  We therefore use the system gfortran here (which is
# presumably older) instead of $FC, which is the just compiled and installed gfortran 10.2 or above.
cat > Makefile.inc <<EOF
LSCOTCHDIR =
ISCOTCH =
LSCOTCH =
LPORDDIR = $SRCDIR/PORD/lib
IPORD = -I$SRCDIR/PORD/include
LPORD = -L\$(LPORDDIR) -lpord
LMETISDIR = $TARGET/lib 
IMETIS = -I$TARGET/include
LMETIS = -L\$(LMETISDIR) -lmetis

ORDERINGSF = -Dpord -Dmetis 
ORDERINGSC = \$(ORDERINGSF)

LORDERINGS = \$(LPORD) \$(LMETIS)
IORDERINGSF =
IORDERINGSC = \$(IPORD) \$(IMETIS)

PLAT    =
LIBEXT  = .a
OUTC    = -o
OUTF    = -o
RM = /bin/rm -f
CC = $CC
FC = /usr/bin/gfortran
FL = /usr/bin/gfortran
AR = ar vr
RANLIB = ranlib
LAPACK = $LAPACKLIB



INCSEQ = -I$SRCDIR/libseq
LIBSEQ  = \$(LAPACK) -L$SRCDIR/libseq -lmpiseq

LIBBLAS = $BLASLIB
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_


#Begin Optimized options
# uncomment -fopenmp in lines below to benefit from OpenMP
OPTF    = -O -fPIC # -fopenmp
OPTL    = -O -fPIC # -fopenmp
OPTC    = -O -fPIC # -fopenmp
#End Optimized options

INCS = \$(INCSEQ)
LIBS = \$(LIBSEQ)
LIBSEQNEEDED = libseqneeded

EOF

# As of MUMPS 5.1.2 - 5.3.3 there's a bug in the Makefiles (missing space)
for i in `find . -name Makefile`; do
  mv $i $i.old
  sed -e 's/$(AR)$@/$(AR) $@/g' < $i.old > $i
done

# Build and install MUMPS and the pseudo MPI stub library.
echo building MUMPS >> $LOGFILE
make
if [ $? -ne 0 ]; then
  echo BUILDING MUMPS FAILED: status $? | tee -a $LOGFILE
  exit
fi

cp lib/*.a $TARGET/lib
cp libseq/libmpiseq.a $TARGET/lib
cp include/*.h $TARGET/include

echo MUMPSINC = -I$TARGET/include/                                            >> $TARGET/Makefile.Local
echo MUMPSLIB = -L$TARGET/lib -ldmumps -lmumps_common -lpord -lmpiseq -lmetis >> $TARGET/Makefile.Local

#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi


