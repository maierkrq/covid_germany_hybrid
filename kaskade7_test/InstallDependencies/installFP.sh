#! /bin/sh
# Copyright 2018 - 2018 Zuse Institute Berlin

cd $SCRATCH
git clone $FP_SRC_PATH

FP_NAME=`basename -s .git $FP_SRC_PATH`

cd $FP_NAME

mkdir -p $TARGET/include/$FP_NAME
cp -r include/* $TARGET/include/$FP_NAME/

echo FPINC = -I$TARGET/include/$FP_NAME    >> $TARGET/Makefile.Local


#------------------------------------------------------------------------------
# Clean up
# if [ -n "$SCRATCH" ]; then
#   echo running rm -r "$SCRATCH"/* in a few seconds
#   sleep 20
#   rm -r "$SCRATCH"/*
# fi
# 

