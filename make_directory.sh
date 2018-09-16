#!/bin/bash

OUTDIR="OUTPUT"
CONFDIR="CONFIG"
EIGENDIR="EIGENS"

if [ ! -d ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

if [ ! -d ${CONFDIR} ]; then
  mkdir ${CONFDIR}
fi

if [ ! -d ${EIGENDIR} ]; then
  mkdir ${EIGENDIR}
fi

