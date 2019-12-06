#!/bin/bash

# This script takes two input arguments:
# 1. The path to the output file of 'ntCard'
# 2. F0 value output by ntCard. 
# It then calculates log of the number of slots needed by Squeakr to count
# k-mers.

F0=$2
f1=`grep -P "\t1\t" $1 | cut -f3- -d$'\t'`
f2=`grep -P "\t2\t" $1 | cut -f3- -d$'\t'`
fgt2=$[ $F0 - $f1 - $f2 ]
ns=$[ $f1 + 2 * $f2 + 3 * $fgt2 ]
echo "x=$ns;l2=l(x)/l(2);s=scale;scale=0;l2ru=l2+1-(l2%1);np2=2^(l2ru/1);scale=s;if(x > (0.9*np2)) { l2ru=l2ru+1; np2=2*np2; }; s=scale; scale=0; print l2ru/1; scale=s;" | bc -l
echo
