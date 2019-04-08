#!/bin/bash

# This script takes as input the path to the output file of 'ntCard'. It then
# calculates log of the number of slots needed by Squeakr to count k-mers.

F0=`grep -P "^F0\t" $1 | cut -f2- -d$'\t'`
f1=`grep -P "^1\t" $1 | cut -f2- -d$'\t'`
f2=`grep -P "^2\t" $1 | cut -f2- -d$'\t'`
fgt2=$[ $F0 - $f1 - $f2 ]
ns=$[ $f1 + 2 * $f2 + 3 * $fgt2 ]
echo "x=$ns;l2=l(x)/l(2);s=scale;scale=0;l2ru=l2+1-(l2%1);np2=2^(l2ru/1);scale=s;if(x > (0.9*np2)) { l2ru=l2ru+1; np2=2*np2; }; s=scale; scale=0; print l2ru/1; scale=s;" | bc -l
echo
