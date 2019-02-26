#!/bin/sh

set -eo pipefail

R1=$1
R2=$2

size=$(stat -c%s $R1)
dir=`dirname $R1`
newR1=`basename $R1`
newR2=`basename $R2`
if [ $size -gt 7000000000 ]; then
	module load seqtk/1.2-r102
	seqtk sample -s100 $R1 70000000 |gzip  >$dir/subsample_${newR1} 
	seqtk sample -s100 $R2 70000000 |gzip  >$dir/subsample_${newR2} 
	mv -f $dir/subsample_${newR1} $R1
	mv -f $dir/subsample_${newR2} $R2
fi

