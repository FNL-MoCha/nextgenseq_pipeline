#!/bin/sh

name=`basename $1 .cns`

echo -e "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean"
grep -v ^chromosome $1 |awk -v name="$name" '{OFS="\t"}{print name,$1,$2,$3,$7,$5}'
