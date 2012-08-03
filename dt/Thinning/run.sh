#!/bin/sh

file=$1

root=`echo $file | sed 's/\./ /g' | awk '{print $1}'`
# size=`echo $file | sed 's/\./ /g' | awk '{print $2}'`
# sx=`echo $size | sed 's/x/ /g' | awk '{print $1}'`
# sy=`echo $size | sed 's/x/ /g' | awk '{print $2}'`
# sz=`echo $size | sed 's/x/ /g' | awk '{print $3}'`

echo "$sx $sy $sz"

skelFile=${root}.skel

./thin $file ${root}.skel
