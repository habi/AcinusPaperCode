#! /bin/sh

echo "test"
file=$1
thr=0.6

root=`echo $file | sed 's/\./ /g' | awk '{print $1}'`
size=`echo $file | sed 's/\./ /g' | awk '{print $2}'`

skelFile=${root}-th${thr}.skel

./dtSkel $file $thr tmp

if [ -e tmp ]; then
    awk '{print $1, $2, $3, 0.8}' tmp > $skelFile
else
    echo "no found"
fi

rm -f tmp