#!/bin/sh

rm rezultat.csv
rm SimulacijaCestica_BezReda*

nodeCount=`python readJsonValue.py conf.json nodeCount`
ppn=`python readJsonValue.py conf.json processorsPerNode`

echo "Node count :"
echo $nodeCount
echo "PPN:"
echo $ppn

qsub run.sub -l nodes=$nodeCount:ppn=$ppn