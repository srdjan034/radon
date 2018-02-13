#!/bin/sh

rm rezultat_*.csv
rm SimulacijaCesticaJob*

nodeCount=`python readJsonValue.py conf.json nodeCount`
ppn=`python readJsonValue.py conf.json processorsPerNode`

echo "Node count :"
echo $nodeCount
echo "PPN:"
echo $ppn

qsub run.sub -l nodes=$nodeCount:ppn=$ppn