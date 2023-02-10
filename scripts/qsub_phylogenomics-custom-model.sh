#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=110G,h_rt=720:00:00
#$ -o tmp/
#$ -e tmp/

# input
fas=$1 # fasta
rid=$2 # run id
nth=$3 # num cpus
mod=$4
mod_file=$5
typ=$6 # DNA,BIN,AA,MORPH

if [ -z "$typ" ] ; then typ="AA" ; fi

echo "data type: ${typ}"

iqtree \
	-s ${fas} \
	-m ${mod} \
	-mdef ${mod_file} \
	-nt ${nth} \
	-st ${typ} \
	-ntmax ${nth} \
	-bb 1000 \
	-wsr \
	-wbt \
	-pre ${fas%%.fasta}.iqt.${rid}
