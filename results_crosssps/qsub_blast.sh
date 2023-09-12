#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7,mem_512_12h,short-sl7
#$ -l virtual_free=20G,h_rt=6:00:00
#$ -o tmp/
#$ -e tmp/

ref=$1 # blast db
que=$2 # query fasta
out=$3 # output
bin=$4 # binary

blastp -db ${ref} -query ${que} -outfmt 6 -out ${out} -num_threads 30
