#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q mem_512_12h,short-sl7
#$ -l virtual_free=10G,h_rt=6:00:00
#$ -o tmp/
#$ -e tmp/
#$ -t 1-3000

# input
l=$1 # list of input alignments
i=$(awk 'NR == '"${SGE_TASK_ID}"'' ${l})
c=$2 # NUM CPUS
path_iqtree="/users/asebe/xgraubove/Programes/iqtree-2.1.0-Linux/bin/iqtree2"

echo "input alignment: ${i}"

function do_alitrimphy {

	local i=$1
	local c=$2
	local il=$3
	local ilt=$4
	local tre=$5

	# align
	if [ -s $il ] ; then
		echo "$il already exists, go to next step"
	else
		mafft --genafpair --thread $c --reorder --maxiterate 10000 $i > $il
	fi
	
	# trim
	if [ -s $ilt ] ; then
        echo "$ilt already exists, go to next step"
	else
		clipkit $il -m kpic-gappy -o $ilt -g 0.7
		bioawk -c fastx '{ B=$2 ; O=gsub("-","",$2) } {if (O != length(B))  { print ">"$1"\n"B }}' ${ilt} > ${ilt}.noallgap
		mv ${ilt}.noallgap ${ilt}
	fi
	
	# # probabilistic scores with zorro
	# if [ -s ${il%.l.fasta}.l.zorro.txt  ] ; then
    #     echo "${il%.l.fasta}.l.zorro.txt already exists, go to next step"
	# else 
	# 	zorro ${i} > ${i%%.l.fasta}.lz.txt
	# fi	


	# iqtree
	${path_iqtree} \
		-s $ilt \
		-m TEST \
		-mset LG,WAG,JTT \
		-nt AUTO \
		-ntmax $c \
		-bb 1000 \
		-pre $tre \
		-nm 10000 \
		-nstop 200 \
		-cptime 1800

}

# do alignments and trees
do_alitrimphy $i $c ${i%%.fasta}.l.fasta ${i%%.fasta}.lt.fasta ${i%%.fasta}.iqtree

# now check if there are outlier sequences in the tree
python /users/asebe/xgraubove/Programes/treeshrink/run_treeshrink.py -c -t ${i%%.fasta}.iqtree.treefile -m per-gene -q 0.05 -s 10,1 -f ${i%%.fasta}.treeshrink.firstpass.txt

# if there's an outlier sequence in the tree, we'll have to re-run the tree...
if [ -s ${i%%.fasta}.treeshrink.firstpass.txt ] ; then

	# backup previous files
	for o in ${i%%.fasta}.l.fasta ${i%%.fasta}.lt.fasta ${i%%.fasta}.iqtree.* ; do
		mv $o ${o}.firstpass
	done

	# remove unwanted seqs
	bioawk -c fastx '{ print $1, $2 }' $i | fgrep -v -w -f ${i%%.fasta}.treeshrink.firstpass.txt | awk '{ print ">"$1"\n"$2 }'  > $i.shrunk

	# redo trees
	do_alitrimphy $i.shrunk $c ${i%%.fasta}.l.fasta ${i%%.fasta}.lt.fasta ${i%%.fasta}.iqtree

fi

# test compositional heterogeneity
python ../scripts/compositional_heterogeneity_test_p4.py \
	-a ${i%%.fasta}.lt.fasta \
	-t ${i%%.fasta}.iqtree.treefile \
	-o ${i%%.fasta}.p4.txt \
	-m lg -g 4 -n 1000 -b
