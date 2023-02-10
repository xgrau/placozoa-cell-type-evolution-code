#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q short-sl7
#$ -l virtual_free=20G,h_rt=4:00:00
#$ -o tmp/
#$ -e tmp/

# input
sid=$1 # sps id
fas=$(readlink --canonicalize $2) # genome fasta file (can be gzipped)
gtf=$(readlink --canonicalize $3) # GTF
our=$(readlink --canonicalize $4) # output folder for reference genome
nth=$5

# retain
timestamp=$(date +%s)
oui=$(pwd)

# path to cellranger
cellranger_path="/users/asebe/xgraubove/Programes/cellranger-6.1.1/cellranger"

# create dirs
mkdir -p ${our}

# functions
# unpack fasta if necessary
prepare_fasta () {
    local fas=$1
	local out=$2
    if [[ $fas == *.gz ]] ; then
        echo "# unpack ${fas} gz to ${out}"
        zcat ${fas} > ${out}
    else
        echo "# unpack ${fas} cat to ${out}"
		cp ${fas} ${out}
    fi
}

# add exon features if necessary
prepare_genes_gtf () {
    local gtf=$1
	local out=$2
	# check if features are indicated as genes, transcripts, etc
	if awk '$3 == "gene"' ${gtf} | grep -q -w "gene" ; then
		awk 'BEGIN { OFS="\t" ; FS="\t" } { if ($3 == "gene") { $3="exon" ; print $0 } }' ${gtf} > ${out}
	elif awk '$3 == "transcript"' ${gtf} | grep -q -w "transcript" ; then
		awk 'BEGIN { OFS="\t" ; FS="\t" } { if ($3 == "transcript") { $3="exon" ; print $0 } }' ${gtf} > ${out}
	fi
}


# prepare reference genome
echo "cellranger | mkref ${sid} start"
prepare_fasta     ${fas} ${our}/tmp.${sid}.${timestamp}.fasta
prepare_genes_gtf ${gtf} ${our}/tmp.${sid}.${timestamp}.gtf
cd ${our}
${cellranger_path} \
	mkref \
	--nthreads=${nth} \
	--genome=gdb_${sid} \
	--fasta=tmp.${sid}.${timestamp}.fasta \
	--genes=tmp.${sid}.${timestamp}.gtf \
	&& rm tmp.${sid}.${timestamp}.fasta tmp.${sid}.${timestamp}.gtf
cd ${oui}
echo "cellranger | mkref ${sid} done"

# finish
echo "cellranger | all ${sid} all done"
