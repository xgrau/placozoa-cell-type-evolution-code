#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q short-sl7,mem_512_12h
#$ -l virtual_free=20G,h_rt=12:00:00
#$ -o tmp/
#$ -e tmp/

# input
sid=$1 # sps id
fas=$(readlink --canonicalize $2) # genome fasta file (can be gzipped)
gtf=$(readlink --canonicalize $3) # GTF
fql=$(readlink --canonicalize $4) # file with a list of paths to cell ranger-formatted fastq files
out=$(readlink --canonicalize $5) # output folder for reference genome
nth=$6

# retain
timestamp=$(date +%s)
oui=$(pwd)

# ucsc tools path
ucsc_path="/users/asebe/xgraubove/Programes/ucsc_tools/"
macs_path="/users/asebe/xgraubove/miniconda3/envs/atacpip/bin/macs2"
star_path=$(which STAR)

# create dirs
mkdir -p ${out}/index
mkdir -p ${out}/mapping
mkdir -p ${out}/macs2

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

# prepare reference genome
echo "star | mkref ${sid} start"
prepare_fasta ${fas} ${out}/tmp.${sid}.${timestamp}.fasta
cp            ${gtf} ${out}/tmp.${sid}.${timestamp}.gtf
echo "star | mkref ${sid} done"
${star_path} --runThreadN 6 \
	--runMode genomeGenerate \
	--genomeDir ${out}/index/${sid} \
	--genomeFastaFiles ${out}/tmp.${sid}.${timestamp}.fasta \
	--sjdbGTFfile ${out}/tmp.${sid}.${timestamp}.gtf \
	--genomeSAindexNbases 11 \
	--sjdbOverhang 60 \
	--limitSjdbInsertNsj 2022663 \
	&& rm ${out}/tmp.${sid}.${timestamp}.fasta ${out}/tmp.${sid}.${timestamp}.gtf

# run STAR with all 2nd reads
r2s=$(while read -a rip ; do rip=$(readlink --canonicalize ${rip}) ; ls -d $rip/* | grep "_R2_" | sort -u ; done < <(grep "^${sid}" ${fql} | cut -f 4 | sort -u) | xargs | tr ' ' ',' )
echo ${r2s}
# r2s="/users/asebe/xgraubove/prova_h2.r2.fq.gz"

echo "star | run ${sid} start"
${star_path} --genomeDir ${out}/index/${sid} \
	--runThreadN ${nth} \
	--readFilesIn ${r2s} \
	--readFilesCommand zcat \
	--outFileNamePrefix ${out}/mapping/pool_${sid}. \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 5 \
	--limitBAMsortRAM 50000000000 \
	--outFilterMismatchNmax 3 \
	--alignIntronMax 5500 \
	--genomeLoad LoadAndRemove \
	--readNameSeparator ' ' \
	--outSAMattributes Standard
echo "star | run ${sid} done"

echo "star | get strand-specific alignments"
sambamba view -f bam --filter "strand=='+'" -t ${nth} -o ${out}/mapping/pool_${sid}.plus.bam  ${out}/mapping/pool_${sid}.Aligned.sortedByCoord.out.bam
sambamba view -f bam --filter "strand=='-'" -t ${nth} -o ${out}/mapping/pool_${sid}.minus.bam ${out}/mapping/pool_${sid}.Aligned.sortedByCoord.out.bam


# run MACS2 to find negative and positive peaks
echo "macs | calculate effective genome size ${fas}"
eff_genome_size=$(awk '{ s=s+$1 } END { print int(s*0.9) }' ${out}/index/${sid}/chrLength.txt)

echo "macs | peak calling on plus"
${macs_path} callpeak \
	-t ${out}/mapping/pool_${sid}.plus.bam \
	-f BAM \
	-g ${eff_genome_size} \
	--keep-dup 20 \
	-q 0.01 \
	--shift 1 \
	--extsize 20 \
	--broad \
	--bdg \
	--nomodel \
	--min-length 30 \
	-n pool_${sid}_plus \
	--outdir ${out}/macs2
	
echo "macs | peak calling on minus"
${macs_path} callpeak \
	-t ${out}/mapping/pool_${sid}.minus.bam \
	-f BAM \
	-g ${eff_genome_size} \
	--keep-dup 20 \
	-q 0.01 \
	--shift 1 \
	--extsize 20 \
	--broad \
	--bdg \
	--nomodel \
	--min-length 30 \
	-n pool_${sid}_minus \
	--outdir ${out}/macs2

# create concatenated bed files with 3' peaks
cat \
	<(awk '{ OFS="\t" } { $6 = "+" ; print $0 }' ${out}/macs2/pool_${sid}_plus_peaks.broadPeak)  \
	<(awk '{ OFS="\t" } { $6 = "-" ; print $0 }' ${out}/macs2/pool_${sid}_minus_peaks.broadPeak)  \
	| sort -k 1,1 -k2,3n > ${out}/pool_${sid}_peaks_3p.broadPeak
	
# finish
echo "macs | $(wc -l ${out}/pool_${sid}_peaks_3p.bed | awk '{print $1 }') peaks found"
echo "macs | all done!"
