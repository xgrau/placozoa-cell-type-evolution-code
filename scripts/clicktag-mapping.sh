# input
bct=$(readlink --canonicalize $1) # clicktag barcode csv file
fq1=$(readlink --canonicalize $2) # path to clicktag fq R1 (can be gzipped)
fq2=$(readlink --canonicalize $3) # path to clicktag fq R2 (can be gzipped)
bcc=$(readlink --canonicalize $4) # 10x cell barcodes whitelist (txt file)
che=$5                            # 10x chemistry (10xv2, 10xv3)
out=$(readlink --canonicalize $6) # path to output folder for kallisto big files
oug=$(readlink --canonicalize $7) # path to output folder for clicktag counts
nth=$8                            # num threads for kallisto and bustools

# retain
timestamp=$(date +%s)

# check input files
echo "# Input files"
for fii in ${bct} ${fq1} ${fq2} ${bcc} ; do
	ls -lh ${fii}
done

# create mutated set of barcodes
echo "# Create extended list of clicktag barcodes"
python ../scripts/clicktag-featuremap.py ${bct} --fa ${bct%%.csv}.fasta --t2g ${bct%%.csv}.t2g.tsv
kallisto index -i ${bct%%.csv}.fasta.idx -k 11 ${bct%%.csv}.fasta

# map with kallisto
echo "# Map clicktag reads with kallisto"
kallisto bus -i ${bct%%.csv}.fasta.idx \
  -o ${out} \
  -x ${che} \
  -t ${nth} \
  ${fq1} \
  ${fq2}

# # prepare output dirs
mkdir -p ${oug}
mkdir -p ${out}/tmp/

# correct with bustools
echo "# Correct clicktag mapping with bustools"
bustools correct \
  -w ${bcc} \
  -p ${out}/output.bus | \
bustools sort \
  -T ${out}/tmp/ \
  -t ${nth} -p - | \
bustools count \
  -o ${oug}/genes \
  -g ${bct%%.csv}.t2g.tsv \
  -e ${out}/matrix.ec \
  -t ${out}/transcripts.txt \
  --genecounts -

echo "# Done!"

