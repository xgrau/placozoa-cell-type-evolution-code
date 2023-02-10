# input
fas=$1 # folder with fasta input
nth=$2 # num threads

path_broc="/home/xavi/Programes/Broccoli/broccoli.py"
python ${path_broc} -threads ${nth} -dir ${fas} -phylogenies ml -kmer_size 10000
