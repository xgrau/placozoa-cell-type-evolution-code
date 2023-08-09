# Phylogenomics of placozoans

## Supermatrix generation and phylogentic analysis

Summary of the steps:

1. Broccoli for OG identification.
2. Select orthogroups as candidate markers, using completeness and clustering coefficient statistics from [***Broccoli***](https://academic.oup.com/mbe/article/37/11/3389/5865275).
3. Select **single orthologs per species** within each orthogroup using [***Possvm***](https://academic.oup.com/mbe/article/38/11/5204/6342420) and branch length-related statistics, and concatenate the markers into a supermatrix.
4. Select **compositionally homogeneous markers** using the [`p4` phylogenetics toolkit](https://p4.nhm.ac.uk/tutorial/tut_compo.html).
5. Select **high-information content gene markers** with [`MARE`](https://p4.nhm.ac.uk/tutorial/tut_compo.html) tree-likeness scores.
6. Apply amino-acid recoding schemes (SR4, SR6, Dayhoff6).

### Steps

Guide to prepare the Metazoa-only (codenamed `meto`) and Metazoa+Choanoflagellata (`metc` datasets):

1. Clean peptides with `cdhit`:

```bash
# input peptides in data/peptides_input
mkdir -p data/peptides_input
mkdir -p data/peptides_clustered
for i in data/peptides_input/*.fasta ; do cdhit -i ${i} -o data/peptides_clustered/$(basename $i) -c 0.99 ; done
cat data/peptides_clustered/*.fasta > concatenated_fasta.fasta
esl-sfetch --index data/concatenated_fasta.fasta
```

2. Run **Broccoli** to find orthogroups for each dataset:

```bash
## Metazoa dataset (63 species), or `meto`
# move input fasta files to results_broccoli_meto/dataset_filtered
mkdir -p results_broccoli_meto/dataset_filtered
mkdir -p results_broccoli_meto/alignments
mkdir -p results_broccoli_meto/alignments_prefiltering
while read i ; do cp data/peptides_input/${i}.fasta results_broccoli_meto/dataset_filtered/${i}.fasta ; done < data/sps_list_meto.txt
# launch broccoli
cd results_broccoli_meto/
bash ../../scripts/qsub_broccoli-ml.sh dataset_filtered 30
cd ..
# select markers present in at least 50% of the dataset, with a clustering coefficient >80%
awk '$2 >= 63*0.5 && $5 >= 0.8' results_broccoli_meto/dir_step3/statistics_per_OG.txt | awk 'NR>1' > results_broccoli_meto/statistics_per_OG.txt

## Metazoa+Choanos dataset (81 species), or `metc`
# move input fasta files to results_broccoli_metc/dataset_filtered
mkdir -p results_broccoli_metc/dataset_filtered
mkdir -p results_broccoli_metc/alignments
mkdir -p results_broccoli_metc/alignments_prefiltering
while read i ; do cp data/peptides_input/${i}.fasta results_broccoli_metc/dataset_filtered/${i}.fasta ; done < data/sps_list_metc.txt
# launch broccoli
cd results_broccoli_metc/
bash ../../scripts/qsub_broccoli-ml.sh dataset_filtered 30
cd ..
# select markers present in at least 75% of the dataset, with a clustering coefficient >90%
awk '$2 > 81*0.75 && $5 > 0.9' results_broccoli_meto/dir_step3/statistics_per_OG.txt | awk 'NR>1' > results_broccoli_meto/statistics_per_OG.txt
```

2. Get **candidate markers** (genes reported as one-to-one orthologs or quasi by `broccoli`, based on taxa occupancy and clustering coefficient):

```bash
# get fasta of markers, `meto`:
while read -a i ; do 
grep -w ${i} results_broccoli_meto/dir_step3/orthologous_groups.txt | cut -f2 | tr ' ' '\n' > results_broccoli_meto/alignments_prefiltering/${i}.txt 
esl-sfetch -f data/concatenated_fasta.fasta results_broccoli_meto/alignments_prefiltering/${i}.txt > results_broccoli_meto/alignments_prefiltering/${i}.fasta
done < results_broccoli_meto/statistics_per_OG.txt

# get fasta of markers, `metc`:
while read -a i ; do 
grep -w ${i} results_broccoli_metc/dir_step3/orthologous_groups.txt | cut -f2 | tr ' ' '\n' > results_broccoli_metc/alignments_prefiltering/${i}.txt 
esl-sfetch -f data/concatenated_fasta.fasta results_broccoli_metc/alignments_prefiltering/${i}.txt > results_broccoli_metc/alignments_prefiltering/${i}.fasta
done < results_broccoli_metc/statistics_per_OG.txt
```

3. Launch initial gene trees, from which we'll select a **single ortholog per species** using [***Possvm***](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8557443/). First, we'll build alignments (`mafft`), trim them (`clipkit`), build gene trees, remove outliers (`iqtree` and `treeshrink`), and rebuild clean alignments if outliers are found. Second, we run `possvm` to refine the groups of orthologs identified by Broccoli, clean each marker group. Finally, we save a set of unaligned proteins with a single ortholog per species (in `results_broccoli_meto/alignments/` and `results_broccoli_metc/alignments/`).

```bash
# for `meto` (metazoa only)
ls -d results_broccoli_meto/alignments_prefiltering/OG_[0-9]*.fasta | grep -v ".lt.fasta" | grep -v ".l.fasta" > results_broccoli_meto/list_prefiltered_alignments.txt
qsub -N phygt-pre-meto -pe smp 2 -t 1-$(grep -c "fasta" results_broccoli_meto/list_prefiltered_alignments.txt) qsub_alignment-array.sh results_broccoli_meto/list_prefiltered_alignments.txt 2

# for `metc` (metazoa + choanos)
ls -d results_broccoli_metc/alignments_prefiltering/OG_[0-9]*.fasta | grep -v ".lt.fasta" | grep -v ".l.fasta" > results_broccoli_metc/list_prefiltered_alignments.metc.txt
qsub -N phygt-pre-metc -pe smp 2 -t 1-$(grep -c "fasta" results_broccoli_metc/list_prefiltered_alignments.metc.txt) qsub_alignment-array.sh results_broccoli_metc/list_prefiltered_alignments.metc.txt 2

# find ortholog groups with possvm
# we will run just two iterations of the iterative midroot procedure, to avoid over-fitting in small trees
mkdir -p results_broccoli_meto/possvm_prefiltering
mkdir -p results_broccoli_metc/possvm_prefiltering
for i in results_broccoli_meto/alignments_prefiltering/*treefile ; do possvm -i ${i} -itermidroot 2 -o results_broccoli_meto/possvm_prefiltering ; done
for i in results_broccoli_metc/alignments_prefiltering/*treefile ; do possvm -i ${i} -itermidroot 2 -o results_broccoli_metc/possvm_prefiltering ; done

# select largest possvm-refined OG per gene tree, drop redundant species-specific paralogs and long-branch sequences
# this will save unaligned proteins to results_broccoli_metc/alignments/ and results_broccoli_meto/alignments/
Rscript s01_parse_possvm_2022-03-21.R
```

4. Relaunch gene trees with the single-ortholog alignments. This will produce accurate single-gene alignments that can be concatenated into a supermatrix. In parallel, it'll also test for the compositional homogeneity of each marker gene using the `p4` phylogenetics toolkit.

```bash
# for `meto` (metazoa only)
ls -d results_broccoli_meto/alignments/OG_[0-9]*.fasta | grep -v ".lt.fasta" | grep -v ".l.fasta" | grep -v ".ltt.fasta" > results_broccoli_meto/list_filtered_alignments.txt
qsub -N phygt-fil-meto -pe smp 2 -t 1-$(grep -c "fasta" results_broccoli_meto/list_filtered_alignments.txt) qsub_alignment-array.sh results_broccoli_meto/list_filtered_alignments.txt 2

# for `metc` (metazoa + choanos)
ls -d results_broccoli_metc/alignments/OG_[0-9]*.fasta > results_broccoli_metc/list_filtered_alignments.metc.txt
qsub -N phygt-fil-metc -pe smp 4 -t 1-$(grep -c "fasta" results_broccoli_metc/list_filtered_alignments.metc.txt) qsub_alignment-array.sh results_broccoli_metc/list_filtered_alignments.metc.txt 4
```

5. Build full **concatenated alignment**, and alignment with markers that pass the **compositional homogeneity tests** from `p4`:

```bash
Rscript s02_filter_and_concatenate_2022-10-14.R
```

6. Use `MARE` to **select high-information content markers** based on tree-likeness scores, and trim down the concatenated alignments:

```bash
# `meto`, metazoa only 63 sps
ali=meto10.all.fasta
cha=meto10.all.charset.txt
mad=2
cd results_broccoli_meto/
MARE ${cha} ${ali} -m -t 100 -d ${mad}
mv results/ mare_reduction_$(basename ${ali%%.fasta})_d${mad}
cp mare_reduction_$(basename ${ali%%.fasta})_d${mad}/${ali}_reduced ${ali%%.fasta}_mare_d${mad}.fasta
cd ..

# `metc` dataset, metazoa+choanos
ali=mat.cho09.all.fasta
cha=mat.cho09.all.charset.txt
mad=2
cd results_broccoli_metc/
MARE ${cha} ${ali} -m -t 100 -d ${mad}
mv results/ mare_reduction_$(basename ${ali%%.fasta})_d${mad}
cp mare_reduction_$(basename ${ali%%.fasta})_d${mad}/${ali}_reduced ${ali%%.fasta}_mare_d${mad}.fasta
cd ..
```

7. **Aminoacid recoding** (SR4, SR6, Dayhoff6) of the concatenated alignments:

```bash
# metc
for fas in results_broccoli_meto/meto10.*.fasta ; do
  Rscript s03_recoding_alignments.R ${fas}
done

# meto
for fas in results_broccoli_metc/mat.cho09.*.fasta ; do
  Rscript s12_recoding_alignments.R ${fas}
done
```

9. Launch **species trees with maximum-likelihood**, using IQ-TREE:

```bash
# unrecoded datasets, AA:
fas="results_broccoli_metc/mat.cho09.comhom.fasta"
fas="results_broccoli_metc/mat.cho09.all_mare_d2.fasta"
fas="results_broccoli_meto/meto10.all_mare_d2.fasta"
fas="results_broccoli_meto/meto10.comhom.fasta"
fai=$(basename ${fas%%.fasta})
rid="LGC60"
bash ../scripts/qsub_phylogenomics-model.sh ${fas} ${rid} 12 LG+F+G+C60

# aminoacid recoding schemes, SR4
fas="results_broccoli_metc/mat.cho09.comhom.recSR4.fasta"
fas="results_broccoli_metc/mat.cho09.all_mare_d2.recSR4.fasta"
fas="results_broccoli_meto/meto10.all_mare_d2.recSR4.fasta"
fas="results_broccoli_meto/meto10.comhom.recSR4.fasta"
fai=$(basename ${fas%%.fasta})
rid="GTRC60"
bash ../scripts/qsub_phylogenomics-custom-model.sh ${fas} ${rid} 12 xmC60SR4 data/xmC60SR4.nex DNA

# aminoacid recoding schemes, SR6
fas="results_broccoli_metc/mat.cho09.comhom.recSR6.fasta"
fas="results_broccoli_metc/mat.cho09.all_mare_d2.recSR6.fasta"
fas="results_broccoli_meto/meto10.all_mare_d2.recSR6.fasta"
fas="results_broccoli_meto/meto10.comhom.recSR6.fasta"
fai=$(basename ${fas%%.fasta})
rid="GTRC60"
bash ../scripts/qsub_phylogenomics-custom-model.sh ${fas} ${rid} 12 xmC60SR6 data/xmC60SR6.nex MORPH

# aminoacid recoding schemes, Day6
fas="results_broccoli_metc/mat.cho09.comhom.recDayhoff6.fasta"
fas="results_broccoli_metc/mat.cho09.all_mare_d2.recDayhoff6.fasta"
fas="results_broccoli_meto/meto10.all_mare_d2.recDayhoff6.fasta"
fas="results_broccoli_meto/meto10.comhom.recDayhoff6.fasta"
fai=$(basename ${fas%%.fasta})
rid="GTRC60"
bash ../scripts/qsub_phylogenomics-custom-model.sh ${fas} ${rid} 12 xmC60Dayhoff6 data/xmC60Dayhoff6.nex MORPH
```

10. Launch **species trees with maximum-likelihood**, using Phylobayes MPI:

```bash
# for each unrecoded dataset, launch 4 chains as follows:
fas="results_broccoli_metc/mat.cho09.comhom.phy"
rid="cho09.comhom"
mpirun -np 32 pb_mpi -d ${fas} -catfix C60 -gtr ${rid}

# sequence diversity and compositional homogeneity tests, to run for each chain:
burnin=2000
mpirun -np 32 readpb_mpi -x ${burnin} -div  <chain name>
mpirun -np 32 readpb_mpi -x ${burnin} -comp <chain name>

# test convergence of chains
bpcomp -x ${burnin} <chain 1> <chain 2>
```

11. Summarise trees:

```bash
mkdir -p results_trees
Rscript s04_paint_trees.R
```

## Supertree approach with with ASTRAL

1. Use Astral (5.7.8) to build species tree from single-gene-per-species gene trees:

```bash
# clean concatenation of gene trees
cat results_broccoli_metc/alignments/OG_*treefile | sed "s/_[^:]*:/:/g" > results_broccoli_metc/trees_filtered.astral.metc.concat.all.newick
cat results_broccoli_meto/alignments/OG_*treefile | sed "s/_[^:]*:/:/g" > results_broccoli_meto/trees_filtered.astral.meto.concat.all.newick
```

2. Use [weighted Astral in hybrid mode](https://academic.oup.com/mbe/article/39/12/msac215/6750035) (1.15.2.3) to build species tree from single-gene-per-species gene trees:

```bash
# run Astral in hybrid weighted mode 
astral-hybrid -t 20 -R -u 2 -i results_broccoli_metc/trees_filtered.astral.metc.concat.all.newick -o results_broccoli_metc/trees_filtered.astral.metc.hybrid_output.all.newick
astral-hybrid -t 20 -R -u 2 -i results_broccoli_meto/trees_filtered.astral.meto.concat.all.newick -o results_broccoli_meto/trees_filtered.astral.meto.hybrid_output.all.newick

3. Use [paralog-aware Astral](https://academic.oup.com/mbe/article/37/11/3292/5850411?login=true) (1.15.1.3) to build the species tree with paralog-ridden gene trees:

```bash
# concatenation of gene trees
ls results_broccoli_metc/alignments/OG_*treefile | sed "s/alignments/alignments_prefiltering/" | xargs cat > results_broccoli_metc/trees_prefiltered.astral.metc.concat.all.newick
ls results_broccoli_meto/alignments/OG_*treefile | sed "s/alignments/alignments_prefiltering/" | xargs cat > results_broccoli_meto/trees_prefiltered.astral.meto.concat.all.newick

# create gene to species mapping files
paste <(cat results_broccoli_metc/trees_prefiltered.astral.metc.concat.all.newick | tr ',' '\n' | tr -d '()' | sed "s/:.*//") <(cat results_broccoli_metc/trees_prefiltered.astral.metc.concat.all.newick | tr ',' '\n' | tr -d '()' | sed "s/:.*//" | cut -f1 -d '_') > results_broccoli_metc/trees_prefiltered.astral.metc.concat.all.mapping.csv
paste <(cat results_broccoli_meto/trees_prefiltered.astral.meto.concat.all.newick | tr ',' '\n' | tr -d '()' | sed "s/:.*//") <(cat results_broccoli_meto/trees_prefiltered.astral.meto.concat.all.newick | tr ',' '\n' | tr -d '()' | sed "s/:.*//" | cut -f1 -d '_') > results_broccoli_meto/trees_prefiltered.astral.meto.concat.all.mapping.csv

# run astral in paralog-aware mode
astral-pro -R -t 20 -u 2 -a results_broccoli_metc/trees_prefiltered.astral.metc.concat.all.mapping.csv -i results_broccoli_metc/trees_prefiltered.astral.metc.concat.all.newick -o results_broccoli_metc/trees_prefiltered.astral.metc.pro_output.newick
astral-pro -R -t 20 -u 2 -a results_broccoli_meto/trees_prefiltered.astral.meto.concat.all.mapping.csv -i results_broccoli_meto/trees_prefiltered.astral.meto.concat.all.newick -o results_broccoli_meto/trees_prefiltered.astral.meto.pro_output.newick
```

## Macrosynteny blocks

1. Align proteins all-to-all, all species of interest:

```bash
# blast all species
mkdir -p results_macrosynteny_all/alignments
mkdir -p results_macrosynteny_all/bins
while read i ; do 
while read j ; do 
  if [ $i != $j -a ! -f results_macrosynteny_all/alignments/dmd.${i}-${j}.csv ] ; then 
    echo $i $j
    diamond blastp -d ../data/reference/${j}_long.pep.fasta -q ../data/reference/${i}_long.pep.fasta -o results_macrosynteny_all/alignments/dmd.${i}-${j}.csv --evalue 1e-9 --more-sensitive --max-target-seqs 10
  fi
done < data/species_list_synteny_blocks.long.txt
done < data/species_list_synteny_blocks.long.txt
```

3. Get clusters:

```bash
# with MCL
while read i ; do 
  while read j ; do 
    if [ $i != $j ] ; then 
      awk '{print $1,$2 }' results_macrosynteny_all/alignments/dmd.${i}-${j}.csv 
    fi
  done < data/species_list_synteny_blocks.long.txt
done < data/species_list_synteny_blocks.long.txt > results_macrosynteny_all/mcl.all.abc.csv
mcl results_macrosynteny_all/mcl.all.abc.csv --abc -I 2.1 -o results_macrosynteny_all/mcl.all.out.csv
awk '{ { for(i = 1; i <= NF; i++) { printf("HG%06d\t%s\n", NR,$i) } } }' results_macrosynteny_all/mcl.all.out.csv > results_macrosynteny_all/mcl.all.out.txt
```

3. Match orthologroups to running windows of genes along each chromosome:

```bash
Rscript s60_define_regions_2022-05-30-whole.R
```

4. Define ancestral linkage groups based on homology to three reference species:

```bash
Rscript s61_find_ALG_2022-06-03.R
```

5. Score ancestral linkage groups in each species:

```bash
Rscript s62_score_ALG_2022-06-03.R
```

6. Match new ALGs with Simakov 2022 based on overlaps with *Ephydatia* (only species in the dataset for which there's clear gene IDs).

```bash
Rscript s64_match_Emue_to_Simakov22.R
```
