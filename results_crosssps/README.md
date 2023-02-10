# Cross-species comparison of cell types

Align and compare cell types in placozoans.

## Pairwise cell type comparisons based on the iterative comparison of coexpression (ICC) procedure

We select gene markers for cross-species comparisons based on the [ICC procedure](https://link.springer.com/article/10.1186/gb-2007-8-4-r50), which we then use to evaluate cross-species cell type similarity based on marker gene expression correlation.

1. Select cross-species markers with ICC, based on `it4` metacells. This produces pairs of homologs for each species pairs (folder `results_alignment_icc`), with their expression conservation value (**EC score**). These pairs of homologs include all one-to-one orthologs plus the best one-to-one pairs selected within each group of many-to-many paralogs for that species pair.

```bash
Rscript s01_icc_marker_selection_2021-08-27.R
```

2. Create cross-species `R` objects for each species pair (`results_metacell_it4`), using ICC-selected markers, and storing quantile-normalised fold-change expression matrices for metacells (`mcs` files), cell types (`cts` files), and broad cell types (`bct` files):

```bash
Rscript s02_icc_celltype_csps_it4_2022-01-28.R
```

3. Compare cell type-cell type similarity for each species pair, using EC score-weighted Pearson correlation. This analysis is run for cell types, broad cell types, and metacells:

```bash
Rscript Rscript s03_icc_expression_fps_it4_2022-07-06.R
```

4. We can use the same ICC-base procedure with the extended dataset with other metazoans (3 cnidarians, 2 bilaterians, *Spongilla*, and *Mnemiopsis*):

```bash
Rscript s11_icc_marker_selection_transphyla_2022-07-20b.R
Rscript s12_icc_celltype_csps_transphyla_it4_2022-07-20.R
Rscript s13_icc_expression_fps_it4_transphyla_2022-07-20.R
Rscript s16_pairwise_similarities_transphyla_2023-01-03.R
```

## Run cell type trees

1. Create cell type trees and matched expression heatmaps for placozoans:

```bash
Rscript s04_cell_type_trees_2022-08-24.R
Rscript s05_expression_markers_csps_2022-08-25.R
```

2. Same, including cnidarians (from the additional dataset with extra metazoans):

```bash
Rscript s14_cell_type_trees_transphyla_2022-08-24.R
Rscript s15_expression_markers_transphyla_2022-08-24.R
```

## Pairwise cell type comparisons using SAMap

1. Create environment:

```bash
# samap 0.3.3 is required to work with more than two species at the same time
conda create -n samap2 -c bioconda -c main -c conda-forge samap=0.3.3
```

2. Get blast database:

```bash
conda activate samap

# blast dbs
sps_list="Tadh Hhon TrH2 HoiH23 Nvec Hvul Spis"
for spi in ${sps_list} ; do
  makeblastdb -dbtype prot -parse_seqids -in ../data/reference/${spi}_long.pep.fasta
done

# blast alignments
for spi in ${sps_list} ; do
for spj in ${sps_list} ; do
if [ ${spi} != ${spj} ] ; then
qsub -N bla-${spi}_${spj} -pe smp 8 qsub_blast.sh ../data/reference/${spi}_long.pep.fasta ../data/reference/${spj}_long.pep.fasta data_samap/blast.${spi}-${spj}.csv
fi
done
done
gzip -v data_samap/blast.*.csv

# rename transcript to gene names...
Rscript s30_rename_tx2ge_samap.R
# prepare samap folder structure
for p in Tadh-TrH2 Tadh-Hhon Tadh-HoiH23 TrH2-Hhon TrH2-HoiH23 Hhon-HoiH23 ; do
s1=$(echo $p | cut -f1 -d '-')
s2=$(echo $p | cut -f2 -d '-')
mkdir -p data_samap/blastdb/${s1}${s2}/
zcat data_samap/blast.${s1}-${s2}.genes.csv.gz | awk '$11 < 1e-18' > data_samap/blastdb/${s1}${s2}/${s2}_to_${s1}.txt
zcat data_samap/blast.${s2}-${s1}.genes.csv.gz | awk '$11 < 1e-18' > data_samap/blastdb/${s1}${s2}/${s1}_to_${s2}.txt
done
```

3. Prepare single-species SAM objects:

```bash
conda activate samap
python s31_samap_prepare_SAM_2021-08-25.py
```

4. Run pairwise and 4-sps SAMap analyses:

```bash
conda activate samap
python s32_samap_objects_2021-09-07.py
```

5. Downstream analyses (pairwise heatmaps, alluvial plots to match cell types, etc.):

```bash
python s33_samap_postalignment_2021-09-07.py
python s34_samap_umap_plot_2021-08-24.R
Rscript s35_samap_heatmap_plot_2021-09-07.R
```
