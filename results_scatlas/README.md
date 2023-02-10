# Cell atlases for placozoans

How to **map scRNA-seq data**, demultiplex Clicktagged libraries, and produce a  **`metacell` clustering** for each placozoan species.

All through this repository, the following short codes are used to identify each species:

* `Tadh` is *Trichoplax adhaerens* H1
* `TrH2` is *Trichoplax* sp. H2
* `Hhon` is *Hoilungia hongkongensis* H13
* `HoiH23` is *Cladtertia collaboinventa* H23 (it's still called `Hoi.,.` because this project started [before this species was formally named](https://www.frontiersin.org/articles/10.3389/fevo.2022.1016357/full)).

## Input data

The following scRNA-seq libraries have been used for each species. To reproduce the analyses below, create a `data_cdna` folder, download the data from [NCBI](XXXX), and save it in this folder using the following folder structure:

```bash
# columns:
# 1. species name
# 2. sample name
# 3. file (two entries per sample, for R1 and R2); some samples appear more than once because they are Clicktag-multiplexed samples with cells from more than one species. See Clicktag instructions below on how to remove cells from the other species.
# 4. folder where fastq files in field 3 should be stored
Hhon	H13_2_ACME_10x	H13_2_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H13_2_ACME_10x
Hhon	H13_2_ACME_10x	H13_2_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H13_2_ACME_10x
Hhon	H13_4_ACME_10x	H13_4_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H13_4_ACME_10x
Hhon	H13_4_ACME_10x	H13_4_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H13_4_ACME_10x
Hhon	H13_5_ACMEMeOH_10x	H13_5_ACMEMeOH_10x_S1_L001_R1_001.fastq.gz	data_cdna/H13_5_ACMEMeOH_10x
Hhon	H13_5_ACMEMeOH_10x	H13_5_ACMEMeOH_10x_S1_L001_R2_001.fastq.gz	data_cdna/H13_5_ACMEMeOH_10x
Hhon	H1H13_1_ACME_CT_10x	H1H13_1_ACME_CT_10x_S1_L001_R1_001.fastq.gz	data_cdna/H1H13_1_ACME_CT_10x
Hhon	H1H13_1_ACME_CT_10x	H1H13_1_ACME_CT_10x_S1_L001_R2_001.fastq.gz	data_cdna/H1H13_1_ACME_CT_10x
Hhon	Plac01_H1_H13_10XscRNAseq	Plac01_H1_H13_10XscRNAseq_S1_R1_001.fastq.gz	data_cdna/Plac01_H1_H13_10XscRNAseq
Hhon	Plac01_H1_H13_10XscRNAseq	Plac01_H1_H13_10XscRNAseq_S1_R2_001.fastq.gz	data_cdna/Plac01_H1_H13_10XscRNAseq
HoiH23	H23_1_ACME_10x	H23_1_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H23_1_ACME_10x
HoiH23	H23_1_ACME_10x	H23_1_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H23_1_ACME_10x
HoiH23	H23_2_ACME_10x	H23_2_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H23_2_ACME_10x
HoiH23	H23_2_ACME_10x	H23_2_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H23_2_ACME_10x
HoiH23	H23_3_ACME_10x	H23_3_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H23_3_ACME_10x
HoiH23	H23_3_ACME_10x	H23_3_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H23_3_ACME_10x
HoiH23	H2H23_4_ACME_10x	H2H23_4_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H2H23_4_ACME_10x
HoiH23	H2H23_4_ACME_10x	H2H23_4_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H2H23_4_ACME_10x
HoiH23	H1H23_5_ACME_10x	H1H23_5_ACME_10x_S1_R1_001.fastq.gz	data_cdna/H1H23_5_ACME_10x
HoiH23	H1H23_5_ACME_10x	H1H23_5_ACME_10x_S1_R2_001.fastq.gz	data_cdna/H1H23_5_ACME_10x
HoiH23	Plac02_H2_H23_10XscRNAseq	Plac02_H2_H23_10XscRNAseq_S2_R1_001.fastq.gz	data_cdna/Plac02_H2_H23_10XscRNAseq
HoiH23	Plac02_H2_H23_10XscRNAseq	Plac02_H2_H23_10XscRNAseq_S2_R2_001.fastq.gz	data_cdna/Plac02_H2_H23_10XscRNAseq
Tadh	H1_06_ACME_10x	H1_06_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H1_06_ACME_10x
Tadh	H1_06_ACME_10x	H1_06_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H1_06_ACME_10x
Tadh	H1H13_1_ACME_CT_10x	H1H13_1_ACME_CT_10x_S1_L001_I_001.fastq.gz	data_cdna/H1H13_1_ACME_CT_10x
Tadh	H1H13_1_ACME_CT_10x	H1H13_1_ACME_CT_10x_S1_L001_R1_001.fastq.gz	data_cdna/H1H13_1_ACME_CT_10x
Tadh	H1H13_1_ACME_CT_10x	H1H13_1_ACME_CT_10x_S1_L001_R2_001.fastq.gz	data_cdna/H1H13_1_ACME_CT_10x
Tadh	H1H23_5_ACME_10x	H1H23_5_ACME_10x_S1_R1_001.fastq.gz	data_cdna/H1H23_5_ACME_10x
Tadh	H1H23_5_ACME_10x	H1H23_5_ACME_10x_S1_R2_001.fastq.gz	data_cdna/H1H23_5_ACME_10x
Tadh	Plac01_H1_H13_10XscRNAseq	Plac01_H1_H13_10XscRNAseq_S1_R1_001.fastq.gz	data_cdna/Plac01_H1_H13_10XscRNAseq
Tadh	Plac01_H1_H13_10XscRNAseq	Plac01_H1_H13_10XscRNAseq_S1_R2_001.fastq.gz	data_cdna/Plac01_H1_H13_10XscRNAseq
TrH2	H2_1_ACME_10x	H2_1_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H2_1_ACME_10x
TrH2	H2_1_ACME_10x	H2_1_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H2_1_ACME_10x
TrH2	H2_2_ACME_10x	H2_2_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H2_2_ACME_10x
TrH2	H2_2_ACME_10x	H2_2_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H2_2_ACME_10x
TrH2	H2_3_ACME_10x	H2_3_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H2_3_ACME_10x
TrH2	H2_3_ACME_10x	H2_3_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H2_3_ACME_10x
TrH2	H2H23_4_ACME_10x	H2H23_4_ACME_10x_S1_L001_R1_001.fastq.gz	data_cdna/H2H23_4_ACME_10x
TrH2	H2H23_4_ACME_10x	H2H23_4_ACME_10x_S1_L001_R2_001.fastq.gz	data_cdna/H2H23_4_ACME_10x
TrH2	Plac02_H2_H23_10XscRNAseq	Plac02_H2_H23_10XscRNAseq_S2_R1_001.fastq.gz	data_cdna/Plac02_H2_H23_10XscRNAseq
TrH2	Plac02_H2_H23_10XscRNAseq	Plac02_H2_H23_10XscRNAseq_S2_R2_001.fastq.gz	data_cdna/Plac02_H2_H23_10XscRNAseq
```

## Data processing

This are the recipes for the steps undertaken prior to obtaining UMI count matrices per species.

### Reannotation of 3' regions

Extend gene annotations for each placozoan towards the 3' ends of genes, to include the regions that are covered by scRNA-seq libraries.

1. Map expression data (STAR) and call 3' peaks (MACS2):

```bash
for i in Tadh TrH2 Hhon HoiH23 ; do
  bash s01_qsub_peaks3p.sh ${i} ../data/reference/${i}_gDNA.fasta.gz ../data/reference/${i}_long.annot.gtf data/list_cDNA_libraries.tsv results_annotate_3p/ 12
done
# inside s01_qsub_peaks3p.sh, you should hange the path to your USCS tools, MACS and STAR binaries if you wish to run it locally
```

2. Assign peaks to nearby genes (extending known genes if available, otherwise saving peaks as "orphan"):

```bash
for i in Tadh Hhon TrH2 HoiH23; do
  Rscript s02_extend3p.R -g ../data/reference/${i}_long.annot.gtf -p results_annotate_3p/pool_${i}_peaks_3p.broadPeak -o results_annotate_3p/reannotate_${i}_genes.gtf -m 5000 -a -r ${i} -q 1e-3
done
```

### Mapping scRNA-seq with Cell Ranger

Mapping and UMI counts with [Cell Ranger 6.1.1](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).

1. Create list of expression datasets for each species:

```bash
for i in Tadh TrH2 Hhon HoiH23 ; do
  grep "^$i" data/list_cDNA_libraries.tsv | cut -f4 | sort -u > data/list_cDNA_libraries.cr_${i}.txt
done
```

2. Launch `CellRanger` for each dataset:

```bash
# launch cell ranger mapping for each species, using reannotated genomes
# these files are not provided in this repository (too big); but the processed UMI matrices are available below
for i in Tadh TrH2 Hhon HoiH23 ; do
  bash s03_cellranger_launcher_2021-11-11.sh ${i} ../data/reference/${i}_gDNA.fasta.gz  results_annotate_3p/reannotate_${i}_genes.gtf data/list_cDNA_libraries.cr_${i}.txt align_reannotated/ align_reannotated/ 12
done
# this bash script executes two supplementary scripts as qsub jobs (qsub_cellranger-db.sh and qsub_cellranger-run.sh); but it can be changed to run them sequentially as bash scripts
# inside s03_cellranger_launcher_2021-11-11.sh, you should hange the path to your CellRanger binary if you wish to run it locally
```

3. Create indexes so that `metacell` can import alignments for each species:

```bash
for i in Tadh TrH2 Hhon HoiH23  ; do 
  echo "Batch.Set.ID mat_fn genes_fn cells_fn" > data/list_cDNA_libraries.mcr_${i}.txt
  f=$(ls -d align_reannotated/map_${i}*10kc/outs/filtered_feature_bc_matrix/) 
  for j in $f ; do 
    echo $(echo ${j} | cut -f2 -d '/' | sed "s/^map_${i}_//") ${j}/matrix.mtx.gz ${j}/features.tsv.gz ${j}/barcodes.tsv.gz 
  done >> data/list_cDNA_libraries.mcr_${i}.txt
done
```

4. Save cell UMI count matrices as `metacell` objects, for downstream processing:

```R
library("metacell")
for (spi in c("Tadh","TrH2","Hhon","HoiH23")) {
  mat = mcell_read_multi_scmat_10x(datasets_table_fn = sprintf("data/list_cDNA_libraries.mcr_%s.txt", spi), base_dir = ".")
  scdb_init(base_dir = "data/scdb", force_reinit=TRUE)
  scdb_add_mat(sprintf("scdr_%s", spi), mat)
}
```

### Mapping and demultiplexing Clicktag libraries

Demultiplexing **Clicktag** libraries to obtain a list of *bona fide* cell barcodes that can be assigned to each cell.

1. Environment:

```bash
# create environment
conda create -n clicktag python=3.7
conda activate clicktag
conda install kallisto bustools
conda install -c conda-forge leidenalg scanpy
```

2. Map Clicktag reads to barcodes, correct mapping and assign reads to cells, with `kallisto` and `bustools`:

```bash
bash ../scripts/clicktag-mapping.sh data_clicktag/barcodes_clicktag_H1H13.csv data_clicktag/run_H1H13/H1H13_CT_2runs_R1.fastq.gz data_clicktag/run_H1H13/H1H13_CT_2runs_R2.fastq.gz <(zcat data_clicktag/barcodes_cells_10x_v3.txt.gz) 10xv3 data_clicktag/run_H1H13/ data_clicktag/map_H1H13/ 4
bash ../scripts/clicktag-mapping.sh data_clicktag/barcodes_clicktag_H1H23.csv data_clicktag/run_H1H23/clean_clicktags_H1H23_R1.fastq.gz data_clicktag/run_H1H23/clean_clicktags_H1H23_R2.fastq.gz <(zcat data_clicktag/barcodes_cells_10x_v3.txt.gz) 10xv3 data_clicktag/run_H1H23/ data_clicktag/map_H1H23/ 10
bash ../scripts/clicktag-mapping.sh data_clicktag/barcodes_clicktag_H2H23.csv data_clicktag/run_H2H23/clean_clicktags_H2H23_R1.fastq.gz data_clicktag/run_H2H23/clean_clicktags_H2H23_R2.fastq.gz <(zcat data_clicktag/barcodes_cells_10x_v3.txt.gz) 10xv3 data_clicktag/run_H2H23/ data_clicktag/map_H2H23/ 10
bash ../scripts/clicktag-mapping.sh data_clicktag/barcodes_clicktag_Plac01H1H13.csv data_clicktag/run_Plac01H1H13/clean_clicktags_H1H13_R1.fastq.gz data_clicktag/run_Plac01H1H13/clean_clicktags_H1H13_R2.fastq.gz <(zcat data_clicktag/barcodes_cells_10x_v3.txt.gz) 10xv3 data_clicktag/run_Plac01H1H13/ data_clicktag/map_Plac01H1H13/ 10
bash ../scripts/clicktag-mapping.sh data_clicktag/barcodes_clicktag_Plac02H2H23.csv data_clicktag/run_Plac02H2H23/clean_clicktags_H2H23_R1.fastq.gz data_clicktag/run_Plac02H2H23/clean_clicktags_H2H23_R2.fastq.gz <(zcat data_clicktag/barcodes_cells_10x_v3.txt.gz) 10xv3 data_clicktag/run_Plac02H2H23/ data_clicktag/map_Plac02H2H23/ 10
```

These tables are imported as `metacell` objects in the metacell clustering step below (`s04_filter_cells_clicktag_samples_2021-11-16.R` script).

## Single cell atlases

### Metacell clustering

Steps for **`metacell`** analyses and annotation transfer from MARS-seq.

1. For the Clicktag samples, find which cells to blacklist (cells from that sample whose Clicktag barcodes denote a different species, putative doublets, etc.):

```bash
Rscript s04_filter_cells_clicktag_samples_2021-11-16.R
```

2. Create an initial iteration of each species' cell atlas (`it0`) after removing blacklisted cells (doublets, etc.) and genes (histones, ribosomal proteins). The second script will produce tentative cell type annotations for each metacell (e.g. table `results_metacell_it0/annotation_mc.Hhon.it0.tsv`), based on similarity to the metacells of the previously published *Trichoplax adhaerens* dataset ([Sebé-Pedrós et al, NEE 2018](https://www.nature.com/articles/s41559-018-0575-6)).

```bash
Rscript s05_scatlas_it0_2021-12-16.R
Rscript s06_postscatlas_it0_2021-12-16.R
```

3. Metacell purging step: remove low-quality metacells from `it0` and obtain blacklists of metacells and single cells that won't be included in subsequent analyses. Criteria used for this:

* Small metacells (low UMI criteria *z*-score based, global and tissue-level)
* Metacells without markers (less than 10 genes with fc>1.5)
* Metacells without TFs (zero TFs with fc>1.5)

```bash
Rscript s07_purge_metacells_it0_2021-12-20.R
```

4. Rerun metacell (`it1`) after low-quality metacell removal. At this point we'll also remove metacells that come almost exclusively from unique Clicktagged samples, and which have first-to-third clicktag count ratios that might suggest they contain doublet cells.

* mean ratio first-to-third ratio `>1`
* `>50%` cells from clicktagged samples
* `>90%` cells within the clicktagged samples come from the same sample
* metacells identical to those above

```bash
Rscript s08_scatlas_it1_2022-01-10.R
Rscript s09_postscatlas_it1_2022-01-10.R
Rscript s10_qualitycontrol_it1_2022-01-19.R

# output
# drop mcs Tadh: 68,69,125,129,130,131,132,133,134
# drop mcs TrH2: none
# drop mcs Hhon: 254
# drop mcs HoiH23: none

# identify cells to drop (manual edition of mc lists in script)
Rscript s11_purge_metacells_it1_2022-01-25.R
```

6. After all quality-control steps are done, we produce an intermediate metacell iteration (`it2`) that contain the final metacells. These will be curated and annotated later (but not changed anymore) to obtain the final iteration (`it4`).

```bash
Rscript s12_scatlas_it2_2022-01-25.R
```

7. After `it2` plots are done, the metacell list is stable. Manual annotation and reordering of metacells is now required: copy the `results_metacell_it2/annotation_mc.Tadh.it2.tsv` to `results_metacell_it4/annotation_mc.Tadh.it4.tsv` and manually edit and reorder metacells (preserving metacell ids), assigning good colors according to your good judgement. This will be the input for the **final `it4` metacell iteration**.

9. Quality control plots, expression maps and final tables for the **final `it4` metacell iteration**, with the reordered metacells:

```bash
# maps and markers
Rscript s13_postscatlas_it4_2022-01-13.R
Rscript s14_plot_markers_it4_2022-03-21.R

# quality control (check if there are any metacells left that fail filters; nothing will happen)
Rscript s15_qualitycontrol_it4_2022-01-26.R

# save cell type-level mat and mc objects
Rscript s16_cell_type_fps_2022-09-02.R
```

10. Analysis of expression profiles of transversal cells (`AUCell` to check if transversal cells have expression signatures intermediate between their respective terminal cell types; and check conservation of their marker genes).

```bash
Rscript s17_trans_cells_profiles_it4_2022-08-24.R
```

11. Reclustering of peptidergic cells into a narrower set of metacells:

```bash
Rscript s20_recluster_peptidergic_it4_2022-07-06.R
Rscript s21_postscatlas_peptidergic_it4_2022-01-25.R
```
