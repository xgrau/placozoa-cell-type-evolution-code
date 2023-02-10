# Gene module analysis

These steps require running the scripts in the `results_scatlas` folder first.

## General placozoan gene modules

Steps:

1. Identify gene modules in each species using `WGCNA`:

```bash
# create modules based on metacell-level coexpression, and annotate them with colors, their constituent TFs, etc.
Rscript s01_gmod_analysis_2022-08-22.R
Rscript s02_gmod_annotation_2022-08-22.R
```

2. Identify multi-species gene modules (termed module "communities") based on shared presence of orthologous genes in single-species gene modules:

```bash
Rscript s03_create_communities_2022-10-17.R
```

3. Annotate multi-species gene communities to specific cell types, and run functional enrichment analyses on each of them. This information will be used to curate gene module-to-cell type annotations (in `results_gmod_it4_ct_annotations/module_communities.curated.csv`)

```bash
Rscript s04_annotate_communities_2022-10-17.R
```

3. Score multi-species gene community activity in species' each cell types, based on the fraction of constituent orthologs from a module that's over-expressed (`fc >= 1.5`) in a particular cell type.

```bash
Rscript s05_community_activity_per_cts_2022-10-27.R
```

## Peptidergic and neural gene modules across metazoans

Identify pan-neural genes in cnidarians and bilaterians, and compare their expression with peptidergic genes in placozoans.

Species list:

```R
c("Nvec","Hvul","Spis", "Dmel","Mmus", "Mlei", "Spolac")
```

Steps:

1. Find genes overexpressed in neural metacells for each species, defined as `fp > 2` in at least 10% of the neural (or peptidergic) metacells.

```bash
Rscript s10_panneural_markers_2022-07-13.R
```

2. Reconstruct presence/absence with Dollo parsimony, and three different species trees:

```bash
# for cnidaria+placozoa tree:
python ../scripts/possvm_reconstruction.py -tree ../data/species_tree_CP.newick -ort results_panneural_markers/matrix.all.long.ogexpression.csv -out results_panneural_markers/anc.all.long.ogexpression.CP
python ../scripts/possvm_reconstruction.py -tree ../data/species_tree_CP.newick -ort results_panneural_markers/matrix.all.long.ogpresence.csv   -out results_panneural_markers/anc.all.long.ogpresence.CP

# for cnidaria+bilateria tree:
python ../scripts/possvm_reconstruction.py -tree ../data/species_tree_CB.newick -ort results_panneural_markers/matrix.all.long.ogexpression.csv -out results_panneural_markers/anc.all.long.ogexpression.CB
python ../scripts/possvm_reconstruction.py -tree ../data/species_tree_CB.newick -ort results_panneural_markers/matrix.all.long.ogpresence.csv   -out results_panneural_markers/anc.all.long.ogpresence.CB

# for full polytomy tree:
python ../scripts/possvm_reconstruction.py -tree ../data/species_tree_PO.newick -ort results_panneural_markers/matrix.all.long.ogexpression.csv -out results_panneural_markers/anc.all.long.ogexpression.PO
python ../scripts/possvm_reconstruction.py -tree ../data/species_tree_PO.newick -ort results_panneural_markers/matrix.all.long.ogpresence.csv   -out results_panneural_markers/anc.all.long.ogpresence.PO
```

3. Visualise evolutionary histories:

```bash
Rscript s11_panneural_evol_2022-07-14.R
```

4. Functional enrichments of neural gene modules, based on *Mus musculus* orthologs:

```bash
Rscript s12_functional_enrichment_panneural_2022-12-05.R
Rscript s13_simple_plots_GOs_2022-12-20.R
```
