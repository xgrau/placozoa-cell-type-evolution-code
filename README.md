# Cell type evolution and diversity in placozoans

This repository contains all necessary scripts to reproduce the single-cell atlases of the placozoans *Trichoplax adhaerens* (H1), *Trichoplax* sp. (H2), *Hoilungia hongkongensis* (H13), and *Cladtertia collaboinventa* (H23), presented in our manuscript [here](https://github.com/sebepedroslab/placozoa-cell-type-evol).

![placozoan tree of life](data/fig_tree.png)

## Analyses

Check instructions in each of these folders.

1. **Single-cell atlases** for each species, in the [`results_scatlas/`](results_scatlas/) folder. This covers all steps from mapping to the final `metacell` clustering and cell type classification for each species, as well as creating expression maps, 2D projections, and dedicated expression plots for certain markers (TFs, GPCRs, etc.).
2. **Cross-species comparisons** within placozoa and with other metazoan species, in [`results_crosssps/`](results_crosssps/).
3. **Gene module analyses**, in [`results_gene_modules/`](results_gene_modules/). This includes the analysis of placozoan-specific gene modules, and the analyses of neural-related gene modules across metazoans.g
4. **Phylogenomics** in [`results_phylogenomics/`](results_phylogenomics). Includes instructions to reproduce our single-gene marker phylogenomics datasets, and run the species trees with IQ-TREE and Phylobayes.

## Browse the data

You can browse the results from our analyses in [this ShinyApp](https://sebelab.crg.eu/placozoa_cell_atlas/).

![snapshot of the database](data/fig_snap.png)

## Access the data

We provide the following [raw data](results_scatlas/data/scdb) for each species:

* **UMI matrices**, which are the files starting with the **`mat` prefix**. For all analyses, we used matrices indicated with the `it2` suffix (previous iterations correspond to intermediate steps). These are sparse matrices stored as R objects, specifically [metacell package](https://tanaylab.github.io/metacell/)-type `mat` objects. See below for an example of how to load these matrices into R.
* **Cell-metacell assignment files** indicated with the **`mc` prefix**. For all analyses, we used matrices indicated with the `it4` suffix (previous iterations correspond to intermediate steps).
* **Two-dimensional porjection** objects indicated in the **`mc2d` prefix**.
* Other objects used by metacell such as gene statistics, gene sets, cell coclustering data... Please refer to the appropriate scripts in the relevant [readme file](results_scatlas/README.md`) for details.

Each species is referred to with a short acronym, all through this project:

* `Tadh` for *Trichoplax adhaerens* H1
* `TrH2` for *Trichoplax* sp. H2
* `Hhon` for *Hoilungia hongkongensis* H13
* `HoiH23` for *Cladtertia collaboinventa* H23 (we came up with the label before it was [recently renamed](https://www.frontiersin.org/articles/10.3389/fevo.2022.1016357/full) from a species in the *Hoilungia* genus to a genus of its own, apologies for the confusion).

We also provide **metacell-cell type assignment tables**, with our cell type annotations for each metacell. The final annotations are stored in this folder: [`results_scatlas/results_metacell_it4/`](results_scatlas/results_metacell_it4/), e.g. [this file](results_scatlas/results_metacell_it4/annotation_mc.Tadh.it4.reordered.tsv).

The [metacell package]([url](https://tanaylab.github.io/metacell/)) contains analogous functions to load other `metacell`-formatted data files, such as cell clusterings in `mc` objects (which specify which cells in the UMI matrix belong to each metacell cluster).

How to load metacell objects in R:

```R
# from the `results_scatlas/` folder in this repository
# load metacell library
library("metacell")
# initialise metacell database using its relative path:
metacell::scdb_init("data/scdb/", force_reinit = TRUE)

# read the mat object using its ID to refer to it (in this case, the `scdr_Tadh_it2` bit); there ara analogous files for metacells, etc
mat = metacell::scdb_mat("scdr_Tadh_it2")
# sparse matrix available in the mat@mat slot, in this case it contains 16386 genes x 13236 cells
dim(mat@mat)

# likewise, you may want to load a mc object containing cell-to-metacell assignments
mc = metacell::scdb_mc("scdr_Tadh_it4")
# the mc@mc slot is a vector with all cells and their associated metacell
length(mc@mc)
# notice that this vector contains only 13151 cells: not all cells in the mat@mat matrix are classified into a metacell

# if you want to know which cell type is annotated to each metacell, check the mc annotation file
# make sure that you load the same version of the metacell clustering as in the mc object, in this case, it4
ctt = read.table("results_metacell_it4/annotation_mc.Tadh.it4.reordered.tsv", header = TRUE)
head(ctt)
# metacell	cell_type	color	broad_cell_type	metacell_it2
# 1	lipophil	khaki2	lipophil	1
# 2	lipophil	khaki2	lipophil	2
# 3	lipophil	khaki2	lipophil	3
# 4	lipophil	khaki2	lipophil	4
# 5	lipophil	khaki2	lipophil	5
```

If you have any queries, feel free to let me know in the [Issues section](https://github.com/xgrau/placozoa-cell-type-evolution/issues).

Enjoy!

```python

    _,,.- ~~--~~~~-.._,
  /                     `---.
  \                         .'
   `~- ,_ ,. ,_, ,., ,., -··
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
```
