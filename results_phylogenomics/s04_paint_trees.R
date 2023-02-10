# libraries
library("ape")
library("phangorn")

# input variables
six_fn = "species_index.csv"
out_fn = "results_trees/"
dir.create(out_fn)


# list of trees to paint
list_trees = c(
#   list.files("results_broccoli_bcps/", pattern="*treefile", full.names=TRUE),
#   list.files("results_broccoli_metc/", pattern="*treefile", full.names=TRUE),
  list.files("results_broccoli_metc/", pattern="*treefile", full.names=TRUE),
  list.files("results_broccoli_meto/", pattern="*treefile", full.names=TRUE),
  list.files("results_broccoli_metc_fastev/", pattern="*treefile", full.names=TRUE),
  list.files("results_broccoli_metc_loci/", pattern="*\\.tree$", full.names=TRUE),
  list.files("results_broccoli_metc_phylobayes/bpcomp_output/", pattern="*\\.con.tre$", full.names=TRUE),
  list.files("results_broccoli_metc_phylobayes/bpcomp_output2/", pattern="*\\.con.tre$", full.names=TRUE),
  list.files("results_broccoli_meto_phylobayes/meto10_bpcomp_out/", pattern="*\\.con.tre$", full.names=TRUE)
)

# colors from taxonomy
six = read.table(six_fn, col.names = c("species", "full_name", "taxonomy"), sep = "\t")
six$color = "gray80"
six [ six$taxonomy == "Bilateria", "color" ] = "cyan4"
six [ six$taxonomy == "Cnidaria",  "color" ] = "cyan3"
six [ six$taxonomy == "Placozoa",  "color" ] = "chocolate2"
six [ six$taxonomy == "Ctenophora","color" ] = "royalblue1"
six [ six$taxonomy == "Porifera",  "color" ] = "royalblue3"
six [ six$taxonomy == "Choanoflagellata", "color" ] = "snow3"
color_per_species = six$color
names(color_per_species) = six$species

# full species names
name_per_species = six$full_name
names(name_per_species) = six$species

# loop and paint
for (tre_fn in list_trees) {

  # read tree  
  tre = ape::read.tree(tre_fn)
  tre = phangorn::midpoint(tre)
  tre = ape::ladderize(tre)
  
  # tree name
  tre_title = basename(tre_fn)
  tre_title = gsub(".treefile","", tre_title)
  message(sprintf("paint %s", tre_title))
  
  # short labels to long labels
  pdf(sprintf("%s/%s.pdf", out_fn, tre_title), width = 5, height = length(tre$tip.label) / 14 + 2)
  tre_tip_label_short = tre$tip.label
  tre$tip.label = name_per_species [ tre_tip_label_short ]
  
  # if support is 100%, do not plot
  if (!is.null(tre$node.label)) {
    tre$node.label [ tre$node.label == "100" ] = "*"
  }
  
  # plot
  ape::plot.phylo(tre, show.node.label = TRUE, tip.color = color_per_species [ tre_tip_label_short ], node.color = "gray50", font = 1, main = tre_title, cex.main = 0.4, cex = 0.4, lwd = 0.5, label.offset = 0.0)
  ape::add.scale.bar(cex=0.4, col = "gray50")
  dev.off()
  
}

