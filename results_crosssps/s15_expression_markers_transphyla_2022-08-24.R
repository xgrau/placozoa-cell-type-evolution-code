#### Input ####

# libraries
source("../scripts/helper.R")
library("umap")
library("phangorn")
library("ape")

# input
heatmap_colors = c("white","orange","orangered2","#520c52")
out_fn = "results_ct_trees/"
out_fn = "results_trees_cell_types_transphyla/"
dir.create(out_fn, showWarnings = FALSE)

# species
list_species = c("Tadh","TrH2","Hhon","HoiH23","Nvec","Spis","Hvul")
# list_species = c("Tadh","TrH2","Hhon","HoiH23","Nvec","Spis")
list_outgroups = c("Nvec","Spis","Hvul")
spr = "TrH2"
spq = list_species [ !list_species %in% spr ]

# fc threshold
fc_thr = 1.5

# fraction of variance required in PCA based dendrogram
fraction_required = 0.4

# analyses
set_list = list(
	ann_fn = "../results_scatlas/results_metacell_it4/",
	icc_fn = "results_alignment_icc/",
	ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv",
	focus_list = list("bct" = "broad_cell_type"))



set_id = "all"
ann_fn = set_list$ann_fn
icc_fn = set_list$icc_fn
ctt_sprtinf_string = set_list$ctt_sprtinf_string
focus_list = set_list$focus_list

for (i in 1:length(focus_list)) {
	
	focid = names(focus_list)[[i]]
	focus = focus_list[[i]]

	# load cell type data for all species
	ann_cts = c()
	for (spi in list_species) {
		if (spi %in% c(list_outgroups)) {
			ctt_spi_fn = sprintf("../results_scatlas/data/scdb_outgroups/annot.%s.tsv", spi)
		} else {
			ctt_spi_fn = sprintf(ctt_sprtinf_string, ann_fn, spi)
		}
		ctt_spi = read.table(ctt_spi_fn, header = TRUE, comment.char = "", sep = "\t")
		# get color annotations for each species
		ann_spi = unique(ctt_spi[,c(focus,"color")])
		ann_spi = ann_spi [ !duplicated(ann_spi[,focus]), ]
		ann_spi_v = ann_spi[,2]
		names(ann_spi_v) = paste(spi, ann_spi[,1], sep = "|")
		# concatenate annotations
		ann_cts = c(ann_cts, ann_spi_v)
	}
	
	# read matrix
	csps_m = readRDS(sprintf("%s/csps.%s.cspsmatrix.%s.rds", out_fn, set_id, focid))
	
	# read tree for order
	ali_f_phy = ape::read.tree(sprintf("%s/csps.%s.dendrogram.%s.UPGMA.newick", out_fn, set_id, focid))
	
	# reorder matrix according to UPGMA tree
	message(sprintf("csps %s | %s UPGMA tree, expression matrix...", set_id, focid))
	ali_t = ali_f_phy
	ali_t_is_tip = ali_t$edge[,2] <= length(ali_t$tip.label)
	ali_t_ordered_tips = ali_t$tip.label [ ali_t$edge[ali_t_is_tip,2] ]
	ali_fp = t(csps_m) [ ali_t_ordered_tips, ]


	for (subset in c("tfs","siggpcr","sig","neuropept","top")) {
		
		# filter matrix
		ali_fp_f = ali_fp
		# ali_fp_f = ali_fp_f [ , apply(ali_fp_f, 2, function(c) length(which(c > 1.5)) > 3 ) ]

		# restrict to subset
		if (subset != "top") {
			
			clas_tfs = read.table(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", subset, spr), col.names = c("transcript","annotation"))
			clas_tfs$gene = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spr), vector_to_fix = clas_tfs$transcript)
			clas_tfs_v = clas_tfs$annotation
			names(clas_tfs_v) = clas_tfs$gene
			list_tfs = unique(clas_tfs$gene)
			list_tfs = list_tfs [ list_tfs %in% colnames(ali_fp_f) ]
			
		} else {
			
			list_tfs = colnames(ali_fp_f)
			
		}
		
		# subset
		ali_fp_f = ali_fp_f [ , list_tfs ]

		# gene order
		gene_ord = order(apply(ali_fp_f, 2, function(x) which.max(rollmean(x, 4))))
		ali_fp_f = ali_fp_f [ , gene_ord ]
		
		# drop genes that are not highly expressed in at least three tissues
		ali_fp_f = ali_fp_f [ , names(which(apply(ali_fp_f, 2, function(r) length(which(r > 1.5)) >= 4))) ]
		
		# # ec value
		# icc_ec_v = csps_i$og_pairs_ec_value
		# names(icc_ec_v) = csps_i$og_pairs[,1]
		# icc_fp_f = data.frame(rownames = colnames(ali_fp_f), ec_value = icc_ec_v [ colnames(ali_fp_f) ])
		# icc_fp_f$ec_value [ is.na(icc_fp_f$ec_value) ] = 0
		
		# add OG name
		colnames(ali_fp_f) = paste(clas_tfs_v [ colnames(ali_fp_f) ], colnames(ali_fp_f))
		
		# plot
		plot_height = ceiling(ncol(csps_m) / 6 + 4)
		plot_width = ceiling(ncol(ali_fp_f) / 6 + 4)
		pdf(sprintf("%s/csps.%s.expression.%s.UPGMA.%s.pdf", out_fn, set_id, focid, subset), height = plot_height, width = plot_width)
		ali_fp_f [ ali_fp_f > 4 ] = 4
		ali_fp_f_b = (ali_fp_f > fc_thr) * 1
		hm = plot_complex_heatmap(
			ali_fp_f_b,
			fontsize = 7,
			cluster_row = FALSE,
			cluster_col = TRUE,
			color_min = 0,
			color_max = 1,
			colors_row = ann_cts,
			use_raster = TRUE,
			color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
			do_dotplot = FALSE,
			cex_dotplot = 0.01,
			# cell_border = gpar(col = "white", lwd = 1, lty = 1),
			heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
			dot_size_min = 0.0, 
			dot_size_max = 4)
		print(hm)
		dev.off()
		
	}
	
	message(sprintf("csps %s | %s expression done!", set_id, focid))
	
}
