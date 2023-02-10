#### Input ####

# libraries
source("../scripts/helper.R")
library("umap")
library("phangorn")
library("ape")

# input
out_fn = "results_trees_cell_types/"
dir.create(out_fn, showWarnings = FALSE)

# species
list_species = c("Tadh","TrH2","Hhon","HoiH23")
spr = "TrH2"
spq = list_species [ !list_species %in% spr ]

# input
set_id = "global"
focid  = "cts"
focus  = "cell_type"
ann_fn = "../results_scatlas/results_metacell_it4/"
icc_fn = "results_alignment_icc/"
ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv"
mdb_fn = "../results_scatlas/data/scdb/"

# gene subsets to plot
list_subsets = c("tfs","siggpcr","sig","neuropept","hypNPs","top")

# log
message(sprintf("csps %s | %s UPGMA tree, expression matrix...", set_id, focid))

# read tree for order
tre_glo = ape::read.tree(sprintf("%s/csps.global.dendrogram.cts.UPGMA.newick", out_fn))
tre_pep = ape::read.tree(sprintf("%s/csps.peptidergic.dendrogram.cts.UPGMA.newick", out_fn))

# read matrix
csps_m = readRDS(sprintf("%s/csps.global.cspsmatrix.cts.rds", out_fn))

# global matrix
alg_t_is_tip = tre_glo$edge[,2] <= length(tre_glo$tip.label)
alg_t_ordered_tips = tre_glo$tip.label [ tre_glo$edge[alg_t_is_tip,2] ]
alg_fp = t(csps_m) [ rev(alg_t_ordered_tips), ]

# ignore peptidergic in global matrix
alg_fp_f = alg_fp [ !grepl("peptidergic", rownames(alg_fp)), ]

# create separate peptidergic matrix with row order from peptidergic tree
alp_fp_f = alg_fp [  grepl("peptidergic", rownames(alg_fp)), ]
alp_fp_f = alp_fp_f [ rev(tre_pep$tip.label) , ]

# load cell type data for all species
ann_cts = c()
for (spi in list_species) {
	ctt_spi_fn = sprintf(ctt_sprtinf_string, ann_fn, spi)
	ctt_spi = read.table(ctt_spi_fn, header = TRUE, comment.char = "", sep = "\t")
	# get color annotations for each species
	ann_spi = unique(ctt_spi[,c(focus,"color")])
	ann_spi = ann_spi [ !duplicated(ann_spi[,focus]), ]
	ann_spi_v = ann_spi[,2]
	names(ann_spi_v) = paste(spi, ann_spi[,1], sep = "|")
	# concatenate annotations
	ann_cts = c(ann_cts, ann_spi_v)
}


### Expression heatmaps ###

# heatmap
for (subset in list_subsets) {
	
	# check which genes have LRR repeats (only meaningful for GPCRs)
	pfa_r = read.table(sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spr), sep = "\t")
	pfa_r_lrr = pfa_r [ grepl("LRR", pfa_r[,2]) | grepl("fn2", pfa_r[,2]) | grepl("Ldl_recept_a",pfa_r[,2]) | grepl("CUB", pfa_r[,2]), ]
	lrr_genes_r = pfa_r_lrr[,1]
	lrr_color_d = rep("snow1", length(pfa_r[,1]))
	names(lrr_color_d) = pfa_r[,1]
	lrr_color_d [ names(lrr_color_d) %in% lrr_genes_r ] = "indianred1"
	names(lrr_color_d) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spr), vector_to_fix = names(lrr_color_d))
	
	## Global matrix ##	
	ali_fp_f = alg_fp
	
	# restrict to subset
	if (subset != "top") {
		
		clas_focus = read.table(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", subset, spr), sep = "\t")
		colnames(clas_focus) = c("transcript","annotation")
		clas_focus$gene = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spr), vector_to_fix = clas_focus$transcript)
		clas_focus = clas_focus [ !grepl("^NFYB_NFYC\\.", clas_focus$annotation), ]
		clas_focus_v = clas_focus$annotation
		names(clas_focus_v) = clas_focus$gene
		list_focus = unique(clas_focus$gene)
		list_focus = list_focus [ list_focus %in% colnames(ali_fp_f) ]
		
	} else {
		
		list_focus = colnames(ali_fp_f)
		
	}
	
	# subset
	ali_fp_f = ali_fp_f [ , list_focus ]

	if (!subset %in% c("neuropept","hypNPs")) {
		bool_variable = apply(ali_fp_f, 2, function(c) length(which(c >= 3)) >= 3 ) & apply(ali_fp_f, 2, function(c) sd(c)) >= 0.5
		# bool_constant = apply(ali_fp_f, 2, function(c) length(which(c >= 2)) >= 2 ) & apply(ali_fp_f, 2, function(c) sd(c)) <  0.5
		ali_fp_f = ali_fp_f [ , bool_variable ]
	}

	# gene order
	max_fc_allowed = ifelse(subset %in% c("neuropept","hypNPs"), 10, 10)
	ali_fp_f [ ali_fp_f > max_fc_allowed ] = max_fc_allowed
	gene_ord = order(apply(ali_fp_f, 2, function(c) which.max(rollmean(c, 3))))
	ali_fp_f = ali_fp_f [ , gene_ord ]
	
	# add OG name
	colnames(ali_fp_f) = paste(colnames(ali_fp_f), stringr::str_trunc(clas_focus_v, 60) [ colnames(ali_fp_f) ])
	lrr_color_d_i = lrr_color_d
	names(lrr_color_d_i) = paste(names(lrr_color_d_i), stringr::str_trunc(clas_focus_v, 60) [ names(lrr_color_d_i) ])
	lrr_color_d_i = lrr_color_d_i [ colnames(ali_fp_f) ]
	
	# plot
	plot_width = ceiling(ncol(ali_fp_f) / 8+ 4)
	plot_height = ceiling(nrow(ali_fp_f) / 8 + 4)
	hm = plot_complex_heatmap(
		ali_fp_f, name = "fp",
		fontsize = 5,
		cluster_row = FALSE,
		cluster_col = FALSE,
		color_min = 1.5,
		color_max = max_fc_allowed,
		colors_row = ann_cts,
		colors_col = lrr_color_d_i,
		use_raster = FALSE,
		color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
		cex_dotplot = 0.01,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		do_dotplot = FALSE,
		dot_size_min = 0.5, 
		dot_size_max = 4)
	pdf(sprintf("%s/expr.cts.all.%s.pdf", out_fn, subset), height = plot_height, width = plot_width)
	print(hm)
	dev.off()


	
	## Global matrix without peptidergic ##	
	ali_fp_f = alg_fp_f
	
	# restrict to subset
	if (subset != "top") {
		
		clas_focus = read.table(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", subset, spr), sep = "\t")
		colnames(clas_focus) = c("transcript","annotation")
		clas_focus$gene = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spr), vector_to_fix = clas_focus$transcript)
		clas_focus = clas_focus [ !grepl("^NFYB_NFYC\\.", clas_focus$annotation), ]
		clas_focus_v = clas_focus$annotation
		names(clas_focus_v) = clas_focus$gene
		list_focus = unique(clas_focus$gene)
		list_focus = list_focus [ list_focus %in% colnames(ali_fp_f) ]
		
	} else {
		
		list_focus = colnames(ali_fp_f)
		
	}
	
	# subset
	ali_fp_f = ali_fp_f [ , list_focus ]

	if (!subset %in% c("neuropept","hypNPs")) {
		bool_variable = apply(ali_fp_f, 2, function(c) length(which(c >= 3)) >= 3 ) & apply(ali_fp_f, 2, function(c) sd(c)) >= 0.5
		# bool_constant = apply(ali_fp_f, 2, function(c) length(which(c >= 2)) >= 2 ) & apply(ali_fp_f, 2, function(c) sd(c)) <  0.5
		ali_fp_f = ali_fp_f [ , bool_variable ]
	}

	# gene order
	max_fc_allowed = ifelse(subset %in% c("neuropept","hypNPs"), 10, 10)
	ali_fp_f [ ali_fp_f > max_fc_allowed ] = max_fc_allowed
	gene_ord = order( apply(ali_fp_f, 2, function(c) which.max(rollmean(c, 3))) )
	ali_fp_f = ali_fp_f [ , gene_ord ]
	
	# add OG name
	colnames(ali_fp_f) = paste(colnames(ali_fp_f), stringr::str_trunc(clas_focus_v, 60) [ colnames(ali_fp_f) ])
	lrr_color_d_i = lrr_color_d
	names(lrr_color_d_i) = paste(names(lrr_color_d_i), stringr::str_trunc(clas_focus_v, 60) [ names(lrr_color_d_i) ])
	lrr_color_d_i = lrr_color_d_i [ colnames(ali_fp_f) ]
	
	# plot
	plot_width = ceiling(ncol(ali_fp_f) / 8 + 4)
	plot_height = ceiling(nrow(ali_fp_f) / 8 + 4)
	pdf(sprintf("%s/expr.cts.no_pep.%s.pdf", out_fn, subset), height = plot_height, width = plot_width)
	hm = plot_complex_heatmap(
		ali_fp_f, name = "fp",
		fontsize = 5,
		cluster_row = FALSE,
		cluster_col = FALSE,
		color_min = 1.5,
		color_max = max_fc_allowed,
		colors_row = ann_cts,
		colors_col = lrr_color_d_i,
		use_raster = FALSE,
		color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
		cex_dotplot = 0.01,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		do_dotplot = FALSE,
		dot_size_min = 0.5, 
		dot_size_max = 4)
	print(hm)
	dev.off()
	


	## Peptidergic matrix ##	
	ali_fp_f = alp_fp_f
	
	# restrict to subset
	if (subset != "top") {
		
		clas_focus = read.table(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", subset, spr), sep = "\t")
		colnames(clas_focus) = c("transcript","annotation")
		clas_focus$gene = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spr), vector_to_fix = clas_focus$transcript)
		clas_focus = clas_focus [ !grepl("^NFYB_NFYC\\.", clas_focus$annotation), ]
		clas_focus_v = clas_focus$annotation
		names(clas_focus_v) = clas_focus$gene
		list_focus = unique(clas_focus$gene)
		list_focus = list_focus [ list_focus %in% colnames(ali_fp_f) ]
		
	} else {
		
		list_focus = colnames(ali_fp_f)
		
	}
	
	# subset
	ali_fp_f = ali_fp_f [ , list_focus ]

	if (!subset %in% c("neuropept","hypNPs")) {
		bool_variable = apply(ali_fp_f, 2, function(c) length(which(c >= 3)) >= 3 ) & apply(ali_fp_f, 2, function(c) sd(c)) >= 0.5
		# bool_constant = apply(ali_fp_f, 2, function(c) length(which(c >= 2)) >= 2 ) & apply(ali_fp_f, 2, function(c) sd(c)) <  0.5
		# markers_variable = names(which(bool_variable))
		# markers_constant = names(which(bool_constant))
		ali_fp_f = ali_fp_f [ , bool_variable ]
	}
	
	# gene order
	max_fc_allowed = ifelse(subset %in% c("neuropept","hypNPs"), 10, 10)
	ali_fp_f [ ali_fp_f > max_fc_allowed ] = max_fc_allowed
	gene_ord = order( as.numeric(apply(ali_fp_f, 2, function(c) sum(c >= 2)) >= ncol(ali_fp_f) / 3), apply(ali_fp_f, 2, function(c) which.max(rollmean(c, 3))) )
	ali_fp_f = ali_fp_f [ , gene_ord ]
	
	# add OG name
	colnames(ali_fp_f) = paste(colnames(ali_fp_f), stringr::str_trunc(clas_focus_v, 60) [ colnames(ali_fp_f) ])
	lrr_color_d_i = lrr_color_d
	names(lrr_color_d_i) = paste(names(lrr_color_d_i), stringr::str_trunc(clas_focus_v, 60) [ names(lrr_color_d_i) ])
	lrr_color_d_i = lrr_color_d_i [ colnames(ali_fp_f) ]
	
	# plot
	plot_width = ceiling(ncol(ali_fp_f) / 8 + 4)
	plot_height = ceiling(nrow(ali_fp_f) / 8 + 4)
	pdf(sprintf("%s/expr.cts.pept_only.%s.pdf", out_fn, subset), height = plot_height, width = plot_width)
	hm = plot_complex_heatmap(
		ali_fp_f, name = "fp",
		fontsize = 5,
		cluster_row = FALSE,
		cluster_col = FALSE,
		color_min = 1.5,
		color_max = max_fc_allowed,
		colors_row = ann_cts,
		colors_col = lrr_color_d_i,
		use_raster = FALSE,
		color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
		cex_dotplot = 0.01,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		do_dotplot = FALSE,
		dot_size_min = 0.5, 
		dot_size_max = 4)
	print(hm)
	dev.off()
	
	
}



### Num markers expressed ###

for (subset in list_subsets [ list_subsets != "top" ]) {
	
	# cell type order
	list_cts_glob = rownames(alg_fp)
	list_cts_pept = rownames(alp_fp_f)
	list_cts_notp = rownames(alg_fp_f)

	message(sprintf("num %s per species and ct", subset))
		
	pdf(sprintf("%s/nummarkers.cts.pept_only.%s.pdf", out_fn, subset), height = 12, width = 2.5)
	
	for (fp_thr in c(1.5, 2, 3, 10)) {
		
		count_per_ct_t = c()
		for (spi in list_species) {
			
			# load focus genes
			clas_focus = read.table(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", subset, spi), sep = "\t")
			colnames(clas_focus) = c("transcript","annotation")
			clas_focus$gene = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = clas_focus$transcript)
			clas_focus = clas_focus [ !grepl("^NFYB_NFYC\\.", clas_focus$annotation), ]
			clas_focus_v = clas_focus$annotation
			names(clas_focus_v) = clas_focus$gene
			list_focus = unique(clas_focus$gene)
			
			# load cell type fps
			metacell::scdb_init(mdb_fn, force_reinit = TRUE)
			ct = metacell::scdb_mc(sprintf("scdr_%s_it4_cts", spi))
			
			# subset and count
			ct_f = ct@mc_fp [ rownames(ct@mc_fp) %in% list_focus , ]
			count_per_ct = apply(ct_f, 2, function(v) {
				length(which(v >= fp_thr))
			})
			names(count_per_ct) = paste(spi, names(count_per_ct), sep = "|")
			count_per_ct_t = c(count_per_ct_t, count_per_ct)

		
		}
		
		b = barplot(count_per_ct_t [ c(rev(list_cts_pept), rev(list_cts_notp)) ], horiz = TRUE, las = 1, cex.names = 0.5, col = "lightblue3")
		text(0, b, sprintf("n=%i",count_per_ct_t [ c(rev(list_cts_pept), rev(list_cts_notp)) ]), cex = 0.5, pos = 4, col = "gray10")
		title(main = sprintf("num %s at fp>%.2f", subset, fp_thr), cex.main = 0.5)
	
	}

	dev.off()
	
}