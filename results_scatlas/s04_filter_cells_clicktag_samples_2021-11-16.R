# libraries
library("metacell")
source("../scripts/helper.R")
library("Seurat")
library("uwot")
par(family  = "Arial")

table_to_matrix = function(table) {
    mat = matrix(table, nrow = nrow(table))
    rownames(mat) = rownames(table)
    colnames(mat) = colnames(table)
    return(mat)
}

# libraries to analyse
libraries_pairs = list(
	list(
		sps_list =    c("Tadh", "Hhon"),
		lib = "H1H13_1_ACME_CT_10x_10kc",
		ctmat_fn = "data_clicktag/map_H1H13/",
		ctbcs_fn = "data_clicktag/barcodes_clicktag_H1H13.pairs.csv"
	),
	list(
		sps_list =    c("Tadh", "HoiH23"),
		lib = "H1H23_5_ACME_10x_10kc",
		ctmat_fn = "data_clicktag/map_H1H23/",
		ctbcs_fn = "data_clicktag/barcodes_clicktag_H1H23.pairs.csv"
	),
	list(
		sps_list =    c("TrH2", "HoiH23"),
		lib = "H2H23_4_ACME_10x_10kc",
		ctmat_fn = "data_clicktag/map_H2H23/",
		ctbcs_fn = "data_clicktag/barcodes_clicktag_H2H23.pairs.csv"
	),
	list(
		sps_list =    c("Tadh", "Hhon"),
		lib = "Plac01_H1_H13_10XscRNAseq_10kc",
		ctmat_fn = "data_clicktag/map_Plac01H1H13/",
		ctbcs_fn = "data_clicktag/barcodes_clicktag_Plac01H1H13.pairs.csv"
	),
	list(
		sps_list =    c("TrH2", "HoiH23"),
		lib = "Plac02_H2_H23_10XscRNAseq_10kc",
		ctmat_fn = "data_clicktag/map_Plac02H2H23/",
		ctbcs_fn = "data_clicktag/barcodes_clicktag_Plac02H2H23.pairs.csv"
	)
)

# thresholds
cz_large_thr = 15000     # cells smaller or larger than this (in terms of expression) are considered non-cells
cz_small_thr = 100       
ct_small_thr = 20        # cells with less than X umis in the most frequent clicktag are marked as unclassified
ct_fc_thr = 1.5          # min threshold of first-to-third clicktag barcodes (ratio)
seurat_purity_frac = 0.7 # this fraction determines the % of cells from a given clicktag class in a seurat cluster below which the cluster is blacklisted

for (i in 1:length(libraries_pairs)) {
	
	#### Load data ####
	
	# get info for this pair of species
	sps_list = libraries_pairs[[i]]$sps_list
	sps_list_ct = libraries_pairs[[i]]$sps_list_ct
	clicktag_lib = libraries_pairs[[i]]$clicktag_lib
	lib = libraries_pairs[[i]]$lib
	sp1 = sps_list[1]
	sp2 = sps_list[2]
	
	# objects
	mat = list()
	mat_cz = list()
	for (spi in sps_list) {
		
		# init database
		metacell::scdb_init("data/scdb/",force_reinit=TRUE)
		
		# load umi table
		run_name = sprintf("scdr_%s", spi)
		mat[[spi]] = metacell::scdb_mat(run_name)
		mat[[spi]]@mat = mat[[spi]]@mat [ , mat[[spi]]@cell_metadata$seq_batch_id %in% lib ]
		mat_cz[[spi]] = Matrix::colSums(mat[[spi]]@mat)
		names(mat_cz[[spi]]) = cell_name_clean(names(mat_cz[[spi]]))
		
	}
	
	# get cell sizes of shared barcodes (cells from this library present in two species)
	all_bcs = unique(c(names(mat_cz[[1]]), names(mat_cz[[2]])))
	mat_cz_s1 = mat_cz[[1]] [ all_bcs ]
	mat_cz_s2 = mat_cz[[2]] [ all_bcs ]
	names(mat_cz_s1) = all_bcs
	names(mat_cz_s2) = all_bcs
	mat_cz_s1 [ is.na(mat_cz_s1) ] = 0
	mat_cz_s2 [ is.na(mat_cz_s2) ] = 0
	
	
	# init cell barcode classification dataframe
	cla = data.frame(row.names = names(mat_cz_s1))
	
	
	#### Open plot ####
	
	pdf(sprintf("data/dual_samples.%s.dist.pdf",lib), height = 10, width = 10)
	par(mfrow=c(2,2))
	
	#### Expression bias classification ####
	
	# get relative cell sizes	
	relative_sizes = log2(mat_cz_s1 / mat_cz_s2)
	
	# histograms to decide thresholds
	# all barcodes
	hh = hist(
		relative_sizes, 
		breaks = 60,
		border = NA,
		col = "gray90",
		las = 1,
		main = sprintf("Cell size %s-%s\n%s", sp1, sp2, lib),
		xlab = sprintf("log2(UMI %s / UMI %s)", sp1, sp2),
		xlim = c(-4,4))
	# cells and non-cells
	hh = hist(
		relative_sizes [ mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr ], 
		breaks = 60,
		border = NA,
		col = alpha("springgreen3", 0.4),
		add = T)
	relative_sizes_noncells = relative_sizes [ !(mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr) ]
	if (length(relative_sizes_noncells) == 0)  { relative_sizes_noncells = 0 }
	hh = hist(
		relative_sizes_noncells, 
		breaks = 60,
		border = NA,
		col = alpha("gray50", 0.7),
		add = T)
	legend("topleft", c("all barcodes",sprintf("cells >%i UMI in either sps", cz_small_thr), "non-cells"), fill = c("gray","springgreen3","gray50"), border = NA, bty = "n", cex = 0.9)
	
	# select cells
	thr_up = log2(1.25)
	thr_do = log2(1.25) * -1
	abline(v = c(thr_do,0,thr_up), lty = 2, col = c("springgreen4","black","springgreen4"))
	whitelist_s1 = names(relative_sizes) [ relative_sizes >= thr_up  & mat_cz_s1 >= cz_small_thr ]
	whitelist_s2 = names(relative_sizes) [ relative_sizes <= thr_do  & mat_cz_s2 >= cz_small_thr ]
	doublet_s1s2 = names(relative_sizes) [ relative_sizes > thr_do & relative_sizes < thr_up & mat_cz_s2 > cz_small_thr ]
	real_cells  = c(whitelist_s1, whitelist_s2, doublet_s1s2)
	
	# add cell size-based classification
	cla$expression_bias_label = "empty"
	cla$expression_bias_label [ rownames(cla) %in% whitelist_s1 ] = sp1
	cla$expression_bias_label [ rownames(cla) %in% whitelist_s2 ] = sp2
	cla$expression_bias_label [ rownames(cla) %in% doublet_s1s2 ] = "doublet"
	cla$expression_bias_label = factor(cla$expression_bias_label, levels = c(sp1,sp2,"doublet","empty"))
	
	#### Clicktag classification based on first and third tags ####
	
	# init database
	dir.create("data/scdb_clicktag/", showWarnings = TRUE)
	metacell::scdb_init("data/scdb_clicktag/",force_reinit=TRUE)
	
	# import clicktag counts
	ctmat_fn = libraries_pairs[[i]]$ctmat_fn
	ctbcs_fn = libraries_pairs[[i]]$ctbcs_fn
	mcell_import_scmat_kallisto(
		mat_nm = basename(ctmat_fn),
		matrix_fn = sprintf("%s/genes.mtx", ctmat_fn), 
		genes_fn =  sprintf("%s/genes.genes.txt", ctmat_fn), 
		cells_fn =  sprintf("%s/genes.barcodes.txt", ctmat_fn) 
	)
	
	# load clicktag counts
	mct = metacell::scdb_mat(basename(ctmat_fn))
	bct = read.table(ctbcs_fn, sep = ",", col.names = c("name","barcode","barcode_pair"))
	bct_pairs = bct$barcode_pair
	names(bct_pairs) = bct$name
	
	# normalise clicktag counts	
	mct_o = as.data.frame(as.matrix(mct@mat) [ all_bcs [ all_bcs %in% rownames(mct@mat) ] , ])
	# mct_o = mct_o [ , !grepl("_XX$", colnames(mct_o)) ]
	mct_o = mct_o [ all_bcs , ]
	rownames(mct_o) = all_bcs
	mct_o = as.matrix(mct_o)
	mct_n = t(t(mct_o) / Matrix::colSums(mct_o, na.rm = TRUE)) * 1e4
	# mct_n = mct_n / rowSums(mct_n, na.rm = TRUE)
	
	# some sort of observed/expected normalisation that isn't working at all
	# mct_o_colsum = Matrix::colSums(mct_o, na.rm = TRUE)
	# mct_o_rowsum = Matrix::rowSums(mct_o, na.rm = TRUE)
	# mct_e = matrix(unlist(lapply(mct_o_rowsum, function(x) x * mct_o_colsum)), nrow = nrow(mct_o))
	# mct_n = mct_o / (mct_e / sum(mct_o, na.rm = T))
	# mct_n [ is.na(mct_n) ] = 0
	
	# get labels and normalised counts from first and third most abundant clicktags
	mct_t = data.frame(
		first_tag   = apply(mct_n, 1,  function(x) colnames(mct_n) [ order(x, decreasing = TRUE) ] [1]), 
		first_val   = apply(mct_n, 1,  function(x) x [ order(x, decreasing = TRUE) ] [1]),
		first_count = apply(mct_o, 1,  function(x) x [ order(x, decreasing = TRUE) ] [1]),
		second_tag  = apply(mct_n, 1, function(x) colnames(mct_n) [ order(x, decreasing = TRUE) ] [2]), 
		second_val  = apply(mct_n, 1, function(x) x [ order(x, decreasing = TRUE) ] [2]),
		third_tag   = apply(mct_n, 1,  function(x) colnames(mct_n) [ order(x, decreasing = TRUE) ] [3]),
		third_val   = apply(mct_n, 1,  function(x) x [ order(x, decreasing = TRUE) ] [3])
	)
	
	# general class of the most abundant clicktag barcode per cell
	mct_t$first_class = gsub(".*_", "", mct_t$first_tag)
	mct_t$second_class = gsub(".*_", "", mct_t$second_tag)
	mct_t$third_class = gsub(".*_", "", mct_t$third_tag)
	# which are the two most common classes for the most abundant barcodes?
	tag_s_i = names(sort(table(mct_t$first_class [ rownames(mct_t) %in% c(whitelist_s1, whitelist_s2) ]), decreasing = TRUE)[1])
	tag_s_j = names(sort(table(mct_t$first_class [ rownames(mct_t) %in% c(whitelist_s1, whitelist_s2) ]), decreasing = TRUE)[2])
	
	# add barcode pair info
	mct_t$first_tag_pair = bct_pairs [ mct_t$first_tag ]
	mct_t$second_tag_pair = bct_pairs [ mct_t$second_tag ]
	mct_t$third_tag_pair = bct_pairs [ mct_t$third_tag ]
	
	# total size
	mct_t$max_counts = apply(mct_o, 1, max)
	mct_t$max_counts [ is.na(mct_t$max_counts) ] = 0
	
	# relative cell sizes
	# they are calculated by adding a constant eps factor to each normalised value
	eps_f = quantile(mct_t$first_val, 0.01, na.rm = TRUE)
	eps_s = quantile(mct_t$second_val, 0.01, na.rm = TRUE)
	eps_t = quantile(mct_t$third_val, 0.01, na.rm = TRUE)
	# eps_f = 0
	# eps_s = 0
	# eps_t = 0
	mct_t$relative_size_ft = (mct_t$first_val + eps_f) / (mct_t$third_val + eps_t)
	mct_t$relative_size_fs = (mct_t$first_val + eps_f) / (mct_t$second_val + eps_s)
	mct_t$fract_ft = mct_t$first_val / (mct_t$third_val + mct_t$first_val)
	
	# do the first, second and third tags come from the same pair of clicktags?
	mct_t$concordant_ft_tags = mct_t$first_tag_pair == mct_t$third_tag_pair
	mct_t$concordant_fs_tags = mct_t$first_tag_pair == mct_t$second_tag_pair
	mct_t$concordant_ft_clas = mct_t$first_class    == mct_t$third_class
	mct_t$concordant_fs_clas = mct_t$first_class    == mct_t$second_class
	
	# grouping vector of clicktag classes (most common ones only), to aggregate count matrices
	bct_v = ifelse(
		grepl(sprintf("_%s$",tag_s_i), colnames(mct_o)), 
		tag_s_i, 
		ifelse(grepl(sprintf("_%s$",tag_s_j), colnames(mct_o)),
					 tag_s_j, 
					 "other"))
	mct_o_a = sca_mc_gene_counts_noobj(mat = mct_o, grouping_vector = bct_v, T_totumi = 0)
	mct_o_a = as.matrix(as.data.frame(mct_o_a) [ all_bcs, ])
	mct_t$fract_agg = mct_o_a[,1] / (mct_o_a[,1] + mct_o_a[,2])
	
	# add clicktag-based annotation
	cla$clicktag_label = "unclassified"
	cla$clicktag_label [ mct_t$concordant_fs_tags &  !mct_t$concordant_ft_tags & mct_t$max_counts > ct_small_thr  & cla$expression_bias_label != "empty" & mct_t$relative_size_ft > ct_fc_thr ] = mct_t$first_class [ mct_t$concordant_fs_tags & !mct_t$concordant_ft_tags & mct_t$max_counts > ct_small_thr & cla$expression_bias_label != "empty" & mct_t$relative_size_ft > ct_fc_thr ]
	cla$clicktag_label [ !mct_t$concordant_fs_tags &  mct_t$max_counts > ct_small_thr &  mct_t$concordant_ft_clas & cla$expression_bias_label != "empty" ] = "doublet_intrasps"
	cla$clicktag_label [ !mct_t$concordant_fs_tags &  mct_t$max_counts > ct_small_thr & !mct_t$concordant_ft_clas & cla$expression_bias_label != "empty" ] = "doublet_intersps"
	cla$clicktag_label [ cla$expression_bias_label == "empty" ] = "empty"
	# as factor
	cla$clicktag_label = factor(cla$clicktag_label, levels = c(tag_s_i, tag_s_j, "doublet_intrasps","doublet_intersps","unclassified","empty"))
	cla$clicktag_label [ is.na(cla$clicktag_label) ] = "unclassified"
	table(cla$clicktag_label, useNA = "ifany")
	
	
	# histograms first to third tag
	# all cells
	hh = hist(
		log10(mct_t$relative_size_ft) [ mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr ] , 
		breaks = 60,
		xlim = c(0, 3),
		border = NA,
		col = alpha("springgreen3", 0.5),
		las = 1,
		main = sprintf("Normalised clicktags %s-%s\n%s", sp1, sp2, basename(ctmat_fn)),
		xlab = sprintf("log10(first ctbc / third ctbc)"))
	# concordant or discordant tags
	hh = hist(
		log10(mct_t$relative_size_ft) [ (mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr ) & ! mct_t$concordant_ft_tags ] , 
		breaks = 60,
		border = NA,
		col = alpha("blue", 0.3),
		add = T)
	hh = hist(
		log10(mct_t$relative_size_ft) [ (mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr ) & ! mct_t$concordant_fs_tags ] , 
		breaks = 60,
		border = NA,
		col = alpha("firebrick3", 0.7),
		add = T)	
	legend("topright", c("cells","cells with discordant 1st/3rd tags", "cells with discordant 1st/2nd tags"), fill = c("springgreen3","blue","firebrick3"), border = NA, bty = "n", cex = 0.9, title = "clicktag class")
	
	# histograms first to third tag
	# all cells
	hh = hist(
		mct_t$fract_agg [ mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr ] , 
		breaks = 60,
		xlim = c(0, 1),
		border = NA,
		col = alpha("springgreen3", 0.5),
		las = 1,
		main = sprintf("Normalised clicktags %s-%s\n%s", sp1, sp2, basename(ctmat_fn)),
		xlab = sprintf("Fraction most common tag in ctbc (%s/%s+%s)", tag_s_i, tag_s_i, tag_s_j))
	# concordant or discordant tags
	hh = hist(
		mct_t$fract_agg [ (mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr ) & ! mct_t$concordant_ft_tags ] , 
		breaks = 60,
		border = NA,
		col = alpha("blue", 0.3),
		add = T)
	hh = hist(
		mct_t$fract_agg [ (mat_cz_s1 >= cz_small_thr | mat_cz_s2 >= cz_small_thr ) & ! mct_t$concordant_fs_tags ] , 
		breaks = 60,
		border = NA,
		col = alpha("firebrick3", 0.7),
		add = T)	
	legend("topleft", c("cells","cells with discordant 1st/3rd tags", "cells with discordant 1st/2nd tags"), fill = c("springgreen3","blue","firebrick3"), border = NA, bty = "n", cex = 0.9, title = "clicktag class")
	
	#### Clicktag classification based on paired tags ####
	
	# # grouping vector of clicktag pairs
	# bct_vp = paste("p", bct$barcode_pair, "_", bct_v, sep = "")
	# mct_o_ap = sca_mc_gene_counts_noobj(mat = mct_o, grouping_vector = bct_vp, T_totumi = 0)
	# mct_o_ap = as.data.frame(mct_o_ap) [ rownames(mct_o),  ]
	# mct_n_ap = t(t(mct_o_ap) / Matrix::colSums(mct_o_ap, na.rm = TRUE) * 1e4)
	# 
	# # get labels and normalised counts from first and third most abundant clicktags
	# mct_tp = data.frame(
	# 	first_tag   = apply(mct_n_ap, 1,  function(x) colnames(mct_n_ap) [ order(x, decreasing = TRUE) ] [1]),
	# 	first_val   = apply(mct_n_ap, 1,  function(x) x [ order(x, decreasing = TRUE) ] [1]),
	# 	first_count = apply(mct_o_ap, 1,  function(x) x [ order(x, decreasing = TRUE) ] [1]),
	# 	second_tag  = apply(mct_n_ap, 1, function(x) colnames(mct_n_ap) [ order(x, decreasing = TRUE) ] [2]),
	# 	second_val  = apply(mct_n_ap, 1, function(x) x [ order(x, decreasing = TRUE) ] [2])
	# )
	# 
	# # general class of the most abundant clicktag barcode per cell
	# mct_tp$first_class = gsub(".*_", "", mct_tp$first_tag)
	# mct_tp$second_class = gsub(".*_", "", mct_tp$second_tag)
	# 
	# # total size
	# mct_tp$max_counts = apply(mct_o_ap, 1, max)
	# mct_tp$max_counts [ is.na(mct_t$max_counts) ] = 0
	# 
	# # relative cell sizes
	# # they are calculated by adding a constant eps factor to each normalised value
	# eps_f = quantile(mct_tp$first_val, 0.01, na.rm = TRUE)
	# eps_s = quantile(mct_tp$second_val, 0.01, na.rm = TRUE)
	# # eps_f = 0
	# # eps_s = 0
	# mct_tp$relative_size_fs = (mct_t$first_val + eps_f) / (mct_t$second_val + eps_s)
	# mct_tp$relative_size_fs [ is.na(mct_tp$relative_size_fs) ] = 0
	# mct_tp$fract_ft = mct_t$first_val / (mct_t$third_val + mct_t$first_val)
	# 
	# #### NEW: based on first-to-second of clicktag pairs
	# cla$clicktag_label = "unclassified"
	# cla$clicktag_label [ mct_tp$relative_size_fs > ct_fc_thr & cla$expression_bias_label != "empty" & mct_tp$max_counts > ct_small_thr ] = mct_tp$first_class [ mct_tp$relative_size_fs > ct_fc_thr & cla$expression_bias_label != "empty" & mct_tp$max_counts > ct_small_thr ]
	# cla$clicktag_label [ cla$expression_bias_label == "empty" ] = "empty"
	# 

	#### Scatter plots ####
	
	# cells sorted by clicktag max value
	mct_o_max = apply(mct_o, 1, max)
	mct_o_max_s = sort(mct_o_max, decreasing = TRUE)
	plot(
		mct_o_max_s,
		col = ifelse(names(mct_o_max_s) %in% real_cells, alpha("springgreen3", 0.2), alpha("gray",0.2)),
		cex = 0.5,
		ylim = c(1,15000),
		ylab = "Max ctbc count",
		xlab = "Cells",
		main = "Cells sorted by clicktag max count",
		sub = sprintf("n=%i", length(mct_o_max_s)),
		las = 1,
		log = "y")
	legend("topright", c("non-cells","cells"), cex = 0.9, col = c("gray","springgreen3"), pch = 1, bty = "n")
	
	# preapre legends
	legend_t = table(cla$expression_bias_label)
	legend_c = table(cla$clicktag_label)
	
	# expression plot with expression class	
	plot(
		mat_cz_s1 , mat_cz_s2, 
		log = "xy", 
		col = alpha(c("blue","orange","magenta3","gray"), 0.4) [ cla$expression_bias_label ],
		cex = 0.5, 
		las = 1,
		main = sprintf("Cell size %s-%s\n%s", sp1, sp2, lib),
		xlim = c(1,cz_large_thr),
		ylim = c(1,cz_large_thr),
		xlab = sp1,
		ylab = sp2)
	abline(a = 0, b = 1, lty = 2)
	legend("topleft", sprintf("%s = %i", names(legend_t), legend_t), col = c("blue","orange","magenta3","gray"), pch = 1, bty = "n", title = "expr bias classification", cex = 0.9)
	
	# expression plot with clicktag class	
	plot(
		mat_cz_s1 , mat_cz_s2, 
		log = "xy", 
		col = alpha(c("blue","orange","chartreuse3", "magenta3","gray50","gray"), 0.4) [ cla$clicktag_label ],
		cex = 0.5, 
		las = 1,
		main = sprintf("Cell size %s-%s\n%s", sp1, sp2, lib),
		xlim = c(1,cz_large_thr),
		ylim = c(1,cz_large_thr),
		xlab = sp1,
		ylab = sp2)
	abline(a = 0, b = 1, lty = 2)
	legend("topleft", sprintf("%s = %i", names(legend_c), legend_c), col = c("blue","orange","chartreuse3", "magenta3","gray50","gray"), pch = 1, bty = "n", title = "clicktag classification", cex = 0.9)
	
	# clicktag plot with expression class	
	plot(
		mct_o_a[,tag_s_i],
		mct_o_a[,tag_s_j],
		col = alpha(c("blue","orange","magenta3","gray"), 0.4) [ cla$expression_bias_label ],
		xlab = sprintf("%s clicktag counts", tag_s_i),
		ylab = sprintf("%s clicktag counts", tag_s_j),
		main = sprintf("Norm counts in most common clicktag bcs (%s&%s)", tag_s_i, tag_s_j),
		xlim = c(1,cz_large_thr),
		ylim = c(1,cz_large_thr),
		log = "xy",
		las = 1,
		cex = 0.5)
	abline(a=0, b=1, lty=2)
	legend("topleft", sprintf("%s = %i", names(legend_t), legend_t), col = c("blue","orange","magenta3","gray"), pch = 1, bty = "n", title = "expr bias classification", cex = 0.9)
	
	# clicktag plot with clicktag class
	plot(
		mct_o_a[,tag_s_i],
		mct_o_a[,tag_s_j],
		col = alpha(c("blue","orange","chartreuse3", "magenta3","gray50","gray"), 0.4) [ cla$clicktag_label ],
		xlab = sprintf("%s clicktag counts", tag_s_i),
		ylab = sprintf("%s clicktag counts", tag_s_j),
		main = sprintf("Norm counts in most common clicktag bcs (%s&%s)", tag_s_i, tag_s_j),
		xlim = c(1,cz_large_thr),
		ylim = c(1,cz_large_thr),
		log = "xy",
		las = 1,
		cex = 0.5)
	abline(a=0, b=1, lty=2)
	legend("topleft", sprintf("%s = %i", names(legend_c), legend_c), col = c("blue","orange","chartreuse3", "magenta3","gray50","gray"), pch = 1, bty = "n", title = "clicktag classification", cex = 0.9)
	
	
	#### Seurat classification ####

	# input data, without non-cells
	mct_o_clean = mct_o [ apply(mct_o, 1, function(x) any(!is.na(x))) & cla$expression_bias_label != "empty" , ]
	seu_o = Seurat::CreateSeuratObject(t(mct_o_clean))
	seu_o = Seurat::NormalizeData(seu_o, normalization.method="LogNormalize", scale.factor=10000)
	Seurat::VariableFeatures(seu_o) = colnames(mct_o)
	seu_o = Seurat::FindVariableFeatures(seu_o)
	seu_o = Seurat::ScaleData(seu_o)

	# find clusters
	seu_o = Seurat::RunPCA(seu_o, approx=FALSE)
	seu_o = Seurat::FindNeighbors(seu_o, dims=1:ncol(mct_o), k.param=50)
	seu_o = Seurat::FindClusters(seu_o, resolution=1, algorithm=1)
	seu_o = Seurat::RunUMAP(seu_o, features = VariableFeatures(seu_o))

	# plot seurat umap and clusters
	seu_umap = seu_o@reductions$umap@cell.embeddings
	seu_umap = as.data.frame(seu_umap) [ rownames(cla) ,  ]
	rownames(seu_umap) = rownames(cla)
	
	# seurat clusters
	# get clusters and reorder them to original cell order
	seu_clu = seu_o@meta.data$seurat_clusters
	names(seu_clu) = rownames(seu_o@meta.data)
	seu_clu = seu_clu [ rownames(cla) ]
	# give them good colors
	color_palette = colorRampPalette(c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen","deepskyblue","cadetblue1","mediumblue","darkviolet","violet"))
	color_palette_cols = color_palette(n=length(unique(seu_clu)))
	seu_clu_col = color_palette_cols [ seu_clu ]
	
	# expression bias class
	plot(seu_umap, cex = 0.5, col = alpha(c("blue","orange","magenta3","gray"), 0.4) [ cla$expression_bias_label ], main = "Seurat UMAP, by expr class (cells only)")
	legend("topleft", sprintf("%s = %i", names(legend_t), legend_t), col = c("blue","orange","magenta3","gray"), pch = 1, bty = "n", title = "expr bias classification", cex = 0.6)

	# clicktag class
	plot(seu_umap, cex = 0.5, col = alpha(c("blue","orange","chartreuse3", "magenta3","gray50","gray"), 0.4) [ cla$clicktag_label ], main = "Seurat UMAP, by clicktag class (cells only)")
	legend("topleft", sprintf("%s = %i", names(legend_c), legend_c), col = c("blue","orange","chartreuse3", "magenta3","gray50","gray"), pch = 1, bty = "n", title = "clicktag classification", cex = 0.6)
	
	# seurat clusters
	plot(seu_umap, cex = 0.5, col = alpha(seu_clu_col,0.2), main = "Seurat UMAP, by seurat cluster (cells only)")
	legend_s = sprintf("c%s, n = %i", levels(seu_clu), table(seu_clu))
	legend("topleft", cex = 0.6, col = color_palette_cols, legend = legend_s, pch = 1, bty = "n", ncol = 3, title = "Seurat clu")
	
	# original bcs
	color_palette_cols = color_palette(n=ncol(mct_o))
	bct_fac = factor(mct_t$first_tag, levels = colnames(mct_o))
	bct_clu_col = color_palette_cols [ bct_fac ]
	plot(seu_umap, cex = 0.5, col = alpha(bct_clu_col, 0.4), main = "Seurat UMAP, by ctbc (cells only)")
	legend("topleft", levels(bct_fac), col = color_palette_cols, pch = 1, bty = "n", title = "bcs", cex = 0.6, ncol = 2)

	
	#### Cell size correlations ####
	
	# cell size correlations
	# cor value
	mat_cz_max = pmax(mat_cz_s1, mat_cz_s2) [ cla$expression_bias_label != "empty" ] 
	mct_cz_max = mct_o_max [ cla$expression_bias_label != "empty" ]
	mat_cz_max = mat_cz_max [ !is.na(mct_cz_max) ]
	mct_cz_max = mct_cz_max [ !is.na(mct_cz_max) ]
	cz_cor = cor(mat_cz_max, mct_cz_max)
	# colored by expression bias label
	plot(
		pmax(mat_cz_s1, mat_cz_s2) [ cla$expression_bias_label != "empty" ], 
		mct_o_max [ cla$expression_bias_label != "empty" ], 
		log = "xy", 
		# col = alpha("blue", 0.4), 
		cex = 0.5,
		xlab = "max count cDNA", ylab = "max count CT",
		las = 1,
		col = alpha(c("blue","orange","magenta3","gray"), 0.4) [ cla [ cla$expression_bias_label != "empty", "expression_bias_label" ] ],
		sub = sprintf("pearson cor = %.2f", cz_cor),
		main = "CT ~ expression cell size correlation"
	)
	# colored by clicktag label
	plot(
		pmax(mat_cz_s1, mat_cz_s2) [ cla$expression_bias_label != "empty" ], 
		mct_o_max [ cla$expression_bias_label != "empty" ], 
		log = "xy", 
		# col = alpha("blue", 0.4), 
		cex = 0.5,
		xlab = "max count cDNA", ylab = "max count CT",
		las = 1,
		col = alpha(c("blue","orange","chartreuse3", "magenta3","gray50","gray"), 0.4) [ cla [ cla$expression_bias_label != "empty", "clicktag_label" ] ],
		sub = sprintf("pearson cor = %.2f", cz_cor),
		main = "CT ~ expression cell size correlation"
	)
	
	
	# #### UMAP metric learning ####
	# 
	# mct_o_clean = mct_o [ apply(mct_o, 1, function(x) any(!is.na(x))) & cla$expression_bias_label != "empty" , ]
	# lab_o_clean = cla [ rownames(mct_o_clean), c("expression_bias_label", "clicktag_label") ]
	# 
	# # lab_o_clean$label = as.character(lab_o_clean$clicktag_label)
	# # lab_o_clean$label [ lab_o_clean$expression_bias_label == "doublet" ] = "doublet"
	# # sort(table(lab_o_clean$label))
	# 
	# lab_o_clean$label = factor(paste(lab_o_clean[,1], lab_o_clean[,2]))
	# 
	# # supervised learning
	# mct_umap = uwot::umap(t(seu_o@assays$RNA@scale.data), y = factor(lab_o_clean$label))
	# # mct_umap = uwot::umap(t(seu_o@assays$RNA@scale.data))
	# rownames(mct_umap) = colnames(seu_o@assays$RNA@scale.data)
	# mct_umap = as.data.frame(mct_umap) [ rownames(cla) ,  ]
	# rownames(mct_umap) = rownames(mct_o)
	# 
	# # seurat clusters
	# plot(mct_umap, cex = 0.5, col = alpha(seu_clu_col,0.2), main = "UMAP learning, by seurat cluster (cells only)", xlab = "UMAP 1", ylab = "UMAP 2")
	# legend("topleft", cex = 0.6, col = color_palette_cols, legend = legend_s, pch = 1, bty = "n", ncol = 3, title = "Seurat clu")
	# 
	# # expression bias class
	# plot(mct_umap, cex = 0.5, col = alpha(c("blue","orange","magenta3","gray"), 0.4) [ cla$expression_bias_label ], main = "UMAP learning, by expr class (cells only)", xlab = "UMAP 1", ylab = "UMAP 2")
	# legend("topleft", sprintf("%s = %i", names(legend_t), legend_t), col = c("blue","orange","magenta3","gray"), pch = 1, bty = "n", title = "expr bias classification", cex = 0.6)
	# 
	# # clicktag class
	# plot(mct_umap, cex = 0.5, col = alpha(c("blue","orange","chartreuse3", "magenta3","gray50","gray"), 0.4) [ cla$clicktag_label ], main = "UMAP learning, by clicktag class (cells only)", xlab = "UMAP 1", ylab = "UMAP 2")
	# legend("topleft", sprintf("%s = %i", names(legend_c), legend_c), col = c("blue","orange","chartreuse3", "magenta3","gray50","gray"), pch = 1, bty = "n", title = "clicktag classification", cex = 0.6)
	# 
	# # original bcts
	# plot(mct_umap, cex = 0.5, col = alpha(bct_clu_col, 0.4), main = "UMAP learning, by ctbc (cells only)", xlab = "UMAP 1", ylab = "UMAP 2")
	# legend("topleft", levels(bct_fac), col = color_palette_cols, pch = 1, bty = "n", title = "bcs", cex = 0.6, ncol = 2)


	dev.off()



	#### Clicktag efficiency ####

	pdf(sprintf("data/dual_samples.%s.clicktag_efficiency.pdf",lib), height = 10, width = 10)
	par(mfrow=c(2,2))

	# distribution of counts per clicktag
	boxplot(
		mct_o [ cla$expression_bias_label != "empty" , ], las = 2, ylim = c(1,1000), log = "y", outline=FALSE,
		ylab = "counts",
		col = rainbow(n=max(bct_pairs), alpha = 1, v = 0.9)[bct_pairs],
		main = sprintf("Clicktag count distribution\n%s", lib))
	abline(h=median(mct_o [ cla$expression_bias_label != "empty" , ], na.rm = T), lty = 2)
	legend_p = sort(unique(bct_pairs))
	legend("topleft", legend = legend_p, fill = rainbow(n=max(bct_pairs), alpha = 1, v = 0.9)[legend_p], bty = "n", ncol = 2, cex = 0.6, title = "CT pair")

	# distribution of max counts per clicktag
	boxplot(
		10 ^ mct_t[ cla$expression_bias_label != "empty" , ]$first_val ~ mct_t[ cla$expression_bias_label != "empty" , ]$first_tag,
		las = 2, outline = FALSE,
		log = "y",
		ylab = "max counts",
		xlab = NULL,
		col = rainbow(n=max(bct_pairs), alpha = 1, v = 0.9)[bct_pairs],
		main = sprintf("Clicktag count distribution\n%s", lib))
	abline(h=median(10 ^ mct_t[ cla$expression_bias_label != "empty" , ]$relative_size_ft, na.rm = T), lty = 2)
	legend_p = sort(unique(bct_pairs))
	legend("topleft", legend = legend_p, fill = rainbow(n=max(bct_pairs), alpha = 1, v = 0.9)[legend_p], bty = "n", ncol = 2, cex = 0.6, title = "CT pair")

	# distribution of first-to-third cell sizes sizes, per most abundant clicktag
	boxplot(
		10 ^ mct_t[ cla$expression_bias_label != "empty" , ]$first_count ~ mct_t[ cla$expression_bias_label != "empty" , ]$first_tag,
		las = 2, outline = FALSE,
		log = "y",
		ylab = sprintf("first ctbc / third ctbc"),
		xlab = NULL,
		col = rainbow(n=max(bct_pairs), alpha = 1, v = 0.9)[bct_pairs],
		main = sprintf("Relative norm counts first-to-third clicktags\n%s", lib))
	abline(h=median(10 ^ mct_t[ cla$expression_bias_label != "empty" , ]$first_count, na.rm = T), lty = 2)
	legend("topleft", legend = legend_p, fill = rainbow(n=max(bct_pairs), alpha = 1, v = 0.9)[legend_p], bty = "n", ncol = 2, cex = 0.6, title = "CT pair")

	dev.off()
	
	
	#### Heatmaps ####
	
	pdf(sprintf("data/dual_samples.%s.heatmap.pdf",lib), height = 6, width = 6)
	
	# heatmap of clicktags
	mm = mct_o [real_cells,]
	mm = mm [ !is.na(mm[,1]), ]
	mm = mm [ order(apply(mm, 1, function(x) which.max(zoo::rollmean(x,2)) )), ]
	hm = plot_complex_heatmap(
		mat = mm,
		name = "count",
		color_min = 0, 
		color_max = quantile(mm, 0.9),
		cluster_col = FALSE,
		cluster_row = FALSE,
		title_col = sprintf("n = %i barcodes", ncol(mm)),
		title_row = sprintf("n = %i cells", nrow(mm)),
		name_row_show = FALSE)
	print(hm)
	
	# heatmap of clicktags
	mm = mct_n [real_cells,]
	mm = mm [ !is.na(mm[,1]), ]
	mm = mm [ order(apply(mm, 1, function(x) which.max(zoo::rollmean(x,2)) )), ]
	hm = plot_complex_heatmap(
		mat = mm,
		name = "norms",
		color_min = 0, 
		color_max = quantile(mm, 0.9), 
		cluster_col = FALSE,
		cluster_row = FALSE,
		title_col = sprintf("n = %i barcodes", ncol(mm)),
		title_row = sprintf("n = %i cells", nrow(mm)),
		name_row_show = FALSE)
	print(hm)
	
	# heatmap of classification, for expression bias and clicktag labels
	tac = table(cla$clicktag_label, cla$expression_bias_label)
	# tacf = tac / rowSums(tac) * 100
	tacf = t(t(tac) / colSums(tac) * 100)
	col_newmc = colorRampPalette(interpolate="l",c("gray95", "deepskyblue","dodgerblue3","midnightblue"))
	# cell counts
	pheatmap(
	    table_to_matrix(tac), 
		color = col_newmc(20), 
		breaks = seq(1,800,length.out = 20 + 1), 
		cellwidth = 10, cellheight = 10, na_col = "grey", 
		cluster_cols = FALSE, cluster_rows = FALSE, 
		number_format = "%i", fontsize = 5, number_color = "orange",
		border_color = "white", 
		display_numbers = TRUE, 
		main = sprintf("CT ~ expr bias\nn=%i", sum(tac)))
	# cell fraction
	pheatmap(
	    table_to_matrix(tacf), 
		color = col_newmc(20), 
		breaks = seq(0,100,length.out = 20 + 1), 
		cellwidth = 10, cellheight = 10, na_col = "grey", 
		cluster_cols = FALSE, cluster_rows = FALSE, 
		number_format = "%.1f", fontsize = 5, number_color = "orange",
		border_color = "white", 
		display_numbers = TRUE, 
		main = sprintf("CT ~ expr bias\nn=%i (col-wise)", sum(tac)))
	
	# heatmap of classification, for clicktags and seurat clusters
	tcs = table(cla$clicktag_label, seu_clu)
	tcsf = t(t(tcs) / colSums(tcs) * 100)
	# cell counts
	pheatmap(
	    table_to_matrix(tcs), 
		color = col_newmc(20), 
		breaks = seq(1,800,length.out = 20 + 1), 
		cellwidth = 10, cellheight = 10, na_col = "grey", 
		cluster_cols = FALSE, cluster_rows = FALSE, 
		number_format = "%i", fontsize = 5, number_color = "orange",
		border_color = "white", 
		display_numbers = TRUE, 
		main = sprintf("CT ~ Seurat\nn=%i", sum(tcs)))
	# cell fraction
	pheatmap(
	    table_to_matrix(tcsf), 
		color = col_newmc(20), 
		breaks = seq(0,100,length.out = 20 + 1), 
		cellwidth = 10, cellheight = 10, na_col = "grey", 
		cluster_cols = FALSE, cluster_rows = FALSE, 
		number_format = "%.1f", fontsize = 5, number_color = "orange",
		border_color = "white", 
		display_numbers = TRUE, 
		main = sprintf("CT ~ Seurat\nn=%i (col-wise)", sum(tcs)))
	
	# stop
	dev.off()
	
	
	#### Output classification ####
	
	# find which clicktag label best represents sps1 and sps2 (expression-defined)
	clicktag_sp1 = names(which.max(tac [ ,sp1 ]))
	clicktag_sp2 = names(which.max(tac [ ,sp2 ]))
	
	# find which seurat clusters should be blacklisted
	# blacklisted clusters will have less than X% from a single clicktag class
	# blacklisted clusters will have an excess of non-properly assigned clicktags (doublets, unclassified, etc)
	seurat_sum = colSums(tcsf [ c(clicktag_sp1, clicktag_sp2) ,])
	# names(seurat_black_sum) = paste("c", names(seurat_black_sum), sep = "")
	seurat_black_clu = names(which(seurat_sum < seurat_purity_frac * 100))
	
	# add seurat classification to table
	cla$seurat_cluster = seu_clu
	cla$seurat_blacklisted_cluster = seu_clu %in% seurat_black_clu
	
	# assign final label
	# cells where expression bias and clicktag coincide, and don't belong to a blacklisted seurat cluster
	cla$valid_cell = ((cla$expression_bias_label == sp1 & cla$clicktag_label == clicktag_sp1) | (cla$expression_bias_label == sp2 & cla$clicktag_label == clicktag_sp2)) & !cla$seurat_blacklisted_cluster
	cla$classification = "other"
	cla$classification [ cla$valid_cell ] = as.character(cla$expression_bias_label) [ cla$valid_cell ]
	
	# rebuild long cell names
	cla$long_cell_name = paste(lib,"_", rownames(cla), "-1", sep = "")
	
	# add relative sizes
	cla$ct_relative_size_ft = mct_t [ rownames(cla), "relative_size_ft" ]
	cla$ct_relative_size_fs = mct_t [ rownames(cla), "relative_size_fs" ]
	cla$ct_total_counts = rowSums(mct_o) [ rownames(cla) ]
	
	# write.table
	write.table(cla, sprintf("data/dual_samples.%s.class.tsv",lib), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
	
}


message("All done!")
