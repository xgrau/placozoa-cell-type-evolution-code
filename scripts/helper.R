# libraries
require("ComplexHeatmap")
require("circlize")
require("metacell")
require("tgconfig")
require("tgstat")
require("tglkmeans")
require("zoo")
require("data.table")
require("dplyr")
require("stringr")
require("ggplot2")
require("scales")
require("RColorBrewer")
require("rasterpdf")
require("vioplot")
require("MASS")

# deactivate complexheatmap warnings
ht_opt$message = FALSE

# from  github/tanaylab/metacell/R/utils.r
# rowFunction is a function that works on rows, like rowMeans
# # much faster than tapply
.row_stats_by_factor = function (data, fact, rowFunction = rowMeans) {
	u = as.character(sort(unique(fact)))
	fact[is.na(fact)] = FALSE
	n=length(u)
	centers = matrix(NA,dim(data)[1], n, dimnames = list(rownames(data), u))
	# much faster than tapply
	for (i in u) {
		if (sum(fact == i, na.rm=TRUE) > 1) {
			centers[,i] = rowFunction(data[,fact == i,drop=FALSE])
		} else {
			centers[,i] = data[,fact == i]
		}
	}
	return(centers)
}

# Helper function to save or return plots
#' 
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param EXPR expression that produces the plot (if multiple lines, enclose it in `{ ... }`)
#' 
#' @return plot (if `output_file` is `NULL`), otherwise `NULL`
#' 
plotting_function <- function(output_file=NULL, width, height, res=NA, EXP) {
	if (!is.null(output_file)) {
		extension <- stringr::str_extract(output_file,"(png|pdf)$")
		if ( is.na(extension) ) {
			extension = "pdf"
		}
		# open graphics device
		if (extension == "png") {
			png(output_file, height = height, width = width, res=res)
		} else if (extension == "pdf" ) {
			pdf(output_file, height = height, width = width, useDingbats=TRUE)
		}
	}
	EXP
	if (!is.null(output_file)) dev.off()
}


# create a transcript to gene dictionary from a GTF annotation file
dictionary_t2g = function(gtf_fn, vector_to_fix, t2g = TRUE, transcript_field = "transcript", transcript_id = "transcript_id", gene_id = "gene_id", return_elements_not_in_gtf = TRUE) {
	
	# import gtf
	gene_txgtf = rtracklayer::import(gtf_fn)
	
	if (t2g) {
		dic = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,gene_id]
		names(dic) = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,transcript_id]
	} else {
		dic = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,transcript_id]
		names(dic) = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,gene_id]
	}
	
	# return object
	
	out = dic [ vector_to_fix ]
	
	# return elements not in GTF dictionary, unaltered
	if (return_elements_not_in_gtf) {
		ixs_to_keep = is.na(out)
		out[ixs_to_keep] = vector_to_fix[ixs_to_keep]
		names(out[ixs_to_keep]) = out[ixs_to_keep]
	}
	
	# return
	out
	
}

# Heatmaps
plot_complex_heatmap = function(
		mat,
		name = "heatmap",
		color_mat = c("white","#d6e72e","#6fb600","#003f4d"),
		color_min = 0,
		color_max = 1,
		fontsize = 10,
		categories_col = NULL,
		categories_row = NULL,
		separate_col = FALSE,
		separate_row = FALSE,
		colors_col = NULL,
		colors_row = NULL,
		title_row = NULL,
		title_col = NULL,
		name_row_show = TRUE,
		name_col_show = TRUE,
		cluster_row = TRUE,
		cluster_col = TRUE,
		use_raster = TRUE,
		raster_quality = 1,
		show_legend_row = FALSE,
		show_legend_col = FALSE,
		both_sides_row = TRUE,
		both_sides_col = TRUE,
		cell_border = gpar(col = NA, lwd = 1, lty = 1),
		heatmap_border = gpar(col = NA, lwd = 1, lty = 1),
		do_dotplot = FALSE,
		dot_size_mat = NULL,
		dot_size_min = NULL,
		dot_size_max = NULL,
		cex_dotplot = 0.02,
		col_dotplot_border = NA
) {
	
	require("ComplexHeatmap")
	require("circlize")
	ht_opt$message = FALSE
	
	# color function
	col_fun = circlize::colorRamp2(seq(color_min, color_max, length.out = length(color_mat)), color_mat)
	
	
	# # vector of clusters (row/col annotation)
	# cluster_vector = as.character(mot_merge_d$cluster)
	# # categorical colors for clusters (get recycled)
	# catcol_vec = rep(catcol_lis, length.out = length(unique(cluster_vector)))
	# names(catcol_vec) = unique(cluster_vector)
	# cluster_colors = catcol_vec [ cluster_vector ]
	# names(cluster_colors) = cluster_colors
	# # quantitative colorscale for heatmap
	# catcol_fun = circlize::colorRamp2(1:length(catcol_vec), catcol_vec)
	# names(cluster_colors) = cluster_vector
	
	if (is.null(title_row)) {
		title_row = sprintf("n = %i rows", nrow(mat))
	}
	if (is.null(title_col)) {
		title_col = sprintf("n = %i columns", ncol(mat))
	}
	
	# left row annotations
	if (is.null(categories_row)) {
		categories_row = rownames(mat)
	}
	if (is.null(colors_row)) {
		colors_row = rep(NA, nrow(mat))
		left_annotation = NULL
	} else {
		if (is.null(names(colors_row))) { names(colors_row) = categories_row }
		left_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_row, name = title_row, col = list(c = colors_row), which = "row", show_legend = show_legend_row)
	}
	
	# top col annotations
	if (is.null(categories_col)) {
		categories_col = colnames(mat)
	}
	if (is.null(colors_col)) {
		colors_col = rep(NA, ncol(mat))
		top_annotation = NULL
	} else {
		if (is.null(names(colors_col))) { names(colors_col) = categories_col }
		top_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_col, name = title_col, col = list(c = colors_col), which = "column", show_legend = show_legend_col)
	}
	
	# add annotations to both sides of the rows or columns?
	if (both_sides_col) {
		bottom_annotation = top_annotation
	} else {
		bottom_annotation = NULL
	}
	if (both_sides_row) {
		right_annotation = left_annotation
	} else {
		right_annotation = NULL
	}
	
	# split columns and rows?
	if (separate_row) {
		split_row = categories_row
	} else {
		split_row = NULL
	}
	if (separate_col) {
		split_col = categories_col
	} else {
		split_col = NULL
	}
	
	# should this be a dot plot?
	if (do_dotplot) {
		if (is.null(dot_size_mat)) {
			dot_size_mat = mat
		}
		if (!is.null(dot_size_max) & !is.null(dot_size_min)) {
			# dot_size_mat [ dot_size_mat > dot_size_max ] = dot_size_max
			# dot_size_min [ dot_size_mat < dot_size_min ] = dot_size_min
			rangesize = function(x) { (x - dot_size_min) / (dot_size_max - dot_size_min) }
		} else {
			rangesize = function(x) { (x - min(x)) / (max(x) - min(x)) }
		}
		cell_fun_dotplot = function(j, i, x, y, width, height, fill) {
			grid.circle(
				x = x, y = y, 
				r = sqrt(rangesize(dot_size_mat)[i, j]) * cex_dotplot, 
				gp = gpar(col = col_dotplot_border, fill = col_fun(mat[i, j]))
			)
		}
		cell_border = gpar(type = "none")
		# cell_border = gpar(col = "gray", lwd = 1, lty = 1, fill = NULL)
	} else {
		cell_fun_dotplot = NULL
	}
	
	
	# plot
	hm = ComplexHeatmap::Heatmap(
		mat,
		name = name,
		cell_fun = cell_fun_dotplot,
		use_raster = use_raster,
		raster_quality = raster_quality,
		cluster_rows = cluster_row,
		cluster_columns = cluster_col,
		row_title = title_row,
		column_title = title_col,
		show_row_names = name_row_show,
		show_column_names = name_col_show,
		row_names_gp = gpar(fontsize = fontsize),
		column_names_gp = gpar(fontsize = fontsize),
		top_annotation = top_annotation,
		left_annotation = left_annotation,
		right_annotation = right_annotation,
		bottom_annotation = bottom_annotation,
		column_split = split_col,
		row_split = split_row,
		row_gap = unit(0.5, "mm"),
		column_gap = unit(0.5, "mm"),
		rect_gp = cell_border,
		border_gp = heatmap_border,
		col = col_fun)
	
	# return heatmap
	return(hm)
	
}




# inflection point of a knee plot
inflection_point <- function(counts_per_cell, min_value = 100) {
	
	counts_per_cell = counts_per_cell [!is.na(counts_per_cell)]
	# get log10 data
	dat = data.frame(
		umis = log10(sort(counts_per_cell, decreasing = TRUE)), 
		rank = log10(1:length(counts_per_cell)))
	
	# ignore cells below min value
	if (!is.null(min_value) & min_value != 0) {
		dat = dat[dat$umis > log10(min_value), ]
	}
	
	# rank difference
	dif_r = diff(dat$umis) / diff(dat$rank)
	min_ix = which.min(dif_r)
	
	# calculate inflection point
	inf = 10 ^ (dat$umis[min_ix])
	return(inf)
	
}


# valleys in a series of data
find_peaks = function (x, m = 10){
	shape <- diff(sign(diff(x, na.pad = FALSE)))
	pks <- sapply(which(shape < 0), FUN = function(i){
		z <- i - m + 1
		z <- ifelse(z > 0, z, 1)
		w <- i + m + 1
		w <- ifelse(w < length(x), w, length(x))
		if (all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
	})
	pks <- unlist(pks)
	pks
}

# wrapper to return peaks or valleys from an histogram
find_peaks_histogram = function(histogram, m = 10, valleys = FALSE) {
	
	x = histogram$counts
	
	if (valleys) {
		x = -x 
	}
	i = find_peaks(x, m = m)
	c = histogram$breaks [ i ]
	return(c)
	
}

# cell name cleaning function
# takes as input a vector of cellranger-style cell names and returns just the barcodes
cell_name_clean = function(v) {
	apply(stringr::str_split(v, "_", simplify = TRUE), 1, function(i) {
		i = c(i,"")
		ix = which(i == "") - 1
		c = i[ix][1]
		c = gsub("-\\d+$", "", c)
	}
	)
}


clean_og_pairs = function(og_pairs_fn, sp1, sp2, t2g = TRUE, t2g_sp1 = NULL, t2g_sp2 = NULL, t2g_sp1_field = "gene_id", t2g_sp2_field = "gene_id") {
	
	# load
	og_pairs = read.table(og_pairs_fn, header = FALSE, sep = "\t", col.names = c("g1", "g2"), stringsAsFactors = FALSE)
	og_pairs = og_pairs [ (grepl(sprintf("^%s_",sp1),og_pairs[,1]) & grepl(sprintf("^%s_",sp2),og_pairs[,2])) | (grepl(sprintf("^%s_",sp1),og_pairs[,2]) & grepl(sprintf("^%s_",sp2),og_pairs[,1])) , ]
	og_pairs = data.frame(
		sp1 = c(og_pairs [ grepl(sprintf("^%s_",sp1),og_pairs[,1]), 1 ] , og_pairs [ grepl(sprintf("^%s_",sp1),og_pairs[,2]), 2 ] ),
		sp2 = c(og_pairs [ grepl(sprintf("^%s_",sp2),og_pairs[,2]), 2 ] , og_pairs [ grepl(sprintf("^%s_",sp2),og_pairs[,1]), 1 ] )
	)
	
	# load gene to transcript dictionary
	if (!is.null(t2g_sp1)) {
		og_pairs$sp1 = dictionary_t2g(gtf_fn = t2g_sp1, vector_to_fix = og_pairs$sp1, gene_id = t2g_sp1_field)
	}
	if (!is.null(t2g_sp2)) {
		og_pairs$sp2 = dictionary_t2g(gtf_fn = t2g_sp2, vector_to_fix = og_pairs$sp2, gene_id = t2g_sp2_field)
	}
	og_pairs = og_pairs [ !is.na(og_pairs$sp1) & !is.na(og_pairs$sp2) , ]
	
	return(og_pairs)
	
}


# subsample matrix
matrix_subsample = function(mat, subsample_rate = 0.1, subsample_rate_row = NULL, subsample_rate_col = NULL) {
	
	# subsample rates
	if (is.null(subsample_rate_row)) { subsample_rate_row = subsample_rate }
	if (is.null(subsample_rate_col)) { subsample_rate_col = subsample_rate }
	
	# num bins and widths
	nbin_r = nrow(mat) * subsample_rate_row
	nbin_c = ncol(mat) * subsample_rate_col
	wbin_r = round(nrow(mat) / nbin_r)
	wbin_c = round(ncol(mat) / nbin_c)
	
	# init subsampled matrix
	smat = matrix(nrow = nbin_r, ncol = nbin_c)
	
	# log
	message(sprintf("subsample matrix at %.2f x %.2f rate: %i x %i to %i x %i", subsample_rate_row, subsample_rate_col, nrow(mat), ncol(mat), nrow(smat), ncol(smat)))
	
	# loop
	for (rs in 0:(nbin_r - 1)) {
		ro_s = (rs * wbin_r) + 1
		ro_e = rs * wbin_r + wbin_r
		for (cs in 0:(nbin_c - 1)) {
			co_s = (cs * wbin_c) + 1
			co_e = cs * wbin_c + wbin_c
			smat_i = mat [ ro_s : ro_e , co_s : co_e ]
			smat[rs + 1, cs + 1]   = mean(smat_i)
		}
	}
	
	# out
	return(smat)
	
}

# Helper function to round values without changing total sum
round_smart <- function(x, digits = 0) {
	up <- 10 ^ digits
	x <- x * up
	y <- floor(x)
	indices <- tail(order(x-y), round(sum(x)) - sum(y))
	y[indices] <- y[indices] + 1
	y / up
}

# jaccard from vectors
jaccard_index = function(v1, v2) {
	
	l_int = length(intersect(v1, v2))
	l_uni = length(union(v1, v2))
	return(l_int / l_uni)
	
}

#' Calculate confusion matrix - needed for plotting function
#' 
mc_compute_norm_confu_matrix <- function(
	mc_id, graph_id, max_deg=NULL
){
	cgraph <- scdb_cgraph(graph_id)
	if(is.null(max_deg))	max_deg <- nrow(cgraph@edges)
	confu <- mcell_mc_confusion_mat(mc_id, graph_id, max_deg,ignore_mismatch=T)
	r_confu <- rowSums(confu)
	c_confu <- colSums(confu) 
	norm <- r_confu %*% t(c_confu)
	confu_n <- confu/norm
	confu_nodiag <- confu_n
	diag(confu_nodiag) <- 0
	confu_n <- pmin(confu_n, max(confu_nodiag))
	confu_n <- pmin(confu_n, quantile(confu_n, 1 - 3/nrow(confu_n)))
	return(confu_n)
}

#' Cluster confusion matrix - needed for plotting function
#' 
mc_confusion_clustering <- function(
	confu_n, clust_method="average"
){
	epsilon <- quantile(confu_n[confu_n != 0], 0.02)
	hc <- hclust(as.dist(-log10(epsilon + confu_n)), clust_method)
	hc$height <- round(hc$height, 6) 
	return(hc)
}


#' Plot metacell size and UMI counts
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param mc_counts UMI counts per metacell (`dataframe` class, genes are rows)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param return_table logical; whether to return a data.frame with per-metacell size stats
#'
scp_plot_mc_size_counts = function(mc_object, mat_object, mc_color="blue", output_file = NULL, width = 24, height = 6, res = NA, return_table = FALSE) {
	
	# UMI counts per metacell
	mc_counts = sca_mc_gene_counts(mc_object,mat_object,0)
	# num cells per metacell
	mc_n_cells = table(mc_object@mc)
	# median cell size
	sc_counts = colSums(as.matrix(mat_object@mat[,names(mc_object@mc)]))
	mc_median_cell_size = tapply(sc_counts, mc_object@mc, median)
	mc_mean_cell_size = tapply(sc_counts, mc_object@mc, mean)
	
	# keep metacell order from mc object
	mc_counts = mc_counts [ , colnames(mc_object@mc_fp) ]
	mc_n_cells = mc_n_cells [ colnames(mc_object@mc_fp) ]
	mc_median_cell_size = mc_median_cell_size [ colnames(mc_object@mc_fp) ]
	mc_object_mc_factor = factor(mc_object@mc, levels = colnames(mc_object@mc_fp))
	
	plotting_function(output_file, width, height, res, EXP = {
		
		# metacell size (num cells)
		barplot(
			mc_n_cells,
			col = mc_color, 
			border="white", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space= 0,
			ylab = "# cells" ,
			main = "# cells per metacell")
		
		boxplot(
			sc_counts + 1 ~ mc_object_mc_factor,
			col = mc_color, 
			border= "black", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space = 0, 
			ylab = "# UMI" ,
			log = "y",
			main = "# UMI per metacell")
		abline(h=0, lty=2, col="gray")
		
		# median umi per metacell
		barplot(
			mc_median_cell_size,
			col = mc_color, 
			border="white", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space=0, 
			ylab = "# UMI" ,
			main = "Median # UMI per metacell")
		
		# total umi per metacell
		barplot(
			colSums(mc_counts),
			col = mc_color, 
			border="white", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space=0, 
			ylab = "# UMI" ,
			main = "Total # UMI per metacell")
		
		# distrubution of foldchange values per metacell
		vioplot::vioplot(
			log2(apply(mc@mc_fp, 2, function(x) x[x > 0] )),
			col = mc_color, 
			border= "black", 
			xlim = c( 1, ncol(mc_counts) + 60 ),
			las = 2,
			space = 0, 
			ylab = "log2 fold change" ,
			main = "fold changes per metacell")
		abline(h=log2(c(1,1.5,2)), lty=2, col="gray")
		
	})
	
	# return table?
	if (return_table) {
		
		sta = data.frame(
			metacell = colnames(mc@mc_fp),
			num_cells = as.vector(mc_n_cells),
			tot_umis = colSums(mc_counts),
			cell_size_mean   = tapply(sc_counts, mc_object@mc, mean),
			cell_size_Q1     = tapply(sc_counts, mc_object@mc, function(i) quantile(i, 0.25)),
			cell_size_Q2     = mc_median_cell_size,
			cell_size_Q3     = tapply(sc_counts, mc_object@mc, function(i) quantile(i, 0.75)),
			log2_fc_mean     = log2(apply(mc@mc_fp, 2, function(i) mean(i[i > 1]))),
			log2_fc_Q1       = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.25))),
			log2_fc_Q2       = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.50))),
			log2_fc_Q3       = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.75))),
			log2_fc_Qp95     = log2(apply(mc@mc_fp, 2, function(i) quantile(i[i > 1], 0.95)))
			
		)
		
	}
	
}


#' Plot single cell/metacell 2D projection with single gene expression
#'
#' @param mc2d 2D projection from `mc2d` object
#' @param mc loaded metacell object (`gMCCov` class)
#' @param mat loaded single cell matrix object (`tgScMat` class)
#' @param gene_id string, gene id (rownames in `mat@mat` object)
#' @param sc_vector vector of single cell-level values to map. If NULL (default), the # of UMIs is taken from the `mat` object.
#' @param sc_scale indicate wether cell-level values are  unidirectional (`unidir`) or bidirectional (`bidir`); defaults to `unidir`
#' @param sc_transform a function to transform sc vector (default is NULL, for no transformation)
#' @param sc_min min value of sc vector (lower values are raised to this threshold); should be in the transformed scale if transformation is used
#' @param sc_max max value of sc vector (higher values are lowered to this threshold); should be in the transformed scale if transformation is used
#' @param sc_max_quant if sc_max is NULL, use this quantile to find the upper threshold (default is 0.75)
#' @param sc_label legend title label for the sc vector. Defaults to "cell UMI"
#' @param sc_zero_color a color to be used for zero values in the sc vector (useful to distinguish zero values from the rest); default to light gray.
#' @param do_umifrac_sc boolean; if set to TRUE, UMI fraction per 10k is calculated, and `sc_label` is set to "UMI/10k"
#' @param mc_vector vector of metacell-level values to map. If NULL (default), footprint values from `mc@mc_fp` are used and log2-transformed.
#' @param mc_scale indicate wether mc-level values are  unidirectional ("unidir") or bidirectional ("bidir"); defaults to "unidir"
#' @param mc_transform a function to transform mc vector (default is `log` transformation)
#' @param mc_min min value of mc vector (lower values are raised to this threshold); should be in the transformed scale if transformation is used
#' @param mc_max max value of mc vector (higher values are lowered to this threshold); should be in the transformed scale if transformation is used
#' @param mc_max_quant if mc_max is NULL, use this quantile to find the upper threshold (default is 0.75)
#' @param mc_label legend title label for the mc vector. Defaults to "mc log(fp)"
#' @param mc_zero_color a color to be used for zero values in the mc vector (useful to distinguish zero values from the rest); default to NULL
#' @param unidir_color_scale vector of colors that will be used for unidirectional color scales
#' @param bidir_color_scale vector of colors that will be used for bidirectional color scales
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param plot_mcs whether to plot metacells (default FALSE)
#' @param plot_mc_name whether to plot metacell names (default FALSE)
#' @param plot_edges whether to plot edges between metacells (default FALSE)
#' @param plot_cells_as_2d_density instead of plotting expression in single cells, create a 2d density plot from the coordinates of single cells (can work well with UMI counts/fraction data; default FALSE).
#' @param cex_mc size of metacell points (default 3)
#' @param cex_sc size of single cell points (default 0.75)
#' @param alpha_mc alpha value for metacell points (default 0.9)
#' @param alpha_sc alpha value for single cell points (default 0.7)
#' @param alpha_sc_2d alpha value for single cell 2d density (default 1)
#' @param do_legend_sc,do_legend_mc whether to add legends with expression scale for metacells and single cells (default TRUE)
#' @param do_axes whether to plot axes for the 2D projection (default FALSE)
#'
scp_plot_sc_2d_gene_exp = function(
	mc2d,
	mc,
	mat,
	gene_id,
	sc_vector = NULL,
	sc_scale = "unidir",
	sc_transform = NULL,
	sc_min = 0,
	sc_max = NULL,
	sc_max_quant = 0.75,
	sc_label = "cell UMI",
	sc_zero_color = "lightblue1",
	do_umifrac_sc = FALSE,
	mc_vector = NULL,
	mc_scale = "unidir",
	mc_transform = log2,
	mc_min = 0,
	mc_max = 2,
	mc_max_quant = 0.75,
	mc_label = "mc log2(fp)",
	mc_zero_color = NULL,
	unidir_color_scale = c("gray90","orange","orangered2","#520c52"),
	bidir_color_scale = c("midnightblue","dodgerblue3","deepskyblue","#b8e0ed","gray90","#eccac0","#ff8d36","#f34312","#8e0631"),
	plot_mcs=FALSE,
	plot_mc_name=FALSE,
	plot_edges=FALSE,
	plot_cells_as_2d_density=FALSE,
	title=NULL,
	width=12,
	height=12,
	res=NA,
	output_file=NULL,
	cex_mc=3,
	cex_sc=0.75,
	alpha_mc=0.9,
	alpha_sc=0.7,
	alpha_sc_2d=1,
	do_legend_sc=TRUE,
	do_legend_mc=TRUE,
	do_axes=FALSE) {

	# get expression vectors (if not already provided)
	if (is.null(sc_vector)) {
		sc_vector = mat@mat[gene_id,]
	}
	if (is.null(mc_vector)) {
		mc_vector = mc@mc_fp[gene_id,]
	}


	# apply transformations
	if (do_umifrac_sc) {
		sc_vector = sc_vector / apply(mat@mat, 2, function(x) sum(x, nar.rm = TRUE)) * 10000
		sc_label = "UMI/10k"
	}
	if (!is.null(sc_transform)) {
		sc_vector = sc_transform(sc_vector)
	}
	if (!is.null(mc_transform)) {
		mc_vector = mc_transform(mc_vector)
	}


	# if sc_max or mc_max values are NULL, get them from a high quantile in the mc_vector
	if (is.null(sc_max)) {
		sc_max = round(quantile(sc_vector[sc_vector > 0], sc_max_quant, na.rm=TRUE))
	}
	if (is.null(mc_max)) {
		mc_max = round(quantile(mc_vector[mc_vector > 0], mc_max_quant, na.rm=TRUE))
	}

	# apply min/max to sc vector
	sc_vector [ sc_vector < sc_min ] = sc_min
	sc_vector [ sc_vector > sc_max ] = sc_max
	# apply min/max to mc vector
	mc_vector [ mc_vector < mc_min ] = mc_min
	mc_vector [ mc_vector > mc_max ] = mc_max

	# create color palettes for single cells
	if (sc_scale == "unidir") {
		sc_color_fun = scales::col_numeric(palette=unidir_color_scale, domain=c(sc_min, sc_max))
	} else if (sc_scale == "bidir") {
		sc_color_fun = scales::col_numeric(palette=bidir_color_scale, domain=c(sc_min, sc_max))
	} else {
		message("`sc_scale` should be either `unidir` or `bidir`")
	}
	sc_color = sc_color_fun(sc_vector)
	sc_color_labels = seq(from=sc_min, to=sc_max, length.out=9)
	sc_color_legend = sc_color_fun(sc_color_labels)
	sc_color_labels = sprintf("%.2f",sc_color_labels)
	sc_color_labels [ length(sc_color_labels) ] = paste(">=", sc_color_labels [ length(sc_color_labels) ], sep="")
	sc_color_labels [ 1 ] = paste("<=", sc_color_labels [ 1 ], sep="")
	# zero color
	if (!is.null(sc_zero_color)) {
		sc_color [ sc_vector == 0 ] = sc_zero_color
		sc_color_legend [ sc_color_labels == 0 ] = sc_zero_color
	}

	# color palettes for metacells
	if (mc_scale == "unidir") {
		mc_color_fun = scales::col_numeric(palette=unidir_color_scale, domain=c(mc_min, mc_max))
	} else if (mc_scale == "bidir") {
		mc_color_fun = scales::col_numeric(palette=bidir_color_scale, domain=c(mc_min, mc_max))
	} else {
		message("`sc_scale` should be either `unidir` or `bidir`")
	}
	mc_color = mc_color_fun(mc_vector)
	mc_color_labels = seq(from=mc_min, to=mc_max, length.out=9)
	mc_color_legend = mc_color_fun(mc_color_labels)
	mc_color_labels = sprintf("%.2f",mc_color_labels)
	mc_color_labels [ length(mc_color_labels) ] = paste(">=", mc_color_labels [ length(mc_color_labels) ], sep="")
	mc_color_labels [ 1 ] = paste("<=", mc_color_labels [ 1 ], sep="")
	# zero color
	if (!is.null(mc_zero_color)) {
		mc_color [ mc_vector == 0 ] = mc_zero_color
		mc_color_legend [ mc_color_labels == 0 ] = mc_zero_color
	}

	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		# determine plot max/min
		if (plot_mcs) {
			xlim=c(min(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE), max(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE), max(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE))
		} else {
			xlim=c(min(mc2d@sc_x, na.rm = TRUE), max(mc2d@sc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, na.rm = TRUE), max(mc2d@sc_y, na.rm = TRUE))
		}
		xlim[1] = xlim[1] - 0.2 * abs(xlim[1] - xlim[2])
		xlim[2] = xlim[2] + 0.2 * abs(xlim[1] - xlim[2])
		ylim[1] = ylim[1] - 0.2 * abs(ylim[1] - ylim[2])
		ylim[2] = ylim[2] + 0.2 * abs(ylim[1] - ylim[2])
		
		# get single cell coordinates
		mc2d_sc_x = mc2d@sc_x [ names(mc2d@sc_x) %in% names(sc_vector) ]
		mc2d_sc_y = mc2d@sc_y [ names(mc2d@sc_y) %in% names(sc_vector) ]
		sc_vector = sc_vector [ names(sc_vector) %in% names(mc2d_sc_x) ]
		
		# plot individual cells
		# by default, plot individual cells
		if (! plot_cells_as_2d_density ) {
			plot(
				mc2d_sc_x,
				mc2d_sc_y,
				pch = 19,lwd = 0,
				cex = cex_sc,
				col = alpha(sc_color, alpha_sc),
				xlim = xlim,
				ylim = ylim,
				xlab = NA, ylab = NA, axes=do_axes
			)
		
		# if else, plot 2d density (works with UMI counts or UMI frac, probably with anythin that's unidirectional too)
		} else {
			
			# create fake mc2d positions where each cell appears in the 
			# same positionas many times as its UMI counts
			den_mc2d_x = unlist(plyr::alply(1:length(sc_vector), 1, function(i) {
				if (sc_vector[i] > sc_min) { rep(mc2d_sc_x[i], sc_vector[i]) }
			}))
			den_mc2d_y = unlist(plyr::alply(1:length(sc_vector), 1, function(i) {
				if (sc_vector[i] > sc_min) { rep(mc2d_sc_y[i], sc_vector[i]) }
			}))
			
			# drop NA?
			den_mc2d_x = den_mc2d_x [ !is.na(den_mc2d_x) ]
			den_mc2d_y = den_mc2d_y [ !is.na(den_mc2d_y) ]
			
			# add a bit of jitter
			den_mc2d_x = jitter(den_mc2d_x, factor = 1)
			den_mc2d_y = jitter(den_mc2d_y, factor = 1)
			
			# create color scale
			den_mc2d_colvec = unidir_color_scale
			den_mc2d_colvec[1] = NA
			den_mc2d_colfun = colorRampPalette(den_mc2d_colvec)
			den_mc2d_col = den_mc2d_colfun(40)
			
			# create interploated 2d density plot
			den_mc2d_k = MASS::kde2d(den_mc2d_x, den_mc2d_y, n = 200, lims = c(xlim, ylim))
			
			# plot
			image(den_mc2d_k, col= alpha(den_mc2d_col, alpha_sc_2d), xlim = xlim, ylim = ylim, axes = do_axes, add = FALSE, useRaster = TRUE)
			points(
				mc2d@sc_x,
				mc2d@sc_y,
				pch = 1, lwd = 0.5,
				cex = cex_sc,
				col = alpha("black", alpha_sc),
				xlim = xlim,
				ylim = ylim,
				xlab = NA, ylab = NA, axes=do_axes
			)
			
		}
		
		# plot legend
		if (do_legend_sc) {
			legend("topleft",fill=sc_color_legend, legend=sc_color_labels, cex=0.4, title=sc_label)
		}
		
		# plot edges between metacells?
		if (plot_edges) {
			fr = mc2d@graph$mc1
			to = mc2d@graph$mc2
			segments(
				x0=mc2d@mc_x[fr],
				y0=mc2d@mc_y[fr],
				x1=mc2d@mc_x[to],
				y1=mc2d@mc_y[to],
				col="gray70")
		}
		
		# plot metacells?
		if (plot_mcs) {
		points(
			mc2d@mc_x,
			mc2d@mc_y,
			pch=19,lwd=0.5,
			cex=cex_mc,
			col=alpha(mc_color,alpha_mc))
			if (do_legend_mc) {
				legend("topright",fill=mc_color_legend, legend=mc_color_labels, cex=0.4, title=mc_label)
			}
		}
		
		# plot metacell ids?
		if (plot_mc_name) {
			text(
				mc2d@mc_x,
				mc2d@mc_y,
				#cex=cex_mc,
				labels=names(mc2d@mc_x))
		}
		
		# plot title
		if (is.null(title)) {
			title = gene_id
		}
		title(
			main = title,
			sub = sprintf(
				"2D projection\nn = %i cells | n = %i metacells", 
				length(mc2d@sc_x),
				length(mc2d@mc_x)
			)
		)
		
	}
	)

}


#' Plot single cell/metacell 2D projection
#'
#' @param mc2d 2D projection from `mc2d` object
#' @param mc loaded metacell object (`gMCCov` class)
#' @param mc_colors optional named vector of metacell colors; if `NULL` (default), `mc@colors` are used
#' @param cell_colors optional named vector of individual cell colors; if `NULL` (default), metacell color is used
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' 
#' @return data.frame with metacells' cell_type annotations. Also saves the plot to disk.
#'
scp_plot_sc_2d = function(
	mc2d,
	mc,
	mc_colors=NULL,
	cell_colors=NULL,
	plot_edges=FALSE,
	plot_mcs=FALSE,
	plot_mc_name=FALSE,
	width=12,
	height=12,
	res=NA,
	output_file=NULL,
	cex_mc=3,
	cex_sc=0.75,
	alpha_mc=0.8,
	alpha_sc=0.4,
	do_axes=FALSE) {
	
	# get colors of metacells
	if (is.null(mc_colors)) mc_colors = mc@colors
	if (is.null(names(mc_colors))) names(mc_colors) <- colnames(mc@mc_fp)
	# get colors of individual cells
	if (is.null(cell_colors)) cell_colors = mc_colors[ as.character(mc@mc[names(mc2d@sc_x)]) ]
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		# determine plot max/min
		if (plot_mcs) {
			
			xlim=c(min(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE), max(mc2d@sc_x, mc2d@mc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE), max(mc2d@sc_y, mc2d@mc_y, na.rm = TRUE))
			
		} else {
			
			xlim=c(min(mc2d@sc_x, na.rm = TRUE), max(mc2d@sc_x, na.rm = TRUE))
			ylim=c(min(mc2d@sc_y, na.rm = TRUE), max(mc2d@sc_y, na.rm = TRUE))
			
		}
		
		# plot individual cells
		plot(
			mc2d@sc_x,
			mc2d@sc_y,
			pch = 19,lwd = 0,
			cex = cex_sc,
			col = alpha(cell_colors,alpha_sc),
			xlim = xlim,
			ylim = ylim,
			xlab = NA, ylab = NA, axes=do_axes
		)
		
		# plot edges between metacells?
		if (plot_edges) {
			fr = as.character(mc2d@graph$mc1)
			to = as.character(mc2d@graph$mc2)
			segments(
				x0=mc2d@mc_x[fr],
				y0=mc2d@mc_y[fr],
				x1=mc2d@mc_x[to],
				y1=mc2d@mc_y[to],
				col="gray70")
		}
		
		# plot metacells?
		if (plot_mcs) {
			points(
				mc2d@mc_x,
				mc2d@mc_y,
				pch=19,lwd=0.5,
				cex=cex_mc,
				col=alpha(mc_colors,alpha_mc))
		}
		
		# plot metacell ids?
		if (plot_mc_name) {
			text(
				mc2d@mc_x,
				mc2d@mc_y,
				#cex=cex_mc,
				labels=names(mc2d@mc_x))
		}
		title(
			sub = sprintf(
				"2D projection\nn = %i cells | n = %i metacells", 
				length(mc2d@sc_x),
				length(mc2d@mc_x)
			)
		)
		
	})
}

#' Prepare heatmap of gene expression: select marker genes & create annotations
#'
#' Prepares heatmap of gene expression fold change for metacells and single cells
#' (no plotting done, returns list with selected markers and prepared annotations),
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param black_list character, blacklisted genes
#' @param sub_list_mc
#' @param gene_list character, list of genes to plot, if NULL (default) ...
#' @param order_genes logical, whether to cluster genes (default: TRUE)
#' @param gene_annot_file charcter, file path to the file containig gene annotations,
#'    it should have three tab separated columns containing gene ID, pfam architecture
#'    and any additional annotation in the last column
#' @param annot_header logical, gene annotation file has column names?
#' @param gene_font_size numeric, size of the gene names plotted as rownames
#'    (default: 4)
#' @param clust_ord character, metacells in the order in which they should be
#'    plotted; if cluster order is not specified (default: NULL), it is
#'    determined by hierarchical clustering
#' @param per_clust_genes integer, how many genes per cluster to aim to show in the heatmap
#'    (default: 20)
#' @param gene_min_fold numeric, minimum fold change for a gene to be considered for plotting
#'    (default: 2)
#' @param transverality_N integer, number of metacells in which a gene can be highly expressed (>transverality_thr, by default >1.4)
#'    to be considered for plotting, by default this is the total number of metacells
#' @param transverality_thr integer, expression threshold for the transverality_N filter (>1.4)
#' @param transv_excluded_mc character, metacells to be excluded in transversality calculation
#'    (default: NULL)
#' @param output_file optionally, a path to RDS file to which the function output will be saved
#'
scp_plot_cmod_markers_select <- function(
	mc_object,
	black_list = c(),
	sub_list_mc = NULL,
	gene_list = NULL,
	order_genes = TRUE,
	gene_annot_file = NULL,
	annot_header = FALSE,
	gene_font_size = 4,
	clust_ord = NULL,
	per_clust_genes = 20,
	gene_min_fold = 2,
	n_marker_order_rollmean = 1,
	transverality_N = ncol(mc_object@mc_fp),
	transverality_thr = 1.4,
	transv_excluded_mc = NULL,
	output_file=NULL
) {
	
	# load gene annotations
	if (!is.null(gene_annot_file) & "character" %in% class(gene_annot_file)) {
		annot = read.table(gene_annot_file, header=annot_header, sep="\t", fill=TRUE, quote="", row.names=1)
	} else if (!is.null(gene_annot_file) & "data.frame" %in% class(gene_annot_file)) {
		annot = gene_annot_file
	}
	
	
	# expression matrix
	if (is.null(sub_list_mc)) {
		niche_geomean_n= mc_object@mc_fp
	} else {
		niche_geomean_n= mc_object@mc_fp[,sub_list_mc]
		clust_ord=sub_list_mc
	}
	
	# exclude genes with fc < gene_min_fold
	genes=unique(
		as.vector(unlist( apply(niche_geomean_n, 2, function(x) names(head(sort(-x[x > gene_min_fold]), n=per_clust_genes ))) ) )
	)
	
	# select genes for plotting
	# if (!is.null(black_list)) 
	#   black_list <- vector("character")
	if (is.null(gene_list)){
		# exclude blacklisted genes
		genes=setdiff(genes, black_list)
		message(sprintf("Excluded %s blacklisted genes", length(black_list)))
		# exclude genes with transversality > transverality_N
		transversal_genes=names(which(
			apply(
				niche_geomean_n[,setdiff(as.character(colnames(niche_geomean_n)),transv_excluded_mc)],
				1,
				function(x) sort(x,decreasing=TRUE)[transverality_N] > transverality_thr
			)
		))
		genes=setdiff(genes, transversal_genes)
	} else {
		# plot only genes in gene list, if it is specified
		genes <- gene_list[gene_list %in% genes]
	}
	genes=genes[genes %in% rownames(niche_geomean_n)]
	message("Will use ",length(genes)," genes")
	
	mat_niche <- niche_geomean_n[genes,]
	
	# if cluster order is not specified, do hierarchical clustering
	if (is.null(clust_ord)) {
		message("Recomputing cell ord")
		hc1 = hclust(dist(cor(mat_niche,method="pearson")), "ward.D2")
		clust_ord = as.character(hc1$order)
		scr_tmp_niche_order <- as.character(hc1$order)
	}
	
	# if gene order is TRUE, order genes
	if (order_genes){
		message("Ordering genes")
		gene_ord = order(apply(mat_niche[,as.character(clust_ord)],1,function(x) which.max(rollmean(x,n_marker_order_rollmean))))
	} else {
		gene_ord= 1:nrow(mat_niche)
	}
	gene_ord <- rev(gene_ord)
	
	# gene labels
	if (!is.null(gene_annot_file)) {
		
		gene_labels_0 <- genes[gene_ord]
		
		message("Genes: ", head(gene_labels_0), "...")
		
		gene_labels_1 <- as.character(annot[genes[gene_ord],2])
		message("Gene labels: ", head(gene_labels_0), "...")
		bad_labels <- gene_labels_1 %in% c("","-"," ") | is.na(gene_labels_1)
		message(sum(bad_labels), " bad gene labels")
		gene_labels_1[bad_labels] <- gene_labels_0[bad_labels]
		# OMIT GENE ANNOTATION SHORTENING: truncation + padding done in the actual plotting functions
		# long_labels <- nchar(gene_labels_1)>gene_chr_limit
		# message(sum(long_labels), " long gene labels")
		# gene_labels_1[long_labels] <- paste0(substr(gene_labels_1[long_labels],1,gene_chr_limit-3),"...")
		names(gene_labels_1) <- genes[gene_ord]
		
		gene_labels_3 <- as.character(annot[genes[gene_ord],1])
		bad_labels <- gene_labels_3 %in% c("","-"," ") | is.na(gene_labels_3)
		gene_labels_3[bad_labels] <- genes[gene_ord][bad_labels]
		# long_labels <- nchar(gene_labels_3)>gene_chr_limit
		# gene_labels_3[long_labels] <- paste0(substr(gene_labels_3[long_labels],1,gene_chr_limit-3),"...")
		gene_labels_2 <- ifelse(gene_labels_0 == gene_labels_3, gene_labels_0, paste(gene_labels_0,gene_labels_3, sep=" "))
		names(gene_labels_2) <- genes[gene_ord]
		
	} else {
		
		gene_labels_0 <- genes[gene_ord]
		gene_labels_1 <- genes[gene_ord]
		gene_labels_2 <- genes[gene_ord]
		
	}
	
	# return objects necessary for plotting
	marker_data_list = list(
		genes = genes,
		gene_ord = gene_ord,
		clust_ord = clust_ord,
		niche_geomean_n = niche_geomean_n,
		gene_labels_1 = gene_labels_1,
		gene_labels_2 = gene_labels_2
	)
	
	if (!is.null(output_file)) saveRDS(marker_data_list, output_file)
	
	return(marker_data_list)
	
}



#' Plot heatmap of gene expression for metacells
#'
#' Plots heatmap of gene expression fold change for metacells. Requires selecting markers and ordering first with `scp_plot_cmod_markers_select`
#'
#' @param marker_data_list character, list object produced by the heatmap pre-processing function `scp_plot_cmod_markers_select`
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param show_gene_names logical, show gene names as heatmap row names
#'    (default: FALSE)
#' @param highlight_genes vector of gene names to highlight in the heatmap (in a different color).
#' @param highlight_color highlight color for genes in `highlight_genes`; default is "blue"
#' @param omit_unhighlighted if TRUE, do not plot any other gene name other than the ones in `highlight_genes` (default: FALSE)
#' @param gene_font_size numeric, size of the gene names plotted as rownames
#'    (default: 4)
#' @param gene_chr_limit numeric, limit gene annotations to given number of charaters
#' @param clust_row data.frame with gene colour annotations in columns,
#'    rownames should be `marker_data_list$gids`; 
#'    if `NULL` (default), row colour annotation bar is not printed
#' @param clust_row_color named list of color mappings for annotations in `clust_row`; 
#'    names should be `colnames(clust_row)`
#' @param clust_col either a character vector or data.frame with cell colour annotations;
#'    if `clust_col` is a vector, it should specify colors to assign to metacells in the order given by `gene_list$clust_ord`,
#'    and if it is named, the names should be metacell names;
#'    if `clust_col` is a data.frame, rownames should be metacell names, annotations (not colors) should be in columns, 
#'    and color mappings should be given by `clust_col_color` argument;
#'    if `NULL` (default), cluster colour annotation bar is not printed
#' @param clust_col_color named list of color mappings for annotations in `clust_col`; 
#'    names should be `colnames(clust_col)`
#' @param clust_bars numeric, optional values to be ploted as bars annotation on top of
#'    heatmap columns, length should be same as `length(mc@mc_fp)` (default:  NULL)
#' @param clust_bars_color` character, color for barplots, can be either single color
#'    or a (named) vector (names should be `gene_list$clust_ord`)
#' @param clust_anno_size unit, height of the column annotation bar (default: `unit(1,"mm")`)
#' @param show_mc_names logical, show metacell names as heatmap column names (default: TRUE) 
#' @param use_raster logical, whether to rasterise heatmap bodies (default: TRUE) 
#' @param mc_font_size numeric, size of the metacell names (default: 4)
#' @param heatmap_colors vector of colors to use for heatmap coloring function (low to high; default: `c("white","gray99","orange","orangered2","#520c52")`)
#' @param max_expression_fc,min_expression_fc numeric, max and minimum expression values to scale metacell heatmap coloring to (default: 5 and 0). It's applied AFTER any transformation of the data.
#' @param add_expression_constant numeric, a psuedocount to add to the expresion matrix (default 1; useful if `use_log2` = TRUE). It's applied BEFORE transformation of the data.
#' @param transformation_fun, the name of a function used to transform the data (default is `log2`; to omit, set to NULL)
#'
scp_plot_cmod_markers_mc <- function(
	marker_data_list,
	mat_label = "FC",
	output_file = NULL,
	height = 10,
	width = 5,
	res = NA,
	show_heatmap_legend = TRUE,
	show_gene_names = FALSE,
	highlight_genes = NULL,
	highlight_color = "blue",
	omit_unhighlighted = FALSE,
	gene_font_size = 5,
	clust_row = NULL, clust_row_color=NULL,
	clust_col = NULL, clust_col_color=NULL,
	clust_bars = NULL, clust_bars_color="grey",
	clust_anno_size = unit(4,"mm"), clust_anno_gap = unit(1,"mm"),
	show_mc_names = TRUE, mc_font_size = 5,
	heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
	max_expression_fc = 5,
	min_expression_fc = 0,
	add_expression_constant = 1,
	transformation_fun = log2,
	gene_chr_limit = 70,
	use_raster = TRUE,
	verbose=FALSE,
	print_border=TRUE,
	show_clust_borders=TRUE
) {
	
	require("ComplexHeatmap")
	require("stringr")
	
	# PLOT METACELL PROFILE
	message("Plotting metacell expression")
	
	# get variables necessary for hm definition (from previous function call)
	genes = marker_data_list$genes
	gene_ord = marker_data_list$gene_ord
	clust_ord = marker_data_list$clust_ord
	niche_geomean_n = marker_data_list$niche_geomean_n
	gene_labels_1 = marker_data_list$gene_labels_1
	gene_labels_2 = marker_data_list$gene_labels_2
	
	# truncate left-side annotations
	gene_labels_1 = stringr::str_trunc(gene_labels_1, gene_chr_limit)
	gene_labels_1 = stringr::str_pad(gene_labels_1, gene_chr_limit, side="right")
	# truncate right annotations
	gene_labels_2 = stringr::str_trunc(gene_labels_2, gene_chr_limit)
	gene_labels_2 = stringr::str_pad(gene_labels_2, gene_chr_limit, side="left")
	
	# define matrix per metacell, based on geometric means
	if (!is.null(transformation_fun)) {
		mat1 = transformation_fun( niche_geomean_n[genes[gene_ord],as.character(clust_ord)] + add_expression_constant )
		mat1 = pmin( mat1, max_expression_fc )
		mat1 = pmax( mat1, min_expression_fc )
	} else {
		mat1 = niche_geomean_n[genes[gene_ord],as.character(clust_ord)] + add_expression_constant
		mat1 = pmin( mat1, max_expression_fc )
		mat1 = pmax( mat1, min_expression_fc )
	}
	
	# create gene annotations
	# if we have a vector of
	if (!is.null(highlight_genes)) {
		highlight_list = scp_highlight_genes_function(gene_labels = names(marker_data_list$gene_labels_1), highlight_genes = highlight_genes, highlight_color = highlight_color)
		gids = highlight_list$highlight_ixs
		gene_font_col = highlight_list$gene_font_col
	} else {
		gene_font_col = rep("black", length(marker_data_list$gene_labels_1))
	}
	
	if (show_gene_names) {
		if (!is.null(highlight_genes) & omit_unhighlighted) {
			# if we have a vector of genes to highlight and we don't want to plot the other genes...
			message("Gene annots highlights")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"), 
				gene = anno_mark(which="row", side="right", at=gids, labels=gene_labels_1[gids], labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"),
				gene = anno_mark(which="row", side="left", at=gids, labels=gene_labels_2[gids],  labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
		} else {
			# default scenario: plot all genes (and if there are genes to highlight, paint them in a different color)
			message("Gene annots")
			if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_1, location = 0, just = "left", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
			if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_2, location = 1, just = "right", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
		}
	} else {
		row_ha_left = ComplexHeatmap::HeatmapAnnotation(
			which = "row", empty = anno_empty(which = "row", border = FALSE)
		)
		row_ha_right = row_ha_left
	}
	
	# gene colour labels
	if (!is.null(clust_row)) {
		message("Row colour annots...")
		if (!is.null(clust_row_color)) {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1),,drop=FALSE], col = clust_row_color,
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		} else {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1),,drop=FALSE], 
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		}
		row_ha_left = c(row_ha_left, row_col_ha, gap = clust_anno_gap)
		row_ha_right = c(row_col_ha, row_ha_right, gap = clust_anno_gap)
	}
	
	
	message("Expression colors...")
	col_fun = circlize::colorRamp2(
		breaks = seq(from = min_expression_fc, to = max_expression_fc, length.out = length(heatmap_colors)), 
		colors = heatmap_colors)
	# shades = col_fun(seq(from = min_expression_fc, to = max_expression_fc, length.out = 40))
	
	# mc labels
	message("Metacell labels...")
	collabs <- colnames(mat1)
	if (show_mc_names) { 
		#collabs <- rep("",length(collabs))
		column_lab_ha = ComplexHeatmap::HeatmapAnnotation(
			which = "column",
			LAB = anno_text(which = "column", collabs, gp = gpar(fontsize = mc_font_size, rot=90)),
			height=unit(2,"mm")
		)
	} else {
		column_lab_ha = ComplexHeatmap::HeatmapAnnotation(
			which = "column", empty = anno_empty(which = "column", border = FALSE),
			height = clust_anno_gap 
		)
	}
	top_column_ha = c(column_lab_ha)
	bottom_column_ha = c(column_lab_ha)
	
	# column color annotation
	if (!is.null(clust_col)) {
		message("Column colour annots...")
		
		if (class(clust_col) == "character") {
			
			if (is.null(names(clust_col))) {
				names(clust_col) <- clust_ord
			} else {
				if (!all(names(clust_col) %in% clust_ord))
					stop("Colour and cluster names do not match!")
				clust_col <- clust_col[clust_ord]
			}
		  annot_labs = as.character(unique(names(clust_col)))
			column_col_ha = ComplexHeatmap::HeatmapAnnotation(
				which = "column",
				"cluster" = colnames(mat1),
				col = list("cluster" = clust_col),
				border = TRUE,
				simple_anno_size = clust_anno_size,
				height = unit(1,"mm"),
				show_annotation_name = TRUE, show_legend = FALSE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size),
				annotation_legend_param = list(cluster = list(
				  at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120))
				))
			)
			
		} else if ("data.frame" %in% class(clust_col)) {
		  annotation_legend_order <- sapply(colnames(clust_col_df), function(cn) {
		    annot_labs = as.character(unique(clust_col_df[,cn]))
		    list(at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120)))
		  }, USE.NAMES = TRUE, simplify = FALSE)
			if (!is.null(clust_col_color)) {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col, col = clust_col_color,
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
					annotation_name_gp = gpar(fontsize = mc_font_size),
					annotation_legend_param = annotation_legend_order
				)
			} else {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col, 
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
					annotation_name_gp = gpar(fontsize = mc_font_size),
					annotation_legend_param = annotation_legend_order
				)
			}
			
		}
		
		top_column_ha <- c(column_col_ha,top_column_ha, gap = clust_anno_gap)
		bottom_column_ha = c(bottom_column_ha,column_col_ha, gap = clust_anno_gap)
	}
	
	# barplots
	if (!is.null(clust_bars)) {
		
		if (is.null(names(clust_bars)))
			names(clust_bars) <- as.character(clust_ord)
		anno_bar <- clust_bars[as.character(clust_ord)]
		baxl <- range(clust_bars)
		
		column_bar_ha <- ComplexHeatmap::HeatmapAnnotation(
			which = "column",
			BAR = anno_barplot(
				anno_bar, height = 3 * clust_anno_size, bar_width = 0.9,
				gp = gpar(fill = clust_bars_color, col = clust_bars_color, fontsize = mc_font_size),
				axis_param = list(gp = gpar(fontsize = mc_font_size), at = baxl, labels = baxl)
			),
			show_annotation_name = FALSE, show_legend = FALSE, gap = clust_anno_gap
		)
		top_column_ha = c(column_bar_ha,top_column_ha, gap = clust_anno_gap)
		
	}
	
	# expression heatmap
	h1 = ComplexHeatmap::Heatmap(
		mat1, name = mat_label, col = col_fun, use_raster = use_raster,
		cluster_rows = FALSE, cluster_columns = FALSE,
		width = width,
		height = height,
		column_title = sprintf( "%i metacells", ncol(mat1) ),
		row_title = sprintf( "%i marker genes", nrow(mat1) ),
		show_column_names = FALSE,
		show_row_names = FALSE,
		right_annotation = row_ha_right,
		left_annotation = row_ha_left,
		top_annotation = top_column_ha,
		bottom_annotation = bottom_column_ha,
		column_names_gp = gpar(fontsize = mc_font_size),
		show_heatmap_legend = show_heatmap_legend,
		border = print_border
	)
	
	# save figure
	plotting_function(output_file, width, height, res, EXP = {
		
		# draw heatmap
		draw(h1)
		
		# add cluster borders
		if (!is.null(clust_col) & show_clust_borders & class(clust_col) == "character") {
			mat2 <- rbind(clust_col[1][match(clust_ord, rownames(clust_col))])
		} else if (!is.null(clust_col) & show_clust_borders & class(clust_col) %in% "data.frame") {
			mat2 <- rbind(clust_col[,1][match(clust_ord, rownames(clust_col))])
		} else {
			mat2 <- rbind(clust_ord)
		}
		if (show_clust_borders) {
			change_clust <- which(sapply(2:ncol(mat2), function(i) mat2[,i] != mat2[,i - 1]))
			decorate_heatmap_body(mat_label, {
				for (i in change_clust) {
					grid.lines(x = i / ncol(mat2), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
				}
			})
		}
	})
	
	return(h1)
	
	message("Metacell heatmap done")
	
}




#' Plot heatmap of gene expression for single cells
#'
#' Plots heatmap of gene expression fold change for metacells. Requires selecting markers and ordering first with `scp_plot_cmod_markers_select`
#'
#' @param marker_data_list character, list object produced by the heatmap pre-processing function `scp_plot_cmod_markers_select`
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param show_gene_names logical, show gene names as heatmap row names (default: FALSE)
#' @param highlight_genes vector of gene names to highlight in the heatmap (in a different color).
#' @param highlight_color highlight color for genes in `highlight_genes`; default is "blue"
#' @param omit_unhighlighted if TRUE, do not plot any other gene name other than the ones in `highlight_genes` (default: FALSE)
#' @param gene_font_size numeric, size of the gene names plotted as rownames
#'    (default: 4)
#' @param gene_chr_limit numeric, limit gene annotations to given number of charaters
#' @param clust_row data.frame with gene colour annotations in columns,
#'    rownames should be `marker_data_list$gids`; 
#'    if `NULL` (default), row colour annotation bar is not printed
#' @param clust_row_color named list of color mappings for annotations in `clust_row`; 
#'    names should be `colnames(clust_row)`
#' @param clust_col either a character vector or data.frame with cell colour annotations;
#'    if `clust_col` is a vector, it should be of the same length and in the same order as `gene_list$clust_ord`, 
#'    and if it is named, the names should be metacell names;
#'    if `clust_col` is a data.frame, it should have single cells as rownames, annotations (not colors) should be in columns, 
#'    and color mappings should be given by `clust_col_color` argument;
#'    if `NULL` (default), cluster colour annotation bar is not printed
#' @param clust_col_color named list of color mappings for annotations in `clust_col` when it is a data.frame; 
#'    names should be `colnames(clust_col)`
#' @param clust_bars numeric, optional values to be ploted as bars annotation on top of
#'    heatmap columns, length should be same as `length(mc@mc)` (default:  NULL)
#' @param clust_bars_color` character, color for barplots, can be either single color
#'    or a (named) vector (names should be `gene_list$clust_ord`)
#' @param clust_anno_size unit, height of the column annotation bar (default: `unit(1,"mm")`)
#' @param show_mc_names logical, show metacell names as heatmap column names
#'    (default: TRUE) ~NOT IMPLEMENTED~
#' @param use_raster logical, whether to rasterise heatmap bodies (default: TRUE) 
#' @param mc_font_size numeric, size of the metacell names (default: 4)
#' @param heatmap_colors vector of colors to use for heatmap coloring function (low to high; default: `c("white","white","orange","red","purple","black")`)
#' @param min_expression_fc,max_expression_fc numeric, max and min expression values to scale single cell heatmap coloring to (default: 0 and 5)
#'
scp_plot_cmod_markers_sc <- function(
	marker_data_list,
	mc_object,
	mat_object,
	output_file = NULL,
	height = 10,
	width = 5,
	res = NA,
	show_heatmap_legend = TRUE,
	show_gene_names = FALSE,
	highlight_genes = NULL,
	highlight_color = "blue",
	omit_unhighlighted = FALSE,
	gene_font_size = 5,
	clust_row = NULL, clust_row_color = NULL,
	clust_col = NULL, clust_col_color = NULL,
	clust_bars = NULL, clust_bars_color="grey",
	clust_anno_size = unit(4,"mm"), clust_anno_gap = unit(1,"mm"),
	show_mc_names = TRUE, mc_font_size = 5,
	heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
	gene_chr_limit = 70,
	verbose=FALSE,
	smoothen = 5,
	min_expression_fc = 0,
	max_expression_fc = 5,
	use_raster = TRUE,
	print_border=TRUE,
	show_clust_borders = TRUE
	
) {
	
	# get variables necessary for hm definition (from previous function call)
	genes = marker_data_list$genes
	gene_ord = marker_data_list$gene_ord
	clust_ord = marker_data_list$clust_ord
	niche_geomean_n = marker_data_list$niche_geomean_n
	gene_labels_1 = marker_data_list$gene_labels_1
	gene_labels_2 = marker_data_list$gene_labels_2
	
	# truncate left-side annotations
	gene_labels_1 = stringr::str_trunc(gene_labels_1, gene_chr_limit)
	gene_labels_1 = stringr::str_pad(gene_labels_1, gene_chr_limit, side="right")
	# truncate right annotations
	gene_labels_2 = stringr::str_trunc(gene_labels_2, gene_chr_limit)
	gene_labels_2 = stringr::str_pad(gene_labels_2, gene_chr_limit, side="left")
	
	# directory to save output files to
	if (is.null(output_file)) {
		outdir <- getwd()
	} else {
		outdir <- dirname(output_file)
	}
	
	
	# PLOT SINGLE-CELL PROFILE
	message("Plotting single cell expression")
	cell_order=c()
	for (niche in clust_ord){
		cells=names(mc_object@mc[which(mc_object@mc == niche)])
		cell_order=c(cell_order,cells)
	}
	cluster_cell_count=as.matrix(table(mc_object@mc))
	n_cells_cluster=cluster_cell_count[clust_ord,1]
	cells_clusts <- unlist(mapply(rep, clust_ord, n_cells_cluster, USE.NAMES=FALSE))
	
	umis=as.matrix(mat_object@mat[,names(mc_object@mc)])
	mat = umis[genes, cell_order]
	totu = colSums(umis[, cell_order])
	mat = t(t(mat) / totu) * 800
	
	lus_1 = log2(1 + 7 * mat[genes[gene_ord], cell_order])
	lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
	lus_smoo = t(apply(lus[genes[gene_ord],cell_order], 1, function(x) rollmean(x, smoothen, fill=0)))
	
	# heatmap per metacell
	mat1sc <- pmin(lus_smoo, max_expression_fc)
	mat1sc <- pmax(lus_smoo, min_expression_fc)
	colnames(mat1sc) <- colnames(lus_smoo)
	
	# define matrix per metacell
	mat1 <- niche_geomean_n[genes[gene_ord],as.character(clust_ord)]
	
	# colors in heatmap
	col_fun = circlize::colorRamp2(
		breaks = seq(from = min_expression_fc, to = max_expression_fc, length.out = length(heatmap_colors)), 
		colors = heatmap_colors)
	# shades = col_fun(seq(from = min_expression_fc, to = max_expression_fc, length.out = 40))
	
	# mc labels
	top_column_ha <- HeatmapAnnotation(
		mclabstop = anno_empty(which = "column", border = FALSE, height = unit(2,"mm"))
	)
	bottom_column_ha <- HeatmapAnnotation(
		mclabsbottom = anno_empty(which = "column", border = FALSE, height = unit(2,"mm"))
	)
	
	# Column colour annotation
	if (!is.null(clust_col)) {
		message("Column colour annots...")
		
		if (class(clust_col) == "character") {
			
			if (is.null(names(clust_col))) {
				names(clust_col) = clust_ord
			}
		  annot_labs = as.character(unique(names(clust_col)))
			column_col_ha <- HeatmapAnnotation(
				which = "column",
				"cluster" = cells_clusts, col = list("cluster" = clust_col),
				border = c(TRUE),
				simple_anno_size = clust_anno_size,
				show_annotation_name = TRUE, show_legend = FALSE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size),
				annotation_legend_param = list(cluster = list(
				  at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120))
				))
			)
			top_column_ha <- c(column_col_ha,top_column_ha, gap = clust_anno_gap)
			bottom_column_ha <- c(bottom_column_ha,column_col_ha, gap = clust_anno_gap)
			
		} else if ("data.frame" %in% class(clust_col)) {
		  annotation_legend_order <- sapply(colnames(clust_col_df), function(cn) {
		    annot_labs = as.character(unique(clust_col_df[,cn]))
		    list(at = annot_labs, labels = annot_labs, ncol = pmax(1,round(length(annot_labs)/120)))
		  }, USE.NAMES = TRUE, simplify = FALSE)
			if (!is.null(clust_col_color)) {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col[colnames(mat1sc),,drop=FALSE], col = clust_col_color,
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
					annotation_name_gp = gpar(fontsize = mc_font_size),
					annotation_legend_param = annotation_legend_order
				)
			} else {
				column_col_ha <- HeatmapAnnotation(
					which = "column",
					df = clust_col[colnames(mat1sc),,drop=FALSE], 
					border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
					show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
					annotation_name_gp = gpar(fontsize = mc_font_size),
					annotation_legend_param = annotation_legend_order
				)
			}
			top_column_ha <- c(column_col_ha,top_column_ha, gap = clust_anno_gap)
			bottom_column_ha = c(bottom_column_ha,column_col_ha, gap = clust_anno_gap)
		}
		
	}	
	# barplot annotation
	if (!is.null(clust_bars)) {
		
		cell_cols <- unlist(mapply(rep, clust_bars_color, n_cells_cluster, USE.NAMES=FALSE))
		cell_col <- structure(cell_cols, names=cell_order)
		
		if (is.null(names(clust_bars)))
			names(clust_bars) <- as.character(cell_order)
		
		anno_bar <- clust_bars[as.character(cell_order)]
		baxl <- range(clust_bars)
		
		column_bar_ha <- ComplexHeatmap::HeatmapAnnotation(
			which = "column",
			BAR = anno_barplot(
				anno_bar, height = 3 * clust_anno_size, bar_width = 1, 
				gp = gpar(fill = cell_col, col = cell_col, fontsize = mc_font_size),
				axis_param = list(gp = gpar(fontsize = mc_font_size), at = baxl, labels = baxl)
			),
			show_annotation_name = FALSE, show_legend = FALSE, gap = clust_anno_gap
		)
		top_column_ha = c(column_bar_ha,top_column_ha)
		
	}
	
	# create gene annotations
	# if we have a vector of genes to highlight...
	if (!is.null(highlight_genes)) {
		highlight_list = scp_highlight_genes_function(gene_labels = names(marker_data_list$gene_labels_1), highlight_genes = highlight_genes, highlight_color = highlight_color)
		gids = highlight_list$highlight_ixs
		gene_font_col = highlight_list$gene_font_col
	} else {
		gene_font_col = rep("black", length(marker_data_list$gene_labels_1))
	}
	
	if (show_gene_names) {
		if (!is.null(highlight_genes) & omit_unhighlighted) {
			# if we have a vector of genes to highlight and we don't want to plot the other genes...
			message("Gene annots highlights")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"), 
				gene = anno_mark(which="row", side="right", at=gids, labels=gene_labels_1[gids], labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row", simple_anno_size = unit(1,"mm"),
				gene = anno_mark(which="row", side="left", at=gids, labels=gene_labels_2[gids],  labels_gp=gpar(fontsize = gene_font_size, col = gene_font_col[gids]), extend=unit(0.5, "mm")))
		} else {
			# default scenario: plot all genes (and if there are genes to highlight, paint them in a different color)
			message("Gene annots")
			if (verbose) message(paste(head(gene_labels_1),collapse=", "), ",...")
			row_ha_right = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_1, location = 0, just = "left", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
			if (verbose) message(paste(head(gene_labels_2),collapse=", "), ",...")
			row_ha_left = ComplexHeatmap::HeatmapAnnotation(
				which = "row",
				gene = anno_text(which = "row", gene_labels_2, location = 1, just = "right", gp = gpar(fontsize = gene_font_size, col = gene_font_col)))
		}
	} else {
		row_ha_left = ComplexHeatmap::HeatmapAnnotation(
			which = "row", empty = anno_empty(which = "row", border = FALSE)
		)
		row_ha_right = row_ha_left
	}
	
	# gene colour labels
	if (!is.null(clust_row)) {
		message("Row colour annots...")
		if (!is.null(clust_row_color)) {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1sc),,drop=FALSE], col = clust_row_color,
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap, 
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		} else {
			row_col_ha <- HeatmapAnnotation(
				which = "row",
				df = clust_row[rownames(mat1sc),,drop=FALSE], 
				border = TRUE, simple_anno_size = clust_anno_size, title=NULL,
				show_annotation_name = TRUE, show_legend = TRUE, gap = clust_anno_gap,
				annotation_name_gp = gpar(fontsize = mc_font_size)
			)
		}
		row_ha_left = c(row_ha_left, row_col_ha, gap = clust_anno_gap)
		row_ha_right = c(row_col_ha, row_ha_right, gap = clust_anno_gap)
	}
	
	# expression heatmap
	h1sc <- Heatmap(
		mat1sc, name = "sc_expression", col = col_fun, use_raster = use_raster,
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		show_column_names = FALSE,
		show_row_names = FALSE,
		width = width,
		height = height,
		column_title = sprintf( "%i cells", ncol(mat1sc) ),
		row_title = sprintf( "%i marker genes", nrow(mat1sc) ),
		right_annotation = row_ha_right,
		left_annotation = row_ha_left,
		top_annotation = top_column_ha,
		bottom_annotation = bottom_column_ha,
		show_heatmap_legend = show_heatmap_legend,
		border = print_border
	)
	hlistsc <- h1sc
	
	# save figure
	
	plotting_function(output_file, width, height, res, EXP={
		# heatmap
		# draw(hlistsc, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
		# drop all this artificial padding...
		draw(hlistsc)
		
		# add labels of metacells
		if (!is.null(clust_col)) {
			if (class(clust_col) == "character") {
				mat2sc <- rbind(clust_col[match(cells_clusts, names(clust_col))])
			} else if ("data.frame" %in% class(clust_col)) {
				mat2sc <- rbind(clust_col[colnames(mat1sc),colnames(clust_col)[1]])
			}
		} else {
			mat2sc <- rbind(cells_clusts)
		}
		if (is.null(colnames(mat2sc)))
			colnames(mat2sc) <- cells_clusts
		change_clust_sc <- which(sapply(2:ncol(mat2sc), function(i) mat2sc[,i] != mat2sc[,i - 1]))
		change_mc_sc <- c(
			which(sapply(2:ncol(mat2sc), function(i) colnames(mat2sc)[i] != colnames(mat2sc)[i - 1])),
			ncol(mat2sc)
		)
		if (show_clust_borders) {
			decorate_heatmap_body("sc_expression", {
				for (i in change_clust_sc) {
					grid.lines(x = i / ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
				}
				for (i in change_mc_sc) {
					grid.lines(x = i / ncol(mat2sc), y = c(0,1), gp = gpar(lty = 1, lwd = 0.5))
				}
			})
		}
		
		.add_mc_labels <- function(pos,labs) {
			for (j in 1:length(pos)) {
				i <- pos[j]
				iprev <- ifelse(j == 1,0,pos[j - 1])
				nt <- ncol(mat2sc)
				grid.text(label = labs[j], x = i / nt - (i / nt - iprev / nt) / 2, y = 0.5, gp = gpar(fontsize = mc_font_size), rot=90)
			}
		}
		if (show_mc_names) {
			decorate_annotation("mclabstop", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
			decorate_annotation("mclabsbottom", .add_mc_labels(pos=change_mc_sc, labs=colnames(mat1)))
		}
	})
	ht_opt(RESET = TRUE)
	
	message("Single-cell heatmap done")
	
}


#' Gene UMI count in metacell
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#' @param T_totumi integer, filter out genes with less than `T_totumi` counts
#' @param grouping_vector a vector with grouping categories for individual cells
#'    (same format and length as `mc@mc`). NULL by default, which means that `mc@mc`
#'    is used as a grouping factor.
#'
#' @return matrix of gene UMI counts in metacells
#'
sca_mc_gene_counts = function(mc_object, mat_object, T_totumi=0, grouping_vector = NULL) {
	
	# get UMIs
	scr_umis = mat_object@mat
	keep_genes_bool = Matrix::rowSums(scr_umis) > T_totumi
	keep_genes = names(which(keep_genes_bool))
	
	# counts per metacell or grouping vector
	if (is.null(grouping_vector)) {
		grouping_vector = mc_object@mc
		grouping_vector_order = colnames(mc_object@mc_fp)
	} else if ("factor" %in% class(grouping_vector)) {
		grouping_vector_order = levels(grouping_vector)
	} else {
		grouping_vector_order = unique(grouping_vector)
	}
	niche_counts = .row_stats_by_factor( as.matrix(scr_umis[ keep_genes, names(mc_object@mc) ]), grouping_vector, rowSums )
	niche_counts = niche_counts[ , grouping_vector_order ]
	return(niche_counts)
	
}


#' Gene UMI fraction in metacell
#'
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mc_counts matrix of gene UMI counts in metacells (output of `sca_mc_gene_counts()`)
#'
#' @return matrix of gene UMI fractions in metacells
#'
sca_mc_gene_umifrac = function(mc_object, mc_counts, multiplying_factor = 1000){
	niche_counts = mc_counts
	niche_totals = Matrix::colSums(niche_counts)
	niche_umifrac = t(apply(niche_counts,1,function(x) x * multiplying_factor / niche_totals))
	niche_umifrac = niche_umifrac[,colnames(mc_object@mc_fp)]
	return(niche_umifrac)
}


#' Gene UMI fraction in metacell (object-independent)
#'
#' @param mc_counts matrix of gene UMI counts in metacells (output of `sca_mc_gene_counts()`)
#' @param mc_columns vector indicating order of metacells (defaults to column names in `mc_counts`, another common value could be `colnames(mc@mc_fp)`)
#'
#' @return matrix of gene UMI fractions in metacells
#'
sca_mc_gene_umifrac_noobj = function(mc_counts, mc_columns = colnames(mc_counts), multiplying_factor = 1000){
	niche_counts = mc_counts
	niche_totals = colSums(niche_counts)
	niche_umifrac = t(apply(niche_counts,1,function(x) x * multiplying_factor / niche_totals))
	niche_umifrac = niche_umifrac[,mc_columns]
	return(niche_umifrac)
}

#' Calculate expression conservation (EC) score for pairs of orthologous genes in a cross-species iterative comparison of coexpression (ICC)
#'
#' @param mat_sp1 expression matrix of species 1 (rows are genes, columns are any conditions)
#' @param mat_sp2 expression matrix of species 2 (rows are genes, columns are any conditions)
#' @param og_pairs either a path to a ortholog pairs table, or one such table. Requires two columns: first are genes from species 1, 
#'    second are genes from species 2
#' @param niter maximum number of ICC iterations (default = 100, min is 2)
#' @param icc_thr similarity threshold to stop ICC iteration (ICC stops when difference betwee iteration is below this value; default = 0.05)
#' @param method ICC correlation method (default is `pearson`)
#' @param do_duplicates bool (default `FALSE`): if `TRUE`, calculate EC values for best reciprocal pairs of duplicate genes in addition 
#'    to the one-to-one ortholog pairs (which is the default)
#' @param num_o2o_ref_genes if `do_duplicates` is TRUE, select up to this number of genes for the duplicate EC scoring step. Less than 1000 is not
#'    recommended (very unstable EC values).
#' @param use_variable_o2o_genes bool; if `TRUE` (default), use only one-to-one orthologs with variable expression in both species as a reference 
#'    matrix for the duplicates' EC calculations. Variable genes will be computed from footprint matrices (default: top quartile of max footprint
#'    across cells). If input matrices are not footprint matrices, you can provide pre-computed variable genes with `vargenes_sp1` and `vargenes_sp1`
#'    instead (see below).
#' @param variable_o2o_thr if use_variable_o2o_genes is `TRUE`, use this quantile threshold to select variable genes based on the distribution of 
#'    max observed footprint across all cells (default is 0.75, i.e. the top quartile of variable genes). Increasing the size of the reference o2o
#'    matrix can result in extremely long runtimes. Decreasing it can result in poor EC score estimation for paralogs. Rough benchmarking: thr=0.75 
#'    is fairly accurate and fast (~1-2h, cor ~0.80 with whole-matrix result. A thr=0.5 is more accurate but slower (~12h, cor ~0.95).
#' @param vargenes_sp1, vargenes_sp2 (default is `NULL`): vectors of variable one-to-one orthologous genes to use as reference for the 
#'    `do_duplicates` step. Only applicable if `do_duplicates` is `TRUE`. 
#' @param nthreads_icc number of threads to use for ICC calculations (i.e. parallelisation of the correlation matrix).
#' @param nthreads_dup number of parallel processes of duplicate scoring (each of them will run with n=`nthreads_icc` threads).
#' @param cross_fc_thrs numeric, fold change threshold (default: 2, set to NULL to skip filtering by fc)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` we require in each species (default: 1, set to NULL to skip filtering by fc)
#'
#' @return dataframe of unique gene pairs from species 1 and 2, and their cross-species expression conservation score
#'
csps_markers_icc = function(
	mat_sp1,
	mat_sp2, 
	og_pairs, 
	niter = 100, 
	icc_thr = 0.05, 
	method = "pearson", 
	do_duplicates = FALSE, 
	num_o2o_ref_genes = NULL,
	use_variable_o2o_genes = TRUE, 
	variable_o2o_thr = 0.75, 
	vargenes_sp1 = NULL, 
	vargenes_sp2 = NULL,
	cross_fc_thrs = 2, 
	do_quantile_normalisation = TRUE,
	cross_n = 1,
	nthreads_icc = 2,
	nthreads_dup = 1) {
	
	require("igraph")
	require("data.table")
	
	# register cores
	registerDoParallel(cores = nthreads_dup)

	# ensure data are matrices
	mat_sp1 = as.matrix(mat_sp1)
	mat_sp2 = as.matrix(mat_sp2)
	
	# load orthologous pairs
	if ("character" %in% class(og_pairs)) {
		og_pairs = data.table::fread(og_pairs, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	}	
	colnames(og_pairs) = c("sp1", "sp2")

	# restrict pairs to expressed genes
	og_pairs = og_pairs [ og_pairs[,1] %in% rownames(mat_sp1) & og_pairs[,2] %in% rownames(mat_sp2) , ]
	
	### Find ec values of one-to-one orthologs
	# get one to one pairs (reference)
	list_o2o_sp1 = names(which(table(og_pairs[,1]) == 1))
	list_o2o_sp2 = names(which(table(og_pairs[,2]) == 1))
	bool_o2o = og_pairs[,1] %in% list_o2o_sp1 & og_pairs[,2] %in% list_o2o_sp2
	og_pairs_o2o = og_pairs [ bool_o2o, ]

	# reference o2o matrix
	mar_sp1 = mat_sp1 [ og_pairs_o2o[,1] , ]
	mar_sp2 = mat_sp2 [ og_pairs_o2o[,2] , ]
	
	# ec values for the one to one orthologs
	message(sprintf("ICC markers | pairs of one-to-one orthologs = %i", nrow(og_pairs_o2o)))
	ecv_o2o = csps_calc_icc(mat_sp1 = mar_sp1, mat_sp2 = mar_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = TRUE, num_cores = nthreads_icc)
	# ecv_o2o = ecv_o2o_list[[1]]
	ecv_o2o [ ecv_o2o$ec_value < 0 | is.na(ecv_o2o$ec_value) , "ec_value" ] = 0
	ecv_o2o$is_o2o = 1
	
	### Find best reciprocal duplicate pairs (highest ec value)
	# loop through pairs of non-o2o homologs, add them to the original o2o matrix, and calculate ec value
	if (do_duplicates) { 
		
		# select genes for the reference o2o matrix
		if (use_variable_o2o_genes & is.null(vargenes_sp1) & is.null(vargenes_sp2)) {

			# o2o orthologs with variable expression
			max_fcs_sp1 = apply(mat_sp1, 1, max, na.rm = TRUE)
			max_fcs_sp2 = apply(mat_sp2, 1, max, na.rm = TRUE)
			vargenes_sp1 = rownames(mat_sp1) [ max_fcs_sp1 > quantile(max_fcs_sp1, variable_o2o_thr) ]
			vargenes_sp2 = rownames(mat_sp2) [ max_fcs_sp2 > quantile(max_fcs_sp2, variable_o2o_thr) ]
			# variable AND o2o orthologs
			ecv_o2o_ref = ecv_o2o [ ecv_o2o[,1] %in% vargenes_sp1 & ecv_o2o[,2] %in% vargenes_sp2 , ]
			refgenes_sp1 = ecv_o2o_ref[,1]
			refgenes_sp2 = ecv_o2o_ref[,2]
		
		} else if (!is.null(vargenes_sp1) & !is.null(vargenes_sp2)) { 
			
			ecv_o2o_ref = ecv_o2o [ ecv_o2o[,1] %in% vargenes_sp1 & ecv_o2o[,2] %in% vargenes_sp2 , ]
			refgenes_sp1 = ecv_o2o_ref[,1]
			refgenes_sp2 = ecv_o2o_ref[,2]
			
		} else if (is.null(num_o2o_ref_genes)) { 
			
			refgenes_sp1 = ecv_o2o[,1]
			refgenes_sp2 = ecv_o2o[,2]
			
		} else {
			
			ixs_ref = sample(1:nrow(ecv_o2o), min(nrow(ecv_o2o), num_o2o_ref_genes), replace = FALSE)
			refgenes_sp1 = ecv_o2o[ixs_ref,1]
			refgenes_sp2 = ecv_o2o[ixs_ref,2]
			
		}
		
		# sanity check: do we have enough genes in the reference matrix?
		message(sprintf("ICC markers | reference one-to-one orthologs for duplicate EC scoring = %i", length(refgenes_sp1)))
		if (length(refgenes_sp1) == 0) { 
			stop(sprintf("Only %i pairs of one-to-one orthologs in the reference matrix for duplicate EC calculations. I can't do this.", length(refgenes_sp1)))
		} else if (length(refgenes_sp1) < 1000 & length(refgenes_sp1) > 0) {
			warning(sprintf("Only %i pairs of one-to-one orthologs in the reference matrix for duplicate EC calculations. Too few?", length(refgenes_sp1)))
		}
		
		# get duplicates (all types)
		list_dup_sp1 = names(which(table(og_pairs[,1]) > 1))
		list_dup_sp2 = names(which(table(og_pairs[,2]) > 1))
		bool_dup = og_pairs[,1] %in% list_dup_sp1 | og_pairs[,2] %in% list_dup_sp2
		og_pairs_dup = og_pairs [ bool_dup , ]
		og_pairs_dup = og_pairs_dup[ order(og_pairs_dup[,1], og_pairs_dup[,2]) , ]
		
		# graph of duplicates
		og_pairs_dup_sps = as.matrix(og_pairs_dup)
		og_pairs_dup_sps[,1] = paste("sp1", og_pairs_dup_sps[,1])
		og_pairs_dup_sps[,2] = paste("sp2", og_pairs_dup_sps[,2])
		graph_dup = igraph::graph_from_edgelist(og_pairs_dup_sps, directed = TRUE)
		graph_dup_components = igraph::components(graph_dup)
		
		# dict of duplicates
		dict_dups = list()
		list_components = unique(graph_dup_components$membership)
		for (component in list_components) {
			dups = names(graph_dup_components$membership) [graph_dup_components$membership == component]
			dict_dups[[component]] = dups
		}

		if (length(dict_dups) > 0) {

			# log progress		
			message(sprintf("ICC markers | clusters with non-o2o homologs = %i (total = %i genes)", length(dict_dups), length(unique(unlist(dict_dups)))))
			
			# loop through duplicates in dict_dups
			# this loop is a function for easy parallelisation with plyr::adply
			loop_dups = function(i) {
			  print(i)
				# get lists of genes from sp1 and sp2
				d = dict_dups[[i]]
				dg = stringr::word(d, 2)
				genes1 = dg [ stringr::word(d, 1) == "sp1" ]
				genes2 = dg [ stringr::word(d, 1) == "sp2" ]
				
				# prepare vector of ec values for this set of paralogs
				gen_dups = og_pairs [ og_pairs[,1] %in% genes1 & og_pairs[,2] %in% genes2 , ]
					
				# add test genes to reference expression matrix, sp1
				mai_sp1 = rbind( mat_sp1[ gen_dups[,1] ,], mar_sp1[refgenes_sp1,] )
				rownames(mai_sp1)[ 1:nrow(gen_dups) ] = gen_dups[,1]
				# add test genes to reference expression matrix, sp2
				mai_sp2 = rbind( mat_sp2[ gen_dups[,2] ,], mar_sp2[refgenes_sp2,] )
				rownames(mai_sp2)[ 1:nrow(gen_dups) ] = gen_dups[,2]

				# recalculate ec values using reference expression matrix and new set of duplicates
				ecv_dup_i = csps_calc_icc(mat_sp1 = mai_sp1, mat_sp2 = mai_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = FALSE, num_cores = nthreads_icc)
				ecv_dup_q = ecv_dup_i [ ecv_dup_i$sp1 %in% gen_dups[,1] & ecv_dup_i$sp2 %in% gen_dups[,2] , ]
						
				return(ecv_dup_q)
				
			}
			
			# apply loop (not parallelised internally)
			ecv_dup = plyr::adply(.data = 1:length(dict_dups), .margins = 1,  .fun = loop_dups, .parallel = TRUE, .id = "cluster", .progress = "none")
			ecv_dup = ecv_dup [ , c("sp1", "sp2", "ec_value", "cluster") ]
			
			# keep best reciprocal pairs of orthologous genes
			ecv_dup_t = data.table::as.data.table(ecv_dup)
			ecv_dup_t$ec_value [ is.na(ecv_dup_t$ec_value) ] = -1
			ecv_dup_f1 = ecv_dup_t[ ecv_dup_t[, .I[ ec_value == max(ec_value) ], by=c("sp1","cluster")]$V1 ]
			ecv_dup_f2 = ecv_dup_t[ ecv_dup_t[, .I[ ec_value == max(ec_value) ], by=c("sp2","cluster")]$V1 ]
			ecv_dup_f  = rbind(ecv_dup_f1,ecv_dup_f2)
			ecv_dup_f  = ecv_dup_f [ duplicated(ecv_dup_f) & ecv_dup_f$ec_value >= 0 , ]
			
			# break ties at random (so as to obtain a truly one-to-one table)
			ecv_dup_f = ecv_dup_f [ !duplicated(ecv_dup_f$sp1), ]
			ecv_dup_f = ecv_dup_f [ !duplicated(ecv_dup_f$sp2), ]
			
			# log progress		
			message(sprintf("ICC markers | pairs of non-o2o homologs kept = %i", nrow(ecv_dup_f)))
			
			# output contains both o2o orthologs & best reciprocal pairs of duplicated genes
			ecv_dup_f = ecv_dup_f [ , c("sp1","sp2","ec_value") ]
			ecv_dup_f$is_o2o = 0
			ecv_out = rbind(ecv_o2o, ecv_dup_f)
		
		} else{
		
			# log progress		
			message(sprintf("ICC markers | clusters with non-o2o homologs = %i (total = %i genes) | Can't score duplicates!", length(dict_dups), length(unique(unlist(dict_dups)))))
			# output will only contain one-to-one orthologs
			ecv_out = ecv_o2o
			ecv_dup = NULL
			
		}

		
	} else {
		
		# output will only contain one-to-one orthologs
		ecv_out = ecv_o2o
		ecv_dup = NULL
	
	}
	
	
	#### Output ####

	# create cross-species object
	csps = list(
		merged = cbind(mat_sp1[ ecv_out[,1], ], mat_sp2[ ecv_out[,2], ]),
		og_pairs = ecv_out[,c(1,2)],
		og_pairs_is_o2o   = ecv_out$is_o2o,
		og_pairs_ec_value = ecv_out$ec_value,
		sp1 = mat_sp1[ ecv_out[,1], ], 
		sp2 = mat_sp2[ ecv_out[,2], ],
		top_cross_sp1 = NULL, 
		top_cross_sp2 = NULL
	)
	
	# do quantile normalisation on merged matrix?
	if (do_quantile_normalisation) {
		csps$merged = quantile_normalisation(csps$merged)
	}
	
	# find covariable genes
	if (!is.null(cross_fc_thrs) & !is.null(cross_n)) {
		
		# get covariable genes
		message(sprintf("ICC markers | select genes that are variable in both species, fc %.2f in %i samples", cross_fc_thrs, cross_n))
		cross_variable_genes = csps_select_covariable_genes(
			sp1_fp = csps$sp1,
			sp2_fp = csps$sp2, 
			merged = csps$merged,
			cross_fc_thrs = cross_fc_thrs,
			cross_n = cross_n)
		
		top_cross_sp1 = cross_variable_genes$sp1
		top_cross_sp2 = cross_variable_genes$sp2
			
	} else {

		# omit getting covariable genes, get them all instead
		top_cross_sp1 = rownames(ecv_out[,1])
		top_cross_sp2 = rownames(ecv_out[,2])
		
	}
	csps$top_cross_sp1 = top_cross_sp1
	csps$top_cross_sp2 = top_cross_sp2
	
	# return data.frame with EC scores for pairs of genes
	message(sprintf("ICC markers | total num of markers = %i", nrow(ecv_out)))
	return(list(csps = csps, ec_markers = ecv_out, ec_duplicates = ecv_dup))
	
	# TODO: regress-out effect of total expression from final ec values?

}


# Select genes that are variable in two species, defined as having fc above a certain threshold in at least n samples/mcs/cell types/etc.
#' 
#' @param sp1_fp,sp2_fp,merged matrices of gene (rownames) expression over mcs/cell types/etc (colnames). Rows must be paired, i.e. correspond to comparable genes in the same order
#' @param method character, one of "min_fc" (default: keep genes with fp>`cross_fc_thrs` in at least `cross_n` cell types/metacells), "top_fc_each" (keep up to n=`genes_n` markers with fp>`cross_fc_thrs` in *each* cell type), or "top_fc_pair" (keep up to n=`genes_n` markers with fp>`cross_fc_thrs` in each pair of cross-sps cell types).
#' @param cross_fc_thrs numeric, fold change threshold (default: 2)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` we require in each species (default: 1); only relevant if `method="min_fc"`
#' @param genes_n integer, keep up to this many genes for each cell type (default: 100); only relevant if `method="top_fc_each"`
#' 
#' @return a list with two vectors of variable genes from species 1 and 2, respectively
#' 
csps_select_covariable_genes = function(
	sp1_fp, 
	sp2_fp,
	merged,
	method = "min_fc",
	cross_fc_thrs = 2,
	cross_n = 1,
	genes_n = 200) {
	
	# get submatrices from merged matrix (it's not the same as the unmerged matrices because this has been quantile normalised)
	merged_sp1 = merged[,1:ncol(sp1_fp)]
	merged_sp2 = merged[,(ncol(sp1_fp) + 1):ncol(merged) ]
	
	if (method == "min_fc") {
		
		top_cross_ix = which(
			apply(merged_sp1, 1, function(x) sort(x, decreasing=TRUE)[cross_n]) > cross_fc_thrs & 
			apply(merged_sp2, 1, function(x) sort(x, decreasing=TRUE)[cross_n]) > cross_fc_thrs
		)
		top_cross_sp1 = rownames(sp1_fp)[top_cross_ix]
		top_cross_sp2 = rownames(sp2_fp)[top_cross_ix]
	
	} else if (method == "top_fc_each") {
		
		# for each column in the merged matrix, get up to `genes_n` marker genes with a fc>thrs
		top_cross_sp1 = unique(unlist(lapply(
			1:ncol(merged), 
			function(i) {
				fpso = sort(merged[,i], decreasing = TRUE)
				gnso = names(fpso) [ fpso > cross_fc_thrs ] [ 1:genes_n ]
				return(gnso)
			}
		)))
		# reorder sp1
		top_cross_sp1 = top_cross_sp1 [ !is.na(top_cross_sp1) ]
		top_cross_sp1 = top_cross_sp1 [ order( match(top_cross_sp1, rownames(sp1_fp)) ) ]
		# get sp2
		top_cross_sp2 = rownames(sp2_fp) [ rownames(sp1_fp) %in% top_cross_sp1 ]
		
	} else if (method == "top_fc_pair") {
		
		# for each pair of columns in the merged matrix, get up to `genes_n` marker genes with a fc>thrs
		top_cross_sp1 = c()
		for (ni in 1:ncol(merged_sp1)) {
			fpso_i = sort(merged_sp1[,ni], decreasing = TRUE)
			gnso_i = names(fpso_i) [ fpso_i > cross_fc_thrs ] [ 1:genes_n ]
			for (nj in 1:ncol(merged_sp2)) {
				fpso_j = sort(merged_sp2[,nj], decreasing = TRUE)
				gnso_j = names(fpso_j) [ fpso_i > cross_fc_thrs ] [ 1:genes_n ]
				gnso_c = intersect(gnso_i, gnso_j)
				top_cross_sp1 = c(top_cross_sp1, gnso_c)
			}
		}
		# reorder sp1
		top_cross_sp1 = top_cross_sp1 [ !is.na(top_cross_sp1) ]
		top_cross_sp1 = unique(top_cross_sp1)
		top_cross_sp1 = top_cross_sp1 [ order( match(top_cross_sp1, rownames(sp1_fp)) ) ]
		# get sp2
		top_cross_sp2 = rownames(sp2_fp) [ rownames(sp1_fp) %in% top_cross_sp1 ]
		
	} else {
		
		stop("`method` has to be either `min_fc` (default) or `top_fc_each`")
		
	}
	
	# return
	return(list(sp1 = top_cross_sp1, sp2 = top_cross_sp2))

}

#' Calculate expression conservation scores with ICC, from symmetrical tables (same genes & order in both; rownames can be different)
#'
#' @param mat_sp1 expression matrix of species 1 (rows are genes, columns are any conditions). Row order in matrix 1 and 2 must be matched (pairs of orthologs), but row names need not be.
#' @param mat_sp2 expression matrix of species 2 (rows are genes, columns are any conditions). Row order in matrix 1 and 2 must be matched (pairs of orthologs), but row names need not be.
#' @param niter maximum number of ICC iterations (default = 100, min is 2)
#' @param icc_thr similarity threshold to stop ICC iteration (ICC stops when difference betwee iteration is below this value; default = 0.05)
#' @param method ICC correlation method (default is `pearson`)
#' @param verbose is an adjective
#'
#' @return dataframe of unique gene pairs from species 1 and 2, and their cross-species expression conservation score
#'
csps_calc_icc = function(mat_sp1, mat_sp2, niter = 100, icc_thr = 0.05, method = "pearson", verbose = TRUE, num_cores = 2) {
	
	require("WGCNA")
	require("tgstat")
	
	# coexpression correlation matrices
	exc_sp1 = WGCNA::cor(t(mat_sp1), method = method, nThreads = num_cores)
	exc_sp2 = WGCNA::cor(t(mat_sp2), method = method, nThreads = num_cores)
	
	# iteration 0: populate with the original cross-species matrix
	# cross-species matrix (single-species matrices must be symmetrical)
	i = 1	
	if (verbose) {
		message(sprintf("ICC iteration %i" , i - 1))
	}
	exc_csp = WGCNA::cor(exc_sp1, exc_sp2, method = method, nThreads = num_cores)
	
	# store expression correlation matrices
	# ecv object
	icc_ecv = list()
	icc_ecv[[i]] = diag(exc_csp)
	# per-species correlations
	exp_sp1 = exc_sp1
	exp_sp2 = exc_sp2
	
	for (i in 2:niter) {
		
		# which per-gene ec values are above zero?
		ixs = which(icc_ecv[[i - 1]] > 0)
		
		# for values above zero, keep matrix values from previous iteration
		exi_sp1 = matrix(0, nrow = nrow(exp_sp1), ncol = ncol(exp_sp1))
		exi_sp1[ixs,ixs] = exp_sp1[ixs,ixs]
		exi_sp2 = matrix(0, nrow = nrow(exp_sp2), ncol = ncol(exp_sp2))
		exi_sp2[ixs,ixs] = exp_sp2[ixs,ixs]
		
		# create vector of weights: previous if positive, zero otherwise
		exi_weights = rep(0, length.out = length(icc_ecv[[i - 1]]))
		exi_weights[ixs] = icc_ecv[[i - 1]][ixs]
		
		# weighted correlation
		icc_eci = WGCNA::cor(x = exi_sp1, y = exi_sp2, weights.x = exi_weights, weights.y = exi_weights, method = method, nThreads = num_cores)
		icc_ecv[[i]] = diag(icc_eci)
		
		# calculate delta icc value between iterations
		icc_delta = sum( (icc_ecv[[i]][ixs] - icc_ecv[[i - 1]][ixs] ) ^ 2 )
		if (verbose) {
			message(sprintf("ICC iteration %i | delta = %.1e" , i - 1, icc_delta))
		}
		
		# store expression correlation matrices for next iteration
		exp_sp1 = exi_sp1
		exp_sp2 = exi_sp2

		# quit loop as soon as delta icc value drops below threshold
		if (icc_delta < icc_thr) {
			break
		}
		
	}
	
	# return dataframe with gene pairs and their ec score
	ec_values = data.frame( sp1 = rownames(exc_sp1), sp2 = rownames(exc_sp2), ec_value = icc_ecv[[i]] )
	return(ec_values)
	# return(list(ec_values, exc_csp))
	
}

#' Compute cell type footprint based on a standard metacell->cell type
#' definition table, using same strategy as mc_fp.
#'
#' @param input_table data.frame with three columns: metacell, cell_type and color,
#'    or a path to tsv file with this table
#' @param mc_object loaded metacell object (`gMCCov` class)
#' @param mat_object loaded single cell matrix object (`tgScMat` class)
#'
#' @return analogous mc_object for cell types
#'
sca_cell_type_fp <- function( input_table, mc_object, mat_object, nbins=10L) {
	
	# load input table
	if (any(class(input_table) %in% "character")) {
		
		cell_type_table = read.table(input_table, header = TRUE, sep="\t", comment.char="")
		colnames(cell_type_table) = c("metacell", "cell_type", "color")
		cell_types_ordered = unique(cell_type_table$cell_type)
		rownames(cell_type_table) = cell_type_table$metacell
		
	} else if (any(class(input_table) %in% "data.frame")) {
		
		cell_type_table = input_table
		class(cell_type_table) <- "data.frame"
		colnames(cell_type_table) = c("metacell", "cell_type", "color")
		cell_types_ordered = unique(cell_type_table$cell_type)
		rownames(cell_type_table) = cell_type_table$metacell
		
	}
	
	sc_ct_label = as.vector(cell_type_table[as.character(mc_object@mc),"cell_type"])
	names(sc_ct_label) = names(mc_object@mc)
	
	cells_cols = cell_type_table[as.character(mc_object@mc),"color"]
	cells_cols = as.character(cells_cols)
	names(cells_cols) = names(mc_object@mc)
	
	# filter low expression genes not included in the mc_fp
	umis = mat_object@mat
	umis = umis[rownames(mc_object@mc_fp),names(mc_object@mc)]
	
	#ct_counts=t(apply(umis,1,function(x) tapply(x, sc_ct_label, sum)))
	#ct_size=colSums(ct_counts)
	#ct_umifrac=t(apply(ct_counts,1,function(x) x*1000/ct_size))
	#ct_umifrac_n=(0.1 + ct_umifrac)/apply(0.1 + ct_umifrac, 1, median)
	
	ct_geomean = tryCatch(
	  t(apply(umis, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1))),
	  error = function(e) {
	    warning(e)
	    message("Calculating geom mean for genes in ",nbins," bins")
	    umis_list <- vector("list",nbins)
	    sl <- split(1:nrow(umis),cut(1:nrow(umis),nbins))
	    for (i in seq_along(sl)) {
	      message(i," / ", nbins)
	      umisub <- umis[sl[[i]],]
	      ct_geomean_sub <- t(apply(umisub, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1)))
	      umis_list[[i]] <- ct_geomean_sub 
	    }
	    do.call(rbind, umis_list)
	  }
	)
	ct_geomean = ct_geomean[ , cell_types_ordered ]
	ct_meansize = tapply(Matrix::colSums(umis), sc_ct_label, mean)
	ideal_cell_size = pmin(1000,median(ct_meansize))
	g_fp = t(ideal_cell_size * t(ct_geomean) / as.vector(ct_meansize))
	fp_reg = 0.05
	g_fp_n = (fp_reg + g_fp) / apply(fp_reg + g_fp, 1, median)
	
	# for compatibility with other functions, return as a MC-like object
	ct_table=mc_object
	ct_table@mc_fp=g_fp_n
	ct_table@mc=sc_ct_label
	ct_table@colors=cells_cols
	ct_table@cell_names=names(sc_ct_label)
	return(ct_table)
	
}

#' Quantile normalization
#' @param x matrix or data.frame
quantile_normalisation <- function(x){
	df_rank <- data.frame(apply(x,2,rank,ties.method="min"))
	df_sorted <- data.frame(apply(x, 2, sort))
	df_mean <- apply(df_sorted, 1, mean)
	
	index_to_mean <- function(my_index, my_mean){
		return(my_mean[my_index])
	}
	
	df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
	rownames(df_final) <- rownames(x)
	return(df_final)
}

#' Create cross-species correlation matrix
#' 
#' Takes a csps object and returns a sp1 (rows) v. sp2 (columns) correlation matrix
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param use_var_genes either a character vector of genes for comparison, or logical;
#'   if TRUE, using the csps-defined top-cross genes; if FALSE, use all genes
#' @param cor_method character, corelation method to use, one of the following: 
#'   `c("pearson","spearman","kendall","jaccard","jsd","kld","wpearson","wspearman","shaferindex")` 
#'   (default: "jaccard")
#' @param gene_weights a numeric vector of non-negative weights used to calculated weighted
#'   correlations when method is "wpearson" o "wspearman". It has to match the order and length 
#    of the vector defined by the `use_var_genes` variable (i.e. either the supplied vector, 
#'   or the csps-defined vector of gene names).
#' @param pow numeric, for plotting, raise the correlation to the power of 
#'   `pow` (default: 1)
#' @param fc_thrs numeric, fold change threshold for binarising gene expression 
#'   (default: 1.2) when calculating Jaccard distance; only used when `cor_method="jaccard" or `"shaferindex"`
#' @param report_overlaps logical, whether to report shared genes between cell types (relies
# '   on the `fc_thrs` binarisation threshold
#' @param quantile_discretisation whether to convert input matrix to a discretised matrix
#'   where gene expression levels are discretised to quantiles (default: FALSE)
#' @param quantile_n if `quantile_discretisation=TRUE`, use this many quantiles (default = 10)
#' 
#' @return a list with following elements: 
#'   1) `cor_matrix` correlation matrix used for plotting
#'   2) `overlap_matrix` matrix with number of overlapping genes
#'   3) `overlapping_genes` list of overlapping genes, nested list where at 
#'   the first level are the columns from the first matrix, and at the second 
#'   level are the columns from the second matrix.
#' 
csps_correlation_matrix = function(
	csps, 
	use_var_genes = TRUE, 
	cor_method = "jaccard", 
	gene_weights = NULL,
	report_overlaps = TRUE,
	fc_thrs = 1.2, 
	pow = 1,
	add_sps_prefix = TRUE,
	prefix_sp1 = NULL,
	prefix_sp2 = NULL,
	quantile_discretisation = FALSE,
	quantile_n = 10) { 
	
	# get matrices
	mm = csps$merged
	m1 = csps$sp1
	m2 = csps$sp2

	# vector of gene names from sp1	
	gene_names = rownames(m1)
	
	if (quantile_discretisation) {
		mm = apply(mm, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
		m1 = apply(m1, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
		m2 = apply(m2, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
	}
	
	# use same rownames across all matrices (from sp1)
	rownames(mm) = gene_names
	rownames(m1) = gene_names
	rownames(m2) = gene_names

	# use unique colnames
	if (add_sps_prefix) {
		if (is.null(prefix_sp1)) prefix_sp1 = "sp1"
		if (is.null(prefix_sp2)) prefix_sp1 = "sp2"
		colnames(m1) = paste(prefix_sp1, colnames(m1), sep = "|")
		colnames(m2) = paste(prefix_sp2, colnames(m2), sep = "|")
		colnames(mm) = c(colnames(m1), colnames(m2))
	}

	# get variable genes: if use_var_genes is a vector, keep it.
	# if it's TRUE or FALSE, retrieve them from the csps object
	if (class(use_var_genes) == "logical" & use_var_genes == TRUE) {
		var_genes = csps$top_cross_sp1
		message(sprintf("csps matrix | use %i genes (from `csps` object)", length(var_genes)))
	} else if (class(use_var_genes) == "logical" & use_var_genes == FALSE) {
		var_genes = rownames(csps$merged)
		message(sprintf("csps matrix | use %i genes (complete matrix)", length(var_genes)))
	} else if (class(use_var_genes) == "character") {
		var_genes = use_var_genes
		message(sprintf("csps matrix | use %i variable genes (given)", length(var_genes)))
	} else {
		stop("`use_var_genes` has to be either logical, or a vector of genes from species 1")
	}
	
	# subset matrices to variable genes
	mm = mm [ var_genes, ]
	m1 = m1 [ var_genes, ]
	m2 = m2 [ var_genes, ]
	var_genes_sp2 = csps$og_pairs$sp2 [ csps$og_pairs$sp1 %in% var_genes ]
	
	# sanity check
	if (nrow(mm) == 0) {
		stop("No rows (genes) left in merged matrix!")
	}
	
	# list of genes for each species
	genes_sp1 = rownames(csps$sp1) [ var_genes %in% rownames(csps$sp1) ]
	genes_sp2 = rownames(csps$sp2) [ var_genes %in% rownames(csps$sp1) ]
	
	# correlation matrix calculation
	message(sprintf("csps matrix | method: %s", cor_method))
	if (cor_method == "jaccard") {
		
		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2, threshold = fc_thrs)
		
		# get jaccard distance
		com = jaccard(m1_bin, m2_bin) ^ pow
		
	} else if (cor_method == "jsd") {

		# calculate similarity based on Jensen-Shannon divergence, for a count-like matrix 
		# `est.prob = "empirical"` converts the count matrix to a probability matrix, like this:
		# counts -> 1:10
		# probs  -> 1:10 / sum(1:10)
		require("philentropy")
		com_d = philentropy::JSD(t(as.matrix(mm)), unit = "log2", est.prob = "empirical")
		# convert Jensen-Shannon distance (i.e. sqrt divergence) to similarity matrix
		com_s = as.matrix(1 - sqrt(com_d))
		colnames(com_s) = colnames(mm)
		rownames(com_s) = colnames(mm)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1) , colnames(com_s) %in% colnames(m2) ]

	} else if (cor_method == "jsdnp") {

		# calculate similarity based on Jensen-Shannon divergence, for a non-count matrix
		require("philentropy")
		com_s = philentropy::JSD(t(as.matrix(mm)), unit = "log2", est.prob = NULL)
		com_s = sqrt(com_s)
		colnames(com_s) = colnames(mm)
		rownames(com_s) = colnames(mm)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1) , colnames(com_s) %in% colnames(m2) ]

	} else if (cor_method == "kld") {

		# old approach: our own KLD function
		# com = calcKLD(mm)$mat ^ pow
		# com = com[ colnames(m1), colnames(m2) ]
		# com = 1 - com
		
		#  new approach: philentropy KLD (faster!)
		require("philentropy")
		com_d = philentropy::KL(t(as.matrix(mm)), unit = "log2", est.prob = "empirical")
		# convert KL distance (i.e. sqrt divergence) to similarity matrix
		com_s = as.matrix(1 - sqrt(com_d))
		colnames(com_s) = colnames(mm)
		rownames(com_s) = colnames(mm)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1) , colnames(com_s) %in% colnames(m2) ]
		
	} else if (cor_method == "shaferindex") {
		
		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2, threshold = fc_thrs)

		# get shaferindex distance
		com = shaferindex(m1_bin, m2_bin) ^ pow
		
	} else if (cor_method == "onls") {
		
		require("pracma")
		cod = matrix(NA, nrow = ncol(m1), ncol = ncol(m2))
		for (i in 1:ncol(m1)) {
			for (j in 1:ncol(m2)) {
				mi = log10(m1[,i])
				mj = log10(m2[,j])
				po = pracma::odregress(mi, mj)
				sq = po$ssq / nrow(m1)
				cod[i,j] = sq
			}
		}
		rownames(cod) = colnames(m1)
		colnames(cod) = colnames(m2)
		com = scale(1 / cod)

	} else if (cor_method %in% c("wpearson","wspearman")) {
		
		require("WGCNA")
		if (is.null(gene_weights)) {
			stop(sprintf("I can't use weighted correlation method %s because no weights were supplied by `gene_weights`!", cor_method))
		} else if (length(gene_weights) != length(var_genes)) {
			stop(sprintf("I can't use weighted correlation method %s because length of `gene_weights` (%i) does not match length of `var_genes` (%i)!", cor_method, length(gene_weights), length(var_genes)))
		}
		
		# create matrix of weights to apply to each footprint (based on gene-level weights applied across all tissues)
		gene_weights_m1 = matrix(rep(gene_weights, ncol(m1)), ncol = ncol(m1))
		gene_weights_m2 = matrix(rep(gene_weights, ncol(m2)), ncol = ncol(m2))
		
		# weighted correlation
		wgcna_cor_method = gsub("^w","", cor_method)
		com = WGCNA::cor(m1, m2, method = wgcna_cor_method, weights.x = gene_weights_m1, weights.y = gene_weights_m2)
		com = com ^ pow
		
	
	} else {
		
		# any other correlation metric
		com = cor(m1, m2, method = cor_method) ^ pow
		
	}

	# get shared genes between each pair of cell types
	# this is a nested named list: first level are cell types in sp1, second level are cell types in sp2
	if (report_overlaps) {

		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2, threshold = fc_thrs)
		
		# empty upper level lists (shi_l1 will contain gene names from sp1, shi_l2 for sp2)
		shi_l1 = vector("list", length = ncol(m1_bin))
		shi_l2 = vector("list", length = ncol(m1_bin))
		names(shi_l1) = colnames(m1_bin)
		names(shi_l2) = colnames(m1_bin)

		# loop through sp1 cell types
		for (i in 1:ncol(m1_bin)) {
			
			# get vector of gene presence in sp1
			v1 = m1_bin[,i]
			
			# empty lower level lists for sp2 cell types
			# create named list of cell types in sp2
			shj_l1 = vector("list", length = ncol(m2_bin))
			shj_l2 = vector("list", length = ncol(m2_bin))
			names(shj_l1) = colnames(m2_bin)
			names(shj_l2) = colnames(m2_bin)
			
			# loop through sp2 cell types
			# populate list with a vector of genes shared in cell type i and j
			for (j in 1:ncol(m2_bin)) {
				v2 = m2_bin[,j]
				shj_l1[[j]] = var_genes     [ which(v1 > 0 & v2 > 0) ]
				shj_l2[[j]] = var_genes_sp2 [ which(v1 > 0 & v2 > 0) ]
			}
			
			# populate upper level entry (sp1 cell type i) with lower level lists (sp2 cell types)
			shi_l1[[i]] = shj_l1
			shi_l2[[i]] = shj_l2
			
		}
		
		# get counts of shared genes between cell types
		sha_m = t(sapply(1:length(shi_l1), function (i) { 
			lengths(shi_l1[[i]])
		}))
		rownames(sha_m) = colnames(m1_bin)
		colnames(sha_m) = colnames(m2_bin)
		
		# shared lists
		sha_l = list(sp1 = shi_l1, sp2 = shi_l2)
		
	} else {
		
		sha_l = list(sp1 = NULL, sp2 = NULL)
		sha_m = NULL
		
	}
	
	# return
	return(list(cor_matrix = com, overlap_matrix = sha_m, overlap_genes = sha_l, method = cor_method, var_genes = var_genes))
	
}


#' Binarise gene expression matrix
#'
#' @param counts a matrix with gene counts (rows) per cell/metacell/else (columns). 
#'    By default, the appropriate threshold is obtained from the first valley in the 
#'    counts distribution (excluding zero values).
#' @param quantile numeric, a quantile of expression used to obtain a threshold for 
#'    binarisation (default is NULL, i.e. threshold is obtained from distribution).
#' @param threshold numeric, a hard threshold for binarisation (default is NULL, 
#'    i.e. threshold is obtained from distribution).
#' @param apply_on_nonzero logical, whether to apply quantiles/thresholds on non-zero
#'    values, or the whole dataset (default is TRUE)
#' 
#' @return matrix with binarised counts
#' 
sca_binarise_expression = function(
	counts,
	quantile = NULL,
	threshold = NULL,
	apply_on_nonzero = TRUE
	
) {
	
	# get distribution of nonzero values
	if (apply_on_nonzero) {
		counts_vec = counts [ counts > 0 ]
	} else {
		counts_vec = counts
	}
	
	if (is.null(threshold) & is.null(quantile)) {
		# if no quantile is specified, get threshold from first
		# valley in the counts distribution
		counts_vec_hist = hist(log10(counts_vec), plot=FALSE)
		threshold_ix = find_peaks(x = 1 / (counts_vec_hist$density + 1), m = 2)[1]
		if (is.null(threshold_ix)) {
			threshold_ix = 1
		}
		threshold = 10 ^ (counts_vec_hist$breaks [ threshold_ix ])
		message(sprintf("Binarising counts at >= %.2f threshold (first valley)", threshold))
	} else if (is.null(threshold) & is.numeric(quantile)) {
		# else, use quantile to define threshold
		threshold = quantile(counts_vec, quantile)
		message(sprintf("Binarising counts at >= %.2f threshold (quantile %.2f)", threshold, quantile))
	} else if (is.numeric(threshold)) {
		message(sprintf("Binarising counts at >= %.2f threshold (given)", threshold))
	}
	
	# apply threshold
	counts_out = (counts >= threshold) * 1
	rownames(counts_out) = rownames(counts)
	colnames(counts_out) = colnames(counts)
	
	# output
	return(counts_out)
	
}


#' Plot annotated matrix
#' 
#' @param mat any type of data matrix
#' @param name name of the type of data in the matrix (default: "data")
#' @param heatmap_colors vector of colors to map to the data
#' @param min_val,max_val min and max values of the colorscale
#' @param use_raster whether to rasterise
#' @param row_title,col_title titles for rows and columns
#' @param row_labels,col_labels ad-hoc labels for rows and columns (if set to NULL, they are taken from `mat` object)
#' @param max_length_labels truncate `row_labels` and `col_labels` to this maximum length, in characters (default 40)
#' @param fontsize size of labels (default 5 pts)
#' @param row_annot,col_annot either dataframes where the 1st column is a vector of categories for each row/column (same order is assumed) and 2nd is a vector of colors, or simply a vector of categories. Default is NULL, i.e. no annotations.
#' @param row_annot_cols,col_annot_cols named vector of colors, where names are categories that match the vector in `row_annot`/`col_annot` (not necessary if `row_annot`/`col_annot` are dataframes).
#' @param row_annot_legend,col_annot_legend whether to plot row/col annotation legends
#' @param row_cluster,col_cluster if TRUE/FALSE, whether to cluster rows/columns with default parameters. Other built-in options are "pearson", "euclidean", or a precomputed `hclust` object.
#' @param cex_dotplot transformation factor for dot size, if `do_dotplot=TRUE`
#' @param do_dotplot draw a dot plot instead of a heatmap
#' 
#' @return a ComplexHeatmap object
#' 
csps_plot_annotated_matrix = function(
	mat,
	name = "data",
	heatmap_colors = c("white","orange","orangered2","#520c52"),
	min_val = 0,
	max_val = 1,
	use_raster = TRUE,
	row_title = NULL,
	col_title = NULL,
	row_labels = NULL,
	col_labels = NULL,
	max_length_labels = 40,
	fontsize = 5,
	row_annot = NULL,
	row_annot_cols = NULL,
	row_annot_legend = FALSE,
	row_cluster = FALSE,
	col_annot = NULL,
	col_annot_cols = NULL,
	col_annot_legend = FALSE,
	col_cluster = FALSE,
	cex_dotplot = 0.02,
	do_dotplot = FALSE) { 
	
	# libraries
	require("ComplexHeatmap")
	require("circlize")
	
	# get colors for heatmap
	col_fun = circlize::colorRamp2(
		breaks = seq(from = min_val, to = max_val, length.out = length(heatmap_colors)), 
		colors = heatmap_colors)
	
	# Titles	
	# get row title
	if (is.null(row_title)) { 
		row_title = sprintf("n = %i", nrow(mat))
	} else {
		row_title = sprintf("%s, n = %i", row_title, nrow(mat))
	}
	# get col title
	if (is.null(col_title)) { 
		col_title = sprintf("n = %i", ncol(mat))
	} else {
		col_title = sprintf("%s, n = %i", col_title, ncol(mat))
	}
	
	# Row and column names
	if (is.null(row_labels)) {
		row_labels = rownames(mat)
	}
	if (is.null(col_labels)) {
		col_labels = colnames(mat)
	}
	# truncate if necessary
	if (!is.null(max_length_labels)) {
		row_labels = stringr::str_trunc(row_labels, max_length_labels)
		col_labels = stringr::str_trunc(col_labels, max_length_labels)
		row_labels = stringr::str_pad(row_labels, max_length_labels, side = "right")
		col_labels = stringr::str_pad(col_labels, max_length_labels)
	}

	# get row annotations
	# by default, get labels
	ha_row_base = ComplexHeatmap::HeatmapAnnotation(lab = anno_text(row_labels, which = "row", gp = gpar(fontsize = fontsize), just = "left"), which = "row")
	# add coloring info if available
	if (!is.null(row_annot)) {
		# ensure that row_annot is a list of dataframes
		if ("data.frame" %in% class(row_annot)) { row_annot = list(row_annot) }
		# loop through list of dataframes, adding colors
		for (ni in 1:length(row_annot)) {
			# get unique named list of colors
			row_annot_u = unique(row_annot[[ni]])
			row_annot_cols = row_annot_u[,2]
			names(row_annot_cols) = row_annot_u[,1]
			# get vector of categories
			row_annot_cats = row_annot[[ni]][,1]
			# add new annotation track
			ha_row_left  = c(ha_row_base, ComplexHeatmap::HeatmapAnnotation(clusters = row_annot_cats, col = list(clusters = row_annot_cols), which = "row", show_annotation_name = FALSE, show_legend = row_annot_legend))
			ha_row_right = c(ComplexHeatmap::HeatmapAnnotation(clusters = row_annot_cats, col = list(clusters = row_annot_cols), which = "row", show_annotation_name = FALSE, show_legend = row_annot_legend), ha_row_base)
		}
	} else {
		ha_row_left  = ha_row_base
		ha_row_right = ha_row_base
	}
		
	# get column annotations
	# by default, get labels
	ha_col_base = ComplexHeatmap::HeatmapAnnotation(lab = anno_text(col_labels, which = "column", gp = gpar(fontsize = fontsize), just = "right"), which = "column")
	# add coloring info if available
	if (!is.null(col_annot)) {
		# ensure that col_annot is a list of dataframes
		if ("data.frame" %in% class(col_annot)) { col_annot = list(col_annot) }
		# loop through list of dataframes, adding colors
		for (ni in 1:length(col_annot)) {
			# get unique named list of colors 
			col_annot_u = unique(col_annot[[ni]])
			col_annot_cols = col_annot_u[,2]
			names(col_annot_cols) = col_annot_u[,1]
			# get vector of categories
			col_annot_cats = col_annot[[ni]][,1]
			# add new annotation track
			ha_col_top = c(ha_col_base, ComplexHeatmap::HeatmapAnnotation(clusters = col_annot_cats, col = list(clusters = col_annot_cols), which = "column", show_annotation_name = FALSE, show_legend = col_annot_legend))
			ha_col_bot = c(ComplexHeatmap::HeatmapAnnotation(clusters = col_annot_cats, col = list(clusters = col_annot_cols), which = "column", show_annotation_name = FALSE, show_legend = col_annot_legend), ha_col_base)
		}
	} else {
		ha_col_top = ha_col_base
		ha_col_bot = ha_col_base
	}

	
	# should this be a dot plot?
	if (do_dotplot) {
		cell_fun_dotplot = function(j, i, x, y, width, height, fill) {
			range01 = function(x) { (x - min(x)) / (max(x) - min(x)) }
			grid.circle(
				x = x, y = y, 
				r = sqrt(range01(mat)[i, j]) * cex_dotplot, 
				gp = gpar(col = col_fun(mat[i, j]), fill = col_fun(mat[i, j]))
			)
		}
		rect_gp = gpar(type = "none")
	} else {
		cell_fun_dotplot = NULL
		rect_gp = gpar(col = "white", lwd = 0.3)
	}
	
	# how to perform row-wise clustering
	if (row_cluster %in% c("pearson","spearman","kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
		row_cluster_method = row_cluster
		row_cluster = TRUE
	} else if (row_cluster == TRUE) {
		row_cluster_method = "euclidean"
	} else {
		row_cluster = FALSE
		row_cluster_method = NULL
	}
	# how to perform column-wise clustering
	if (col_cluster %in% c("pearson","spearman","kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
		col_cluster_method = col_cluster
		col_cluster = TRUE
	} else if (col_cluster == TRUE) {
		col_cluster_method = "euclidean"
	} else {
		col_cluster = FALSE
		col_cluster_method = NULL
	}			
	
	# heatmap object
	hm = ComplexHeatmap::Heatmap(
		mat, 
		name = name,
		cell_fun = cell_fun_dotplot,
		rect_gp = rect_gp,
		use_raster = use_raster,
		col = col_fun,
		border = TRUE,
		cluster_rows = row_cluster,
		clustering_distance_rows = row_cluster_method,
		cluster_columns = col_cluster,
		clustering_distance_columns = col_cluster_method,
		row_title = row_title,
		column_title = col_title,
		show_row_names = FALSE,
		show_column_names = FALSE,
		# column annotations
		top_annotation = ha_col_top,
		bottom_annotation = ha_col_bot,
		# row annotations
		left_annotation = ha_row_left,
		right_annotation = ha_row_right
	)

	# return heatmap object
	return(hm)		
	
}


## Gene module functions ##

#' @param data data.frame of normalised gene expression values with samples as columns and rows as genes.
#' @param output_file name of output file
#' @param propGenes integer, default is 1
#' @return power estimate
#' @description A wrapper function for constructing a WGCNA network
#'
gmod_determineSoftPowerWGCNA = function(
	data,
	propGenes = 1,
	output_file = NULL, width = 8, height = 4, res = NA
) {
	
	options(stringsAsFactors = FALSE)
	
	# Remove bad genes (missing values, ...)
	nGenesInput = dim(data)[1]
	data = data[ WGCNA::goodSamplesGenes(datExpr = t(data))$goodGenes, ]
	nGenesGood = dim(data)[1]
	nGenesRemoved = nGenesInput - nGenesGood
	message(paste(nGenesRemoved, " genes filtered from dataset"))
	propGenes = round(propGenes * dim(data)[1])
	
	# Filter genes based in variance
	keepGenesExpr1 = rank( - matrixStats::rowVars(data) ) <= propGenes
	data = data[keepGenesExpr1, ]
	genesRetained = dim(data)[1]
	message(paste(genesRetained, " genes retained in dataset"))
	
	# plot powers to work out soft-power
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
	
	# Call the network topology analysis function
	sft = WGCNA::pickSoftThreshold(t(data), powerVector = powers, verbose = 5)
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		# side by side
		par(mfrow = c(1, 2))
		
		# Scale-free topology fit index as a function of the soft-thresholding power
		plot(
			sft$fitIndices[, 1],
			-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
			xlab = "Soft Threshold (power)",
			ylab = "Scale Free Topology Model Fit, signed R^2",
			col = "gray",
			main = paste("Scale independence"))
		
		text(
			sft$fitIndices[, 1],
			-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
			labels = powers, cex = 0.8, col = "blue")
		
		# this line corresponds to using an R^2 cut-off of h
		abline(h = 0.80, col = "red", lty = 2)
		
		# Mean connectivity as a function of the soft-thresholding power
		plot(
			sft$fitIndices[, 1], sft$fitIndices[, 5],
			xlab = "Soft Threshold (power)",
			ylab = "Mean Connectivity",
			col = "gray",
			main = paste("Mean connectivity"))
		
		text(
			sft$fitIndices[, 1], sft$fitIndices[, 5],
			labels = powers, cex = 0.8, col = "blue")
		
	})
	return(sft$powerEstimate)
}

#' @param data data.frame of normalised gene expression values with samples as columns and rows as genes.
#' @param propGenes A numeric value indicating the proportion of most variable genes to be retained in the analysis between 0 and 1.
#' @param softPower Integer. The soft thresholding power used to construct the WGCNA network. This can be determined using the determineSoftPowerWGCNA function.
#' @param cor_method Character. One of "spearman", "pearson", or "bicor", defines which correlation metric to calculate.
#' @param hclust_method Character. Passed to `?flashClust()`, defines which clustering method to use.
#' @param signedNetwork Logical. Keep track of the sign of the correlation in network construction?
#' @return A list of 3. data is a data.frame of the input data after filtering. geneTreeA1 is the gene tree constructed by WGCNA. dissTOMA1 is the distance matirx.
#' @description A wrapper function for constructing a WGCNA network
gmod_runWGCNA = function(
    data, 
    propGenes = 1, 
    softPower = 10, 
    cor_method = "pearson", 
    hclust_method = "average",
    signedNetwork = TRUE
) {
	
	# Try to always use unsigned network and average hclust for gene modules.
	options(stringsAsFactors = FALSE)
	
	type = ifelse(test = signedNetwork == TRUE, yes = "signed", no = "unsigned")
	# Remove bad genes (missing values, ...)
	nGenesInput = dim(data)[1]
	data = data[goodSamplesGenes(datExpr = t(data))$goodGenes, ]
	nGenesGood = dim(data)[1]
	nGenesRemoved = nGenesInput - nGenesGood
	message(paste(nGenesRemoved, " genes filtered from dataset"))
	propGenes = round(propGenes * dim(data)[1])
	
	# Filter genes based in variance
	keepGenesExpr1 = rank( -matrixStats::rowVars(data) ) <= propGenes
	data = data[keepGenesExpr1, ]
	genesRetained = dim(data)[1]
	message(paste(genesRetained, " genes retained in dataset"))
	
	# Run WGCNA on the datasets
	datExprA1g = data
	if (cor_method == "pearson") {
		message("Using pearson correlation...")
		adjacencyA1 = adjacency(t(datExprA1g), power = softPower, type = type)
	} else if (cor_method == "spearman") {
		message("Using spearman correlation...")
		adjacencyA1 = adjacency(t(datExprA1g), power = softPower, type = type, corOptions = list(use = "p", method = "spearman"))
	} else if (cor_method == "bicor") {
		message("Using weighted bicorrelation...")
		adjacencyA1 = adjacency(t(datExprA1g), power = softPower, type = type, corFnc = "bicor", corOptions = list(maxPOutliers = 0.5))
	}else{
		message("Unkown correlation function!")
		break
	}
	
	diag(adjacencyA1)=0
	dissTOMA1 = 1 - TOMsimilarity(adjacencyA1, TOMType = type)
	
	geneTreeA1 = flashClust::flashClust(as.dist(dissTOMA1), method = hclust_method)
	
	# Return the relevant objects
	return(list(data = data, geneTreeA1 = geneTreeA1, dissTOMA1 = dissTOMA1))
}

#' Plot dendogram cuts
#' @param minClusterSize Integer. Target minimum size of a gene module, see `?cutreeHybrid`. Default 10.
#' @param cutHeight Numeric. Cut dendogram at this height, see `?cutreeHybrid`. Default 0.99.
#' @param maxds Numeric. Maximum dynamic splits to plot. Default 4.
#' @param output_file
#' @param width,height,res 
gmod_plotModulesCut = function(
	referenceDataset, minClusterSize = 10, cutHeight = 0.99, maxds = 4,
	output_file = NULL, width = 12, height = 6, res = NA
) {
	
	# Seperate list objects
	mColorh = NULL
	for (ds in 0:maxds) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE, minClusterSize = minClusterSize,
			cutHeight = cutHeight, deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh = cbind(mColorh, labels2colors(tree$labels))
	}
	
	for (ds in 0:maxds) {
		tree = cutreeDynamic(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE, minClusterSize = minClusterSize,
			cutHeight = cutHeight, deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh = cbind(mColorh, labels2colors(tree))
	}
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		WGCNA::plotDendroAndColors(
			referenceDataset[["geneTreeA1"]], mColorh,
			c(paste("Hybrid_dpSplt =", 0:maxds), paste("dynamic_dpSplt =", 0:maxds)), main = "",
			dendroLabels = FALSE)
			
	})
}

gmod_plotTOMheatmap = function(
	referenceDataset,
	minClusterSize = minClusterSize, splitdepth = c(2, 3),
	output_file = NULL, width = 6, height = 6, res = NA
) {
	
	mColorh = list()
	for (ds in splitdepth) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE,
			minClusterSize = minClusterSize,
			cutHeight = cutHeight,
			deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh[[as.character(ds)]]=labels2colors(tree$labels)
	}
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		dissim = referenceDataset[["dissTOMA1"]]
		diag(dissim)=NA
		TOMplot(dissim = dissim, dendro = referenceDataset[["geneTreeA1"]], Colors = as.character(unlist(mColorh[[1]])), ColorsLeft = as.character(unlist(mColorh[[2]])))
		
	})
}


gmod_calculateModuleEigengenes = function(referenceDataset, split, minClusterSize = 10, cutHeight = 0.99) {
	#split is to be selected from examining plotModulesCut output. TBA: why default I use hybrid not dynamic?
	# Seperate list objects
	mColorh = NULL
	for (ds in 0:4) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE,
			minClusterSize = minClusterSize,
			cutHeight = cutHeight,
			deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh = cbind(mColorh, labels2colors(tree$labels))
	}
	modulesA1 = mColorh[ , split] # (Chosen based on plot)
	PCs = moduleEigengenes(t(referenceDataset$data), colors = modulesA1)
	ME = PCs$eigengenes
	rownames(ME)=colnames(referenceDataset$data)
	
	return(ME)
}

gmod_moduleHubGenes = function(referenceDataset, MEs, nGenes, split = 1) {
	
	message("ranking genes")
	kMEs = WGCNA::signedKME(datExpr = t(referenceDataset$data), datME = MEs)
	
	# rank the genes for each module on kMEs
	rankGenes= function(x) {
		kMErank = rank(-kMEs[ , x])
		genes = rownames(kMEs)
		genes = genes[order(kMErank)]
		genes[1:nGenes]
	}
	
	topGenes = lapply(1:ncol(kMEs), rankGenes)
	
	# Get the top results in a data.frame
	topGenes = do.call(cbind, topGenes)
	colnames(topGenes)=substr(colnames(kMEs), start = 4, stop = 30)
	return(topGenes)
}

gmod_moduleMembership = function(referenceDataset, MEs, kME_threshold = 0.7) {
	
	# KEY FUNCTION. Increase threshold for tighter, conservative modules.
	message("calculating kME's")
	
	kMEs = WGCNA::signedKME(datExpr = t(referenceDataset$data), datME = MEs)
	# kMEs <<- kMEs
	
	# Get gene names with kME > kME_threshold
	modGenes= function(x) {
		rownames(kMEs[kMEs[ , x] > kME_threshold, ])  # IF you want unsigned, use abs(kMEs[, x]) thresholding
	}
	gmods = lapply(X = 1:ncol(kMEs), FUN = modGenes)
	names(gmods)=substr(colnames(kMEs), start = 4, stop = 20)
	
	gene_module_membership = apply(kMEs, 1, max)
	gmods = lapply(gmods, function(x) x[order(gene_module_membership[x], decreasing = TRUE)])
	return(gmods)
	
}

gmod_moduleNonoverlapingMembership = function(referenceDataset, MEs, kME_threshold = 0.7) {
	
	# Force genes to be included ONLY in 1 module (the one with best correlation)
	message("calculating kME's")
	
	kMEs = WGCNA::signedKME(datExpr = t(referenceDataset$data), datME = MEs)
	# kMEs <<- kMEs
	
	g_to_keep = rownames(kMEs)[which(apply(kMEs, 1, max) > kME_threshold)]  # we you want unsigned, use abs(kMEs) thresholding
	kMEs = kMEs[g_to_keep, ]
	g_to_gmod = apply(kMEs, 1, which.max)  # we you want unsigned, use abs(kMEs) thresholding
	gmods = split(names(g_to_gmod), f = g_to_gmod)
	names(gmods)=substr(colnames(kMEs), start = 4, stop = 20)
	
	gene_module_membership = apply(kMEs, 1, max)
	gmods = lapply(gmods, function(x) x[order(gene_module_membership[x], decreasing = TRUE)])
	return(gmods)
	
}

gmod_moduleGeneOntology = function(gmods, go_annotation = "../Annot_gene_lists/Sros_Gene_Ontology", n_cpu = 32, p.thrs = 0.05, min_genes_in_GO = 5, p.adj = TRUE, out_folder = "./module_GO_analysis/") {
	print("I AM SLOW AS HELL, DAMN YOU TopGO!!")
	thread_num = n_cpu
	doMC::registerDoMC(thread_num)
	
	go = readMappings(go_annotation)
	dir.create(out_folder, showWarnings = FALSE)
	
	go_analysis= function(i) {
		
		fg_genes = as.factor(as.integer(names(go) %in% gmods[[i]]))
		names(fg_genes)=names(go)
		
		GOdataBP = new("topGOdata", ontology = c("BP"), allGenes = fg_genes, annot = annFUN.gene2GO, gene2GO = go)
		fisher = runTest(GOdataBP, "classic", "fisher")
		parent_child = runTest(GOdataBP, "parentchild", "fisher")
		#suppressMessages(ks = runTest(GOdata, "classic", "ks"))
		r1 = GenTable(GOdataBP, fisher = fisher, parent_child = parent_child, topNodes= length(usedGO(GOdataBP)), orderBy = "fisher")
		r1 = r1[which(r1$Annotated > min_genes_in_GO), ]
		r1$fisher = as.numeric(r1$fisher)
		r1$parent_child = as.numeric(r1$parent_child)  #retarded TopGO reports values as characters...
		if (p.adj) {
			r1$fisher = p.adjust(r1$fisher, method = "BH")
			r1$parent_child = p.adjust(r1$parent_child, method = "BH")
		}
		r1[is.na(r1)]=1e-30 # TopGO doesn't report smaller evals
		r1 = r1[which(r1$parent_child < p.thrs | r1$fisher < p.thrs), ]
		r1 = cbind.data.frame(r1[, 1:2], rep("BP", nrow(r1)), r1[, 3:ncol(r1)])
		colnames(r1)[3]="GO_cat"
		
		GOdataMF = new("topGOdata", ontology = c("MF"), allGenes = fg_genes, annot = annFUN.gene2GO, gene2GO = go)
		fisher = runTest(GOdataMF, "classic", "fisher")
		parent_child = runTest(GOdataMF, "parentchild", "fisher")
		#suppressMessages(ks = runTest(GOdata, "classic", "ks"))
		r2 = GenTable(GOdataMF, fisher = fisher, parent_child = parent_child, topNodes= length(usedGO(GOdataMF)), orderBy = "fisher")
		r2 = r2[which(r2$Annotated > min_genes_in_GO), ]
		r2$fisher = as.numeric(r2$fisher)
		r2$parent_child = as.numeric(r2$parent_child)
		if (p.adj) {
			r2$fisher = p.adjust(r2$fisher, method = "BH")
			#r2$parent_child = p.adjust(r2$parent_child, method = "BH")
		}
		r2[is.na(r2)]=1e-30 # TopGO doesn't report smaller evals.
		r2 = r2[which(r2$parent_child < p.thrs | r2$fisher < p.thrs), ]
		r2 = cbind.data.frame(r2[, 1:2], rep("MF", nrow(r2)), r2[, 3:ncol(r2)])
		colnames(r2)[3]="GO_cat"
		
		r_final = rbind.data.frame(r1, r2)
		r_final = r_final[order(r_final$fisher), ]
		r_final$fisher = formatC(r_final$fisher, format = "e", digits = 2) #let's reports few dec positions
		r_final$parent_child = formatC(r_final$parent_child, format = "e", digits = 2)
		write.table(r_final, file = paste0(out_folder, "module_", i, "_GO.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	}
	
	plyr::alply(names(gmods), 1, function(x) go_analysis(x), .parallel = TRUE)
	
}

#' function to visualize WGCNA modules on metacells and to export annotated gene module table
#' 
#' @param expr_matrix expression matrix (e.g. `mc@mc_fp`)
#' @param gmods expression matrix
#' @param me eigenvalue matrix
#' @param expr_matrix_colors vector of metacell colors (same order as mc@mc_fp)
#' @param an_output_file,me_output_file,ex_output_file output files for the gene module annotations, the eigen values heatmap, and the expression heatmap
#' @param me_width,me_height size of the eigenvalues heatmap
#' @param do_expression boolean, whether to output the huge expression map (it's usually unadvisable, it's too big)
#' @param ex_width,ex_height size of the expression heatmap
#' @param resolution_rate downsample matrix to this rate, by average matrix subsampling
#' @param highlight_genes vector of genes to highlight in the expression heatmap
#' @param highlight_genes_annot vector of annotations for the genes to highlight in the expression heatmap (same order)
#' @param cor_cutoff_min,cor_cutoff_max low and high values for the correlation heatmap (default: 0 and 0.5). NULL sets max to quantile 0.99
#' @param eigen_min,eigen_ax low and high values for the eigenvalue heatmap (default: -0.1 to 0.2). NULL sets max to quantile 0.98
#' @param heatmap_colors,heatmap_colors_cor vector of colors to map to the heatmaps
#' 
#' @return nothing, only output plots
#' 
gmod_plotMCheatmap_annotate_modules = function(
	expr_matrix,
	gmods,
	me,
	expr_matrix_colors = NULL,
	an_output_file = "wgcna.gmod_annotation.csv", 
	me_output_file = "wgcna.gmod_eigenvalues.pdf", 
	ex_output_file = "wgcna.gmod_expression.pdf", 
	do_expression = TRUE,
	me_width = 10,  me_height = 5,
	ex_width = 20, ex_height = 10,
	res = NA,
	resolution_rate = 0.1,
	annotation = NULL,
	highlight_genes = NULL, 
	highlight_genes_annot = NULL,
	cor_cutoff_min = 0, 
	cor_cutoff_max = 0.5,
	eigen_min = -0.1,
	eigen_max = 0.2,
	heatmap_colors = c("white","orange","orangered2","#520c52"),
	heatmap_colors_cor = c("white","#d6e72e","#6fb600","#003f4d")
	
) {

	# me derive from calculateModuleEigengenes
	me = t(me) 
	
	# plot eigenvalues and reorder gmods based on them
	me = me[order(apply(me, 1, function(x) which.max(rollmean(x, 1)))), ]
	
	# do we have colors for the metacells?
	if (!is.null(expr_matrix_colors)) {
		mc_colors = expr_matrix_colors
		names(mc_colors) = 1:length(mc_colors)
	} else {
		mc_colors = FALSE
	}
	
	# plot
	# browser()
	plotting_function(me_output_file, width = me_width, height = me_height, res, EXP = {
		
		# color vector
		me_plot = me [ rev(rownames(me)) , ]
		gmod_color_me = gsub("^ME", "", rownames(me_plot))
		names(gmod_color_me) = gsub("^ME", "", rownames(me_plot))
		
		# min and max values for the eigen heatmap
		if (is.null(eigen_max)) {
			eigen_max = quantile(me_plot, 0.98)
		}
		if (is.null(eigen_min)) {
			eigen_min = min(me_plot)
		}
		
		# plot
		hm = gmod_plot_complex_heatmap(
			mat = me_plot,
			name = "eigenvalue",
			color_mat = heatmap_colors_cor,
			color_min = eigen_min, 
			color_max = eigen_max,
			cluster_row = FALSE, 
			cluster_col = FALSE, 
			colors_row = gmod_color_me,
			colors_col = mc_colors,
			fontsize = 5, 
			use_raster = FALSE,
			title_row = "gene modules")
			
		print(hm)
		
	})
	gmods = gmods[gsub("ME", "", rownames(me))]
	gmods = sapply(gmods, function(x) x[x %in% rownames(expr_matrix)], USE.NAMES = TRUE, simplify = FALSE)
	unlist_gmods = unlist(gmods)
	
	# sort genes inside gmod by membership score?
	gene_module_membership = apply(WGCNA::signedKME(t(expr_matrix[unlist_gmods, ]), t(me)), 1, max)
	gmods = lapply(gmods, function(x) x[order(gene_module_membership[x], decreasing = TRUE)])

	# save gene module annotation table	
	tab = data.frame(
		gene = unlist_gmods, 
		gene_module = rep(names(gmods), lengths(gmods)), 
		membership_score = round(gene_module_membership[unlist_gmods], 3)
	)
	# add gene annotations, if available
	if (!is.null(annotation)) {
		tab = cbind(tab, annotation[unlist(gmods), ])
		colnames(tab) [ ( ncol(tab) - ncol(annotation) + 1 ) : ncol(tab) ] = colnames(annotation)
	}
	write.table(tab, file = an_output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	
	# compute pearson correlation matrix (^ 5) and scale it
	if (do_expression) {
		x_cor = tgs_cor(t(expr_matrix[unlist(gmods), ]), spearman = FALSE)
		x_cor[is.na(x_cor)]=0
		diag(x_cor)=0
		x_cor = (( 1 + x_cor ) / 2) ^ 5
	
		plotting_function(ex_output_file, width = ex_width, height = ex_height, res, EXP = {
			
			# subsample x_cor matrix
			if (resolution_rate < 1) {
				x_cor_s = matrix_subsample(x_cor, subsample_rate = resolution_rate)
			} else {
				x_cor_s = x_cor
			}
			rownames(x_cor_s) = 1:nrow(x_cor_s)
			colnames(x_cor_s) = 1:ncol(x_cor_s)
			
			# subsample gmod
			g_cor_v = lengths(gmods)
			g_cor_s = g_cor_v * resolution_rate
			g_cor_s = round_smart(g_cor_s)
			gmod_color = names(g_cor_s)
			names(gmod_color) = names(g_cor_s)
			gmod_gene_color_category = rep(names(g_cor_s), g_cor_s)
			names(gmod_gene_color_category) = unlist(gmods)
			
			# colnames are names of gmods (located at the start position of each module)
			v = as.numeric(as.factor(gmod_color))
			ixs_names = c(1, 1 + which(diff( v ) != 0))
			ixs_other = c(1 + which(diff( v ) == 0))
			ixs_names = pmin(ixs_names, dim(x_cor_s)[1])
			ixs_names = pmax(ixs_names, 1)
			ixs_other = pmin(ixs_other, dim(x_cor_s)[1])
			ixs_other = pmax(ixs_other, 1)
			colnames(x_cor_s) [ ixs_names ] = names(g_cor_s)
			colnames(x_cor_s) [ ixs_other ] = ""
			
			# min and max values for the correlation map
			if (is.null(cor_cutoff_max)) {
				cor_cutoff_max = quantile(x_cor_s, 0.99)
			}
			if (is.null(cor_cutoff_min)) {
				cor_cutoff_min = 0
			}
			
			# first heatmap: gene-gene correlation
			hm1 = gmod_plot_complex_heatmap(
				mat = x_cor_s,
				name = "cor^5",
				color_mat = heatmap_colors_cor,
				color_min = cor_cutoff_min,
				color_max = cor_cutoff_max,
				cluster_row = FALSE,
				cluster_col = FALSE,
				name_row_show = FALSE,
				name_col_show = TRUE,
				colors_row = gmod_color,
				colors_col = gmod_color,
				fontsize = 5,
				use_raster = TRUE,
				raster_quality = 1,
				title_col = "genes sorted by gene module",
				title_row = sprintf("%i genes in %i modules", nrow(x_cor), length(gmods)))
			
			# second heatmap: footprints
			if (resolution_rate < 1) {
				expr_s = matrix_subsample(expr_matrix[rownames(x_cor), ], subsample_rate_row = resolution_rate, subsample_rate_col = 1)
			} else {
				expr_s = expr_matrix[rownames(x_cor), ]
			}
			rownames(expr_s) = rep("", nrow(expr_s))
			colnames(expr_s) = colnames(expr_matrix)
			
			# rownames are names of selected genes (if needed)
			if (!is.null(highlight_genes)) {
				ixs_tfs_o = which(rownames(x_cor) %in% highlight_genes)
				vec_tfs_o = rownames(x_cor) [ ixs_tfs_o ]
				ixs_tfs_s = round(ixs_tfs_o * resolution_rate)
				ixs_tfs_s = pmin(ixs_tfs_s, dim(expr_s)[1])
				ixs_tfs_s = pmax(ixs_tfs_s, 1)
				ixs_oth_s = which(!as.character(1:nrow(expr_s)) %in% as.character(ixs_tfs_s))
				ixs_oth_s = pmin(ixs_oth_s, dim(expr_s)[1])
				ixs_oth_s = pmax(ixs_oth_s, 1)
				rownames(expr_s) [ ixs_tfs_s ] = vec_tfs_o
				rownames(expr_s) [ is.na(rownames(expr_s)) ] = ""
				
				if (!is.null(highlight_genes_annot)) {
					names(highlight_genes_annot) = highlight_genes
					vec_ann_o = highlight_genes_annot [ rownames(expr_s) [ ixs_tfs_s ] ]
					rownames(expr_s) [ ixs_tfs_s ] = paste(vec_tfs_o, stringr::str_trunc(vec_ann_o, width = 40), sep = " | ")
				}
				
			}
			
			hm2 = gmod_plot_complex_heatmap(
				mat = expr_s,
				name = "fp",
				color_mat = heatmap_colors,
				color_min = 0,
				color_max = 4,
				cluster_row = TRUE,
				cluster_col = FALSE,
				name_row_show = TRUE,
				name_col_show = TRUE,
				colors_row = gmod_color,
				colors_col = mc_colors,
				fontsize = 3.5,
				use_raster = TRUE,
				raster_quality = 1,
				title_row = "genes")
			
			# relative sizes of each heatmap
			hm1@matrix_param$width = unit(0.7, "npc")
			hm2@matrix_param$width = unit(0.5, "npc")

			print(hm1 + hm2)

		})
	}
}


gmod_geneModuleColors = function(
	referenceDataset, 
	split,
	minClusterSize = 20) {
	# Seperate list objects
	mColorh = NULL
	for (ds in 0:3) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE, minClusterSize = minClusterSize,
			cutHeight = 0.99, deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]]
		)
		mColorh = cbind(mColorh, labels2colors(tree$labels))
	}
	
	modules = mColorh[ , split] # (Chosen based on plot)
	return(modules)
}

gmod_getModuleLabels = function(
	referenceDataset,
	split,
	minClusterSize = 20) {
	tree = cutreeHybrid(
		dendro = referenceDataset[["geneTreeA1"]],
		pamStage = FALSE, minClusterSize = minClusterSize,
		cutHeight = 0.99, deepSplit = split,
		distM = referenceDataset[["dissTOMA1"]]
	)
	return(tree)
}

#compares the preservation of WGCNA solution in another dataset.
gmod_preservationStatsWGCNA = function(
	referenceDataset, data2, split = 3, type = "unsigned",
	nPermutations = 100, minClusterSize = 20,
	greyName = "grey", qVal = FALSE
) {
	
	colors = geneModuleColors(referenceDataset, split, minClusterSize)
	# Remove bad genes (missing values, ...)
	data2 = data2[goodSamplesGenes(datExpr = t(data2))$goodGenes, ]
	
	geneModules = cbind(rownames(referenceDataset$data), colors)
	
	# Quantify module preservation
	data = referenceDataset[["data"]]
	multiExpr = list(A1 = list(data = t(data)), A2 = list(data = t(data2)))
	multiColor = list(A1 = colors)
	
	mp = modulePreservation(
		multiData = multiExpr, multiColor = multiColor,
		referenceNetworks = 1, verbose = 3,
		calculateQvalue = qVal,
		networkType = type,
		nPermutations = nPermutations,
		maxGoldModuleSize = 100,
		greyName = greyName
		#maxModuleSize = 400
	)
	
	stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
	stats = stats[order(-stats[, 2]), ]
	print(stats[ , 1:3])
	return(list(mp = mp, geneModules = geneModules))
}


# Convert module colors to module names (so thing look more professional)
gmod_convertColorToLabel = function(colors, n = 100, prefix = "mod") {
	
	# Get the WGCNA colors
	cols = standardColors(n)
	labels = 1:n
	labels = paste(prefix, labels, sep = "")
	colToLabel = data.frame(colors = cols, labels = labels)
	
	# Add grey as label 0
	grey = c("grey", paste(prefix, "0", sep = ""))
	colToLabel = rbind(grey, colToLabel)
	
	# Function to return corresponding label for a given color
	convertFunc= function(x) {
		colToLabel[colToLabel$colors == x, ]$labels
	}
	
	# Apply the function to a list of colors
	unlist(lapply(colors, convertFunc))
}

gmod_compare_module_solutions = function(mods1, mods2) {
	# input must be a vector of module membership in same order.
	f = overlapTable(mods1, mods2)
	return(f)
}


gmod_compare_gene_module_membership = function(list1, list2) {
	common = intersect(unlist(list1), unlist(list2))
	
	perc_list1 = round(length(common) / length(unlist(list1)), 3)
	perc_list2 = round(length(common) / length(unlist(list2)), 3)
	
	message("Fraction of list1 genes in common: ", perc_list1)
	message("Fraction of list2 genes in common: ", perc_list2)
	
	list1_red = lapply(list1, function(x) intersect(x, common))
	list2_red = lapply(list2, function(x) intersect(x, common))
	m1 = matrix(rep(0, length(common) * length(common)), ncol = length(common), nrow = length(common))
	colnames(m1)=rownames(m1)=unlist(list1_red)
	m2 = matrix(rep(0, length(common) * length(common)), ncol = length(common), nrow = length(common))
	colnames(m2)=rownames(m2)=unlist(list2_red)
	
	# sapply(names(list1_red), function(x) m1[unlist(list1_red[[x]]), unlist(list1_red[[x]])]=1)
	# sapply(names(list2_red), function(x) m2[unlist(list2_red[[x]]), unlist(list2_red[[x]])]=1)
	
	for (module in names(which(lengths(list1_red) > 0))) {
		m1[list1_red[[module]], list1_red[[module]]]=1
	}
	
	
	for (module in names(which(lengths(list2_red) > 0))) {
		m2[list2_red[[module]], list2_red[[module]]]=1
	}
	
	m1_focused = m1
	m2_reord = m2[unlist(list1_red), unlist(list1_red)]
	m1_focused[lower.tri(m1_focused)]=m2_reord[lower.tri(m2_reord)]
	
	m2_focused = m2
	m1_reord = m1[unlist(list2_red), unlist(list2_red)]
	m2_focused[upper.tri(m2_focused)]=m1_reord[lower.tri(m1_reord)]
	
	png("test.pdf", height = 4, width = 4)
	par(mar = c(0, 0, 0, 0))
	par(fig = c(0.1, 0.5, 0.1, 0.95))
	pheatmap(m1_focused, cluster_cols = FALSE, cluster_rows = FALSE)
	
	par(fig = c(0.55, 0.9, 0.1, 0.95))
	pheatmap(m2_focused, cluster_cols = FALSE, cluster_rows = FALSE)
	dev.off()
}







# Generic function to plot complex heatmaps
gmod_plot_complex_heatmap = function(
		mat,
		name = "heatmap",
		color_mat = c("white","#d6e72e","#6fb600","#003f4d"),
		color_min = 0,
		color_max = 1,
		fontsize = 10,
		categories_col = NULL,
		categories_row = NULL,
		separate_col = FALSE,
		separate_row = FALSE,
		colors_col = NULL,
		colors_row = NULL,
		title_row = NULL,
		title_col = NULL,
		name_row_show = TRUE,
		name_col_show = TRUE,
		cluster_row = TRUE,
		cluster_col = TRUE,
		use_raster = TRUE,
		raster_quality = 1,
		show_legend_row = FALSE,
		show_legend_col = FALSE,
		both_sides_row = TRUE,
		both_sides_col = TRUE,
		cell_border = gpar(col = NA, lwd = 1, lty = 1),
		heatmap_border = gpar(col = NA, lwd = 1, lty = 1),
		do_dotplot = FALSE,
		dot_size_mat = NULL,
		dot_size_min = NULL,
		dot_size_max = NULL,
		cex_dotplot = 0.02
) {
	
	require("ComplexHeatmap")
	require("circlize")
	ht_opt$message = FALSE
	
	# color function
	col_fun = circlize::colorRamp2(seq(color_min, color_max, length.out = length(color_mat)), color_mat)
	
	
	# # vector of clusters (row/col annotation)
	# cluster_vector = as.character(mot_merge_d$cluster)
	# # categorical colors for clusters (get recycled)
	# catcol_vec = rep(catcol_lis, length.out = length(unique(cluster_vector)))
	# names(catcol_vec) = unique(cluster_vector)
	# cluster_colors = catcol_vec [ cluster_vector ]
	# names(cluster_colors) = cluster_colors
	# # quantitative colorscale for heatmap
	# catcol_fun = circlize::colorRamp2(1:length(catcol_vec), catcol_vec)
	# names(cluster_colors) = cluster_vector
	
	if (is.null(title_row)) {
		title_row = sprintf("n = %i rows", nrow(mat))
	}
	if (is.null(title_col)) {
		title_col = sprintf("n = %i columns", ncol(mat))
	}
	
	# left row annotations
	if (is.null(categories_row)) {
		categories_row = rownames(mat)
	}
	if (is.null(colors_row)) {
		colors_row = rep(NA, nrow(mat))
		left_annotation = NULL
	} else {
		if (is.null(names(colors_row))) { names(colors_row) = categories_row }
		left_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_row, name = title_row, col = list(c = colors_row), which = "row", show_legend = show_legend_row)
	}
	
	# top col annotations
	if (is.null(categories_col)) {
		categories_col = colnames(mat)
	}
	if (is.null(colors_col)) {
		colors_col = rep(NA, ncol(mat))
		top_annotation = NULL
	} else {
		if (is.null(names(colors_col))) { names(colors_col) = categories_col }
		top_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_col, name = title_col, col = list(c = colors_col), which = "column", show_legend = show_legend_col)
	}
	
	# add annotations to both sides of the rows or columns?
	if (both_sides_col) {
		bottom_annotation = top_annotation
	} else {
		bottom_annotation = NULL
	}
	if (both_sides_row) {
		right_annotation = left_annotation
	} else {
		right_annotation = NULL
	}
	
	# split columns and rows?
	if (separate_row) {
		split_row = categories_row
	} else {
		split_row = NULL
	}
	if (separate_col) {
		split_col = categories_col
	} else {
		split_col = NULL
	}
	
	# should this be a dot plot?
	if (do_dotplot) {
		if (is.null(dot_size_mat)) {
			dot_size_mat = mat
		}
		if (!is.null(dot_size_max) & !is.null(dot_size_min)) {
			# dot_size_mat [ dot_size_mat > dot_size_max ] = dot_size_max
			# dot_size_min [ dot_size_mat < dot_size_min ] = dot_size_min
			rangesize = function(x) { (x - dot_size_min) / (dot_size_max - dot_size_min) }
		} else {
			rangesize = function(x) { (x - min(x)) / (max(x) - min(x)) }
		}
		cell_fun_dotplot = function(j, i, x, y, width, height, fill) {
			grid.circle(
				x = x, y = y, 
				r = sqrt(rangesize(dot_size_mat)[i, j]) * cex_dotplot, 
				gp = gpar(col = NA, fill = col_fun(mat[i, j]))
			)
		}
		cell_border = gpar(type = "none")
		# cell_border = gpar(col = "gray", lwd = 1, lty = 1, fill = NULL)
	} else {
		cell_fun_dotplot = NULL
	}
	
	
	# plot
	hm = ComplexHeatmap::Heatmap(
		mat,
		name = name,
		cell_fun = cell_fun_dotplot,
		use_raster = use_raster,
		raster_quality = raster_quality,
		cluster_rows = cluster_row,
		cluster_columns = cluster_col,
		row_title = title_row,
		column_title = title_col,
		show_row_names = name_row_show,
		show_column_names = name_col_show,
		row_names_gp = gpar(fontsize = fontsize),
		column_names_gp = gpar(fontsize = fontsize),
		top_annotation = top_annotation,
		left_annotation = left_annotation,
		right_annotation = right_annotation,
		bottom_annotation = bottom_annotation,
		column_split = split_col,
		row_split = split_row,
		row_gap = unit(0.5, "mm"),
		column_gap = unit(0.5, "mm"),
		rect_gp = cell_border,
		border_gp = heatmap_border,
		col = col_fun)
	
	# return heatmap
	return(hm)
	
}

# Helper function to round values without changing total sum
round_smart <- function(x, digits = 0) {
	up <- 10 ^ digits
	x <- x * up
	y <- floor(x)
	indices <- tail(order(x - y), round(sum(x)) - sum(y))
	y[indices] <- y[indices] + 1
	y / up
}


# subsample matrix at a certain rate (to reduce resolution of a matrix in a brute-force way)
matrix_subsample = function(mat, subsample_rate = 0.1, subsample_rate_row = NULL, subsample_rate_col = NULL) {
	
	# subsample rates
	if (is.null(subsample_rate_row)) { 
		subsample_rate_row = subsample_rate
	}
	if (is.null(subsample_rate_col)) { 
		subsample_rate_col = subsample_rate
	}
	
	# num bins and widths
	nbin_r = nrow(mat) * subsample_rate_row
	nbin_c = ncol(mat) * subsample_rate_col
	wbin_r = round(nrow(mat) / nbin_r)
	wbin_c = round(ncol(mat) / nbin_c)
	
	# init subsampled matrix
	smat = matrix(nrow = nbin_r, ncol = nbin_c)
	
	# log
	message(sprintf("subsample matrix at %.2f x %.2f rate: %i x %i to %i x %i", subsample_rate_row, subsample_rate_col, nrow(mat), ncol(mat), nrow(smat), ncol(smat)))
	
	# loop
	for (rs in 0:(nbin_r - 1)) {
		ro_s = (rs * wbin_r) + 1
		ro_e = rs * wbin_r + wbin_r
		for (cs in 0:(nbin_c - 1)) {
			co_s = (cs * wbin_c) + 1
			co_e = cs * wbin_c + wbin_c
			smat_i = mat [ ro_s : ro_e , co_s : co_e ]
			smat[rs + 1, cs + 1]   = mean(smat_i)
		}
	}
	
	# out
	return(smat)
	
}


#' Map modules to metacells based on eigenvalue matrix
#' 
#' @param me eigenvalue matrix
#' @param ct_vector a named vector: names are metacells (corresponding to the rownames of the `me` matrix) and contents are cell types (default is NULL, i.e. the function simply maps the top mc to each module, but doesn't report cell types)
#' @param clean_module_names whether to remove the `ME` prefix in the module names (default is TRUE)
#' 
#' @return dataframe with module to mc and cell type annotations
#' 
gmod_annotate_modules_to_mc_and_ct = function(me, ct_vector = NULL, clean_module_names = TRUE, output_fn = NULL) {
	
	# clean module names? (remove ME prefix)
	if (clean_module_names) {
		colnames(me) = gsub("^ME","",colnames(me))
	}
	
	# find best mc per module
	top_mc_per_module = rownames(me) [ apply(me, 2, which.max) ]
	names(top_mc_per_module) = colnames(me)

	# if cell type table is given, map mcs to cell types
	if (!is.null(ct_vector)) {
		top_ct_per_module = ct_vector [ top_mc_per_module ]
	} else {
		top_ct_per_module = rep(NA, length(top_mc_per_module))
	}

	# output
	module_annots = data.frame(
		module = names(top_mc_per_module),
		mc_top = top_mc_per_module,
		ct_top = top_ct_per_module
	)
	
	if (!is.null(output_fn)) {
		write.table(module_annots, output_fn, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
	
	return(module_annots)
	
}