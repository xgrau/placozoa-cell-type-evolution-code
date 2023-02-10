# libraries
library("metacell")
source("../scripts/Downstream_functions.R")
source("../scripts/Modified_functions.R")
source("../scripts/helper.R")
par(family  = "Arial")


# where to store output?
out_fn = "results_metacell_it4/"
dir.create(out_fn, showWarnings = FALSE)
override_params("data/metacell_parameters.yaml","metacell")

# list of species
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

for (spi in sps_list) {
	
	# init database
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	
	# first load old mc solution
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it3",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it2",run_name))
	
	# load cell type annotations for it4
	ctt_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.tsv", spi)
	ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")
	
	# update order with it4 annotations 
	ctt$metacell_it2 = ctt$metacell
	ctt$metacell = 1:nrow(ctt)
	mc_r = metacell::mc_reorder(mc = mc, ord = ctt$metacell_it2)
	mc_r@colors = ctt$color
	names(mc_r@colors) = ctt$metacell
	ctt_new = sprintf("results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
	write.table(ctt, ctt_new, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

	# save reordered metacells as it4
	metacell::scdb_add_mc(sprintf("%s_it4",run_name), mc = mc_r)
	mc = metacell::scdb_mc(sprintf("%s_it4",run_name))
	
	# COUNTS
    # metacell size and counts
    mc_stat = scp_plot_mc_size_counts(
        mc_object = mc, 
        mc_color = mc@colors, 
        mat_object = mat, 
        output_file = sprintf("%s/%s.mc_size_counts.pdf", out_fn, run_name),
		width = 32,
        return_table = TRUE)
    write.table(mc_stat, sprintf("%s/%s.mc_size_counts.csv", out_fn, run_name), quote = FALSE, sep = "\t", row.names = FALSE)
    gc()
    
    # check batch effects
    scp_batch_effects_per_mc(
        mc, 
        mat, 
        mc_color = mc@colors, 
        output_file = sprintf("%s/%s.batch_effects_mc.pdf", out_fn, run_name))

    # clean
    gc()

	
	# 2D PROJECTION
	
	# create projection
	metacell::mcell_mc2d_force_knn(
		mc2d_id = sprintf("%s_it4",run_name),
		mc_id = sprintf("%s_it4",run_name), 
		feats_gset = sprintf("%s_it2",run_name), 
		graph_id = sprintf("%s_it2",run_name))
	
	# load projection
	mc2d = metacell::scdb_mc2d(sprintf("%s_it4",run_name))
	
	# plot 2D projection of reclustered subset
	scp_plot_sc_2d(
		mc2d = mc2d,
		mc = mc,
		plot_mcs=TRUE,
		plot_mc_name=TRUE,
		plot_edges=TRUE,
		output_file=sprintf("%s/%s.2dproj.pdf", out_fn, run_name),
		width=8, height=8)

	# plot 2D projection, only metacells 
	scp_plot_sc_2d(
		mc2d = mc2d,
		mc = mc,
		plot_mcs=TRUE,
		cex_sc = 0,
		plot_mc_name=TRUE,
		plot_edges=TRUE,
		output_file=sprintf("%s/%s.2dproj_mcs.pdf", out_fn, run_name),
		width=8, height=8)
	
	# plot 2D projection of reclustered subset
	scp_plot_sc_2d(
		mc2d = mc2d,
		mc = mc,
		plot_mcs=TRUE,
		plot_mc_name=TRUE,
		plot_edges=TRUE,
		output_file=sprintf("%s/%s.2dproj_large.pdf", out_fn, run_name),
		width=16, height=16)
	
	# load gene annotations (based on transcripts, need to be changed to genes)		
	gene_annot = read.table(sprintf("../data/reference/%s_long.pep.annotations.csv", spi), sep = "\t", row.names = 1)
	rownames(gene_annot) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = rownames(gene_annot))
	# load TF-specific annotations
	dic_genetf = read.table(sprintf("../data/gene_annotations/tfs.%s_genes.curated.csv", spi), sep = "\t", row.names = 1, col.names = c("gene","OG"))
	dic_genetf_gene_name = gsub("like:","like_",dic_genetf$OG)
	dic_genetf_gene_name = stringr::str_split(dic_genetf_gene_name, pattern = ":", simplify = TRUE)[,2]
	rownames(dic_genetf) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = rownames(dic_genetf))
	names(dic_genetf_gene_name) = rownames(dic_genetf)
	gene_annot_tfs = merge(gene_annot, dic_genetf_gene_name, by.x = 0, by.y = 0, all.x = TRUE, all.y = FALSE)[,"y"]
	gene_annot [ !is.na(gene_annot_tfs) , "V2" ] = gene_annot_tfs [ !is.na(gene_annot_tfs) ]
	
	
	# HEATMAP TABLES
	# first, get marker data (select genes, colors, etc.)
	marker_data_list = scp_plot_cmod_markers_select(
		mc_object = mc,
		gene_annot_file = gene_annot,
		per_clust_genes = 20,
		gene_min_fold = 1.5,
		gene_font_size = 5,
		clust_ord = colnames(mc@mc_fp))
	
	# draw metacell heatmap
	scp_plot_cmod_markers_mc(
		marker_data_list = marker_data_list,
		output_file = sprintf("%s/%s.exp_global_mc.pdf", out_fn, run_name),
		heatmap_colors = c("gray98","orange","orangered2","#520c52"),
		clust_col = mc@colors,
		width = max(8, round(length(marker_data_list$clust_ord) / 10)),
		height = round(length(marker_data_list$genes) / 15) + 1,
		show_gene_names = TRUE,
		highlight_genes = names(dic_genetf_gene_name),
		gene_chr_limit = 60,
		show_clust_borders = TRUE,
		use_raster=TRUE,
		max_expression_fc = 4,
		min_expression_fc = 1
	)
	
	# compressed metacell map
	scp_plot_cmod_markers_mc(
		marker_data_list = marker_data_list,
		output_file = sprintf("%s/%s.exp_global_mc-c.pdf", out_fn, run_name),
		heatmap_colors = c("gray98","orange","orangered2","#520c52"),
		clust_col = mc@colors,
		width = 10,
		height = 18,
		show_gene_names = FALSE,
		highlight_genes = names(dic_genetf_gene_name),
		gene_chr_limit = 60,
		show_clust_borders = TRUE,
		use_raster = TRUE,
		max_expression_fc = 4,
		min_expression_fc = 1
	)
	
	# compressed metacell map
	scp_plot_cmod_markers_mc(
		marker_data_list = marker_data_list,
		output_file = sprintf("%s/%s.exp_global_mc-c2.pdf", out_fn, run_name),
		# heatmap_colors = c("gray98","orange","orangered2","#520c52"),
		heatmap_colors = c("gray98","#ceaed3","#b57dbe","#503b92"),
		clust_col = mc@colors,
		width = 10,
		height = 18,
		show_gene_names = FALSE,
		highlight_genes = names(dic_genetf_gene_name),
		gene_chr_limit = 60,
		show_clust_borders = TRUE,
		use_raster = TRUE,
		max_expression_fc = 4,
		min_expression_fc = 1
	)
	
	# compressed metacell map
	scp_plot_cmod_markers_mc(
		marker_data_list = marker_data_list,
		output_file = sprintf("%s/%s.exp_global_mc-c3.pdf", out_fn, run_name),
		heatmap_colors = c("gray99","#accbcc","#508490","#004066","#000738"),
		clust_col = mc@colors,
		width = 10,
		height = 18,
		show_gene_names = FALSE,
		highlight_genes = names(dic_genetf_gene_name),
		gene_chr_limit = 60,
		show_clust_borders = TRUE,
		use_raster = TRUE,
		max_expression_fc = 4,
		min_expression_fc = 1
	)
	
	# draw single cell heatmap
	scp_plot_cmod_markers_sc(
		marker_data_list = marker_data_list,
		mc_object = mc,
		mat_object = mat,
		output_file = sprintf("%s/%s.exp_global_sc.pdf", out_fn, run_name),
		heatmap_colors = c("gray98","orange","orangered2","#520c52"),
		clust_col = mc@colors,
		width = max(12, round(length(marker_data_list$clust_ord) / 2)),
		height = round(length(marker_data_list$genes) / 15) + 1,
		show_gene_names = TRUE,
		highlight_genes = names(dic_genetf_gene_name),
		smoothen = 5,
		gene_chr_limit = 60,
		show_clust_borders=TRUE,
		use_raster=TRUE,
		max_expression_fc = 4,
		min_expression_fc = 1
	)
	
	# get mc counts, umifrac
	mc_counts = sca_mc_gene_counts(mc_object = mc, mat_object = mat, T_totumi = 0)
	mc_umifrac = sca_mc_gene_umifrac(mc, mc_counts)
	
	# save mc-level tables
	write.table(mc_umifrac, file = sprintf("%s/%s.matrix.mc_umifrac.csv", out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	write.table(mc_counts,  file = sprintf("%s/%s.matrix.mc_counts.csv",  out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	write.table(mc@mc_fp,   file = sprintf("%s/%s.matrix.mc_fp.csv",      out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)

	# save sc-level tables
	# counts as sparse matrix
	wmat = mat@mat [ , names(mc@mc) ]
	Matrix::writeMM(wmat, file = sprintf("%s/%s.matrix.sc_umi.mtx", out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	R.utils::gzip(sprintf("%s/%s.matrix.sc_umi.mtx", out_fn, run_name), overwrite = TRUE, remove = TRUE)
	write.table(colnames(wmat), file = sprintf("%s/%s.matrix.sc_umi.annot_cells.txt", out_fn, run_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(rownames(wmat), file = sprintf("%s/%s.matrix.sc_umi.annot_genes.txt", out_fn, run_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	
	# cell-level it4 annotations
	sc_annot = data.frame(cell = names(mc@mc), metacell = mc@mc)
	sc_annot = merge(sc_annot, ctt, by = "metacell", reorder = FALSE)
	sc_annot = sc_annot [ , c("cell","metacell","cell_type","color") ]
	sc_annot = sc_annot [ match(names(mc@mc), sc_annot$cell ) , ]
	write.table(sc_annot, file = sprintf("%s/%s.matrix.sc_annot.csv", out_fn, run_name), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
	
	# gc
	gc()
	
	#### Dendrogram and order ####
	
	# obtain confusion matrix
	rec_confu_norm = mc_compute_norm_confu_matrix(
		mc_id = sprintf("%s_it4", run_name),
		graph_id = sprintf("%s_it2", run_name),
		max_deg = 100)

	# dendrogram
	rec_hc = mc_confusion_clustering(rec_confu_norm, clust_method = "ward.D2")
	cut_rec_hc = min(as.dist(rec_hc)) + (max(as.dist(rec_hc)) - min(as.dist(rec_hc))) / 2
	mc_clusts = cutree(rec_hc, h=cut_rec_hc)

	# DENDROGRAM
	# visualize tree and decide on a cut height for color assignment
	cluster_colors_per_mc = ctt$color
	names(cluster_colors_per_mc) = ctt$metacell
	cluster_colors_per_mc = cluster_colors_per_mc [ rec_hc$order ]
	pdf(sprintf("%s/%s.mc_dendrogram.pdf", out_fn, run_name), height = 8, width = length(cluster_colors_per_mc) / 3)
	rec_phy = ape::as.phylo(rec_hc)
	rec_phy$tip.label = rec_hc$order 
	ape::plot.phylo(rec_phy, direction="downwards", las=2, font=1, show.tip.label = FALSE)
	if ( length(unique(cutree(rec_hc, h = cut_rec_hc))) > 1  ) {
		rect.hclust(rec_hc,h=cut_rec_hc, border="darkgray")
	}
	text(x=1:length(cluster_colors_per_mc),y=0, names(cluster_colors_per_mc), col = cluster_colors_per_mc, srt = 0)
	dev.off()
	
	# CONFUSION MATRIX
	# obtain confusion matrix (with reordered metacells)
	pdf(sprintf("%s/%s.mc_confumatrix.pdf", out_fn, run_name), width = 8, height = 8)
	mc_order = 1:nrow(ctt)
	print(plot_complex_heatmap(
		sqrt(rec_confu_norm[ mc_order, mc_order ]),
		name = "sqrt(norm confu)",
		color_mat = c("gray98","#d6e72e","#6fb600","#003f4d"),
		color_min = 0, color_max = quantile(sqrt(rec_confu_norm),0.98),
		cluster_row = FALSE, cluster_col = FALSE,
		colors_row = mc@colors, colors_col = mc@colors,
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		use_raster = TRUE))
	dev.off()
	
}

message("All done!")


