# libraries
library("metacell")
source("../scripts/helper.R")
par(family  = "Arial")

# initial metacell iteration stored here
out_fn = "results_metacell_it0/"

# species list
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

for (spi in sps_list) {
	
	# first load reordered recluster mc solution
	metacell::scdb_init("data/scdb/", force_reinit=TRUE)
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it0",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it0",run_name))
	mc_order = 1:ncol(mc@mc_fp)
	
	# COUNTS
	# metacell size and counts
	mc_stat = scp_plot_mc_size_counts(
		mc_object = mc, 
		mc_color = mc@colors, 
		mat_object = mat, 
		output_file = sprintf("%s/%s.mc_size_counts.pdf", out_fn, run_name),
		return_table = TRUE)
	write.table(mc_stat, sprintf("%s/%s.mc_size_counts.csv", out_fn, run_name), quote = FALSE, sep = "\t", row.names = FALSE)
	
	
	# 2D PROJECTION
	# get 2D projection of reclustered subset
	# override_params("data/amphi_metacell_params.yaml","metacell")
	metacell::mcell_mc2d_force_knn(
		mc2d_id = sprintf("%s_it0",run_name),
		mc_id = sprintf("%s_it0",run_name), 
		feats_gset = sprintf("%s_it0",run_name), 
		graph_id = sprintf("%s_it0",run_name))
	mc2d = metacell::scdb_mc2d(sprintf("%s_it0",run_name))
	
	# plot 2D projection of reclustered subset
	scp_plot_sc_2d(
		mc2d = mc2d,
		mc = mc,
		plot_mcs=TRUE,
		plot_mc_name=TRUE,
		plot_edges=TRUE,
		output_file=sprintf("%s/%s.2dproj.pdf", out_fn, run_name),
		width=8, height=8)
	

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
		per_clust_genes = 30,
		gene_min_fold = 1.5,
		gene_font_size = 5,
		clust_ord = mc_order)
	
	# draw metacell heatmap
	scp_plot_cmod_markers_mc(
		marker_data_list = marker_data_list,
		output_file = sprintf("%s/%s.exp_global_mc.pdf", out_fn, run_name),
		heatmap_colors = c("white","orange","orangered2","#520c52"),
		clust_col = mc@colors,
		width = max(8, round(length(marker_data_list$clust_ord) / 10)),
		height = round(length(marker_data_list$genes) / 15) + 1,
		show_gene_names = TRUE,
		highlight_genes = names(dic_genetf_gene_name),
		gene_chr_limit = 60,
		show_clust_borders = TRUE,
		use_raster=TRUE,
		max_expression_fc = 3,
		min_expression_fc = 1
	)
	
	# compressed metacell map
	scp_plot_cmod_markers_mc(
		marker_data_list = marker_data_list,
		output_file = sprintf("%s/%s.exp_global_mc-c.pdf", out_fn, run_name),
		heatmap_colors = c("white","orange","orangered2","#520c52"),
		clust_col = mc@colors,
		width = 4,
		height = 4,
		show_gene_names = FALSE,
		highlight_genes = names(dic_genetf_gene_name),
		gene_chr_limit = 60,
		show_clust_borders = TRUE,
		use_raster = TRUE,
		max_expression_fc = 3,
		min_expression_fc = 1
	)
	
	# draw single cell heatmap
	scp_plot_cmod_markers_sc(
		marker_data_list = marker_data_list,
		mc_object = mc,
		mat_object = mat,
		output_file = sprintf("%s/%s.exp_global_sc.pdf", out_fn, run_name),
		heatmap_colors = c("white","orange","orangered2","#520c52"),
		clust_col = mc@colors,
		width = max(12, round(length(marker_data_list$clust_ord) / 2)),
		height = round(length(marker_data_list$genes) / 15) + 1,
		show_gene_names = TRUE,
		highlight_genes = names(dic_genetf_gene_name),
		smoothen = 5,
		gene_chr_limit = 60,
		show_clust_borders=TRUE,
		use_raster=TRUE,
		max_expression_fc = 3,
		min_expression_fc = 1
	)
	
	# get mc counts, umifrac
	mc_counts = sca_mc_gene_counts(mc,mat,0)
	mc_umifrac = sca_mc_gene_umifrac(mc, mc_counts)
	
	# COUNTS
	# metacell size and counts
	mc = metacell::scdb_mc(id = sprintf("%s_it0", run_name))
	scp_plot_mc_size_counts(
		mc_object = mc, 
		mc_color = mc@colors, 
		mat_object = mat, 
		output_file = sprintf("%s/%s.mc_size_counts.pdf", out_fn, run_name))
	
	# check batch effects
	scp_batch_effects_per_mc(
		mc, 
		mat, 
		mc_color = mc@colors, 
		output_file = sprintf("%s/%s.batch_effects_mc.pdf",out_fn, run_name))

	# MARKERS	
	# plot tfs
	gene_subset = "tfs"
	gene_subset_v = read.table(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", gene_subset, spi), sep = "\t", row.names = 1)
	rownames(gene_subset_v) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = rownames(gene_subset_v))
	gene_subset_m = scp_barplot_heatmap_markers(
		mc_object = mc,
		mat_object = mat,
		mc_counts = mc_counts,
		markers_file = gene_subset_v,
		heatmap_colors =  c("white","orange","orangered2","#520c52"),
		output_file_heatmap = sprintf("%s/%s.markers_%s.heatmap.pdf", out_fn, run_name, gene_subset),
		output_file_barplot = sprintf("%s/%s.markers_%s.barplot.pdf", out_fn, run_name, gene_subset),
		T_totumi = 10,
		width = 12,
		height = NULL,
		use_raster = TRUE,
		min_gene_fc = 1.5,
		min_expression_fc = 1,
		max_expression_fc = 3,
		mc_color = mc@colors,
		print_barplots = TRUE
	)
	
	# cluster names
	mc_clusters = as.numeric(factor(mc@colors, levels = unique(mc@colors)))
	mc_clusters_colname = unlist(sapply(mc@colors, function(c) plotrix::color.id(c)[1]))
	mc_clusters = paste(mc_clusters, mc_clusters_colname, sep = "_")
	
	
	
	#### Putative annotations from previous dataset ####
	
	# load original metacells and annotations
	metacell::scdb_init("data/scdb_nee18/",force_reinit=TRUE)
	mo  = metacell::scdb_mc("mc_reord")
	cto = read.table("data/scdb_nee18/cell_type_table.txt", header = TRUE, sep = "\t") 
	rownames(cto) = cto$metacell
	# wrangle names
	rownames(mo@mc_fp) = sub("_P","_TriadG", rownames(mo@mc_fp))
	
	# load og pairs to match with original dataset
	if (spi != "Tadh") {
		og_pairs = clean_og_pairs(
			og_pairs_fn = "../data/results_broccoli_ml/dir_step4/orthologous_pairs.txt", 
			sp1 = "Tadh", 
			sp2 = spi, 
			t2g = TRUE, 
			t2g_sp1 = "../data/reference/Tadh_long.annot.gtf",
			t2g_sp2 = sprintf("../data/reference/%s_long.annot.gtf", spi)
		)
		list_o2o_sp1 = names(which(table(og_pairs[,1]) == 1))
		list_o2o_sp2 = names(which(table(og_pairs[,2]) == 1))
		bool_o2o = og_pairs[,1] %in% list_o2o_sp1 & og_pairs[,2] %in% list_o2o_sp2
		og_pairs_o2o = og_pairs [ bool_o2o, ]
		og_pairs_o2o_f = og_pairs_o2o [ og_pairs_o2o[,1] %in% rownames(mo@mc_fp) & og_pairs_o2o[,2] %in% rownames(mc@mc_fp),  ]
	} else {
		og_pairs_o2o = data.frame(sp1 = rownames(mc@mc_fp), sp2 = rownames(mc@mc_fp))
		og_pairs_o2o_f = og_pairs_o2o [ og_pairs_o2o[,1] %in% rownames(mo@mc_fp) & og_pairs_o2o[,2] %in% rownames(mc@mc_fp),  ]
	}
	
	# order genes in original matrix
	mo_fp = mo@mc_fp [ og_pairs_o2o_f[,1], ]
	mc_fp = mc@mc_fp [ og_pairs_o2o_f[,2], ]
	
	cor_mo_mc = cor(mc_fp, mo_fp, method = "spearman")
	top_mo_mc = apply(cor_mo_mc, 1, function(r) {
		best_mcs_l = names(sort(r [ r > 0 ], decreasing = TRUE)[1:5])
		best_mcs_c = sort(r [ r > 0 ], decreasing = TRUE)[1:5]
		best_cto_l = cto [as.numeric(best_mcs_l), "broad_cell_type"]
		best_mcs_cs = aggregate(best_mcs_c, list(best_cto_l), sum)
		best_cto_t = best_mcs_cs$Group.1 [ order(best_mcs_cs$x, decreasing = TRUE)[1] ]
		return(list(annot = best_cto_t, mcs = best_cto_l, mc_cors = best_mcs_c))
	})
	top_mo_mc_annot = unlist(lapply(top_mo_mc, function(i) i$annot))
	
	
	
	#### Save metacell table ####
	# iteration 0, treat transferred annotations with care!
	mc_t = data.frame(
		metacell = colnames(mc@mc_fp),
		cell_type = mc_clusters,
		color = mc@colors,
		transferred_annotation = top_mo_mc_annot
	)
	write.table(mc_t, file = sprintf("%s/annotation_mc.%s.it0.tsv", out_fn, spi), sep = "\t", quote = FALSE, row.names = FALSE)
	
	# save mc-level tables
	write.table(mc_umifrac, file = sprintf("%s/%s.matrix.mc_umifrac.csv", out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	write.table(mc_counts,  file = sprintf("%s/%s.matrix.mc_counts.csv",  out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	write.table(mc@mc_fp,   file = sprintf("%s/%s.matrix.mc_fp.csv",      out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	
}

message("All done!")



