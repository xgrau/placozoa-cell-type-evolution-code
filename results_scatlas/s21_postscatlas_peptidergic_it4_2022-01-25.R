# libraries
library("metacell")
source("../scripts/helper.R")
par(family  = "Arial")

# where to store output?
out_fn = "results_metacell_it4_peptidergic/"

# list of species
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

for (spi in sps_list) {
	
	# init database
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	
	# first load reordered recluster mc solution
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it4_pep",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it4_pep",run_name))
	
	# reorder reclustered metacells according to general clustering
	# get mc annotations from previous iteration
	# first, metacells
	mco = metacell::scdb_mc(sprintf("%s_it4",run_name))
	mc_annot = sca_metacell_annotation_freq(mc@mc, mco@mc [ names(mc@mc) ], freq_min = 0.05)
	# second, cell types
	cto_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
	cto = read.table(cto_fn, header = TRUE, comment.char = "", sep = "\t")
	cto_sc = merge(cto, data.frame(metacell = mco@mc, cell = names(mco@mc)), by = "metacell")
	cto_sc_v = cto_sc$cell_type
	names(cto_sc_v) = cto_sc$cell
	ct_annot = sca_metacell_annotation_freq(mc@mc, cto_sc_v [ names(mc@mc) ], freq_min = 0.05)
	ct_annot_u = sca_metacell_annotation_freq(mc@mc, cto_sc_v [ names(mc@mc) ], collapse = FALSE)
	ct_annot_u$assigned_annot = factor(ct_annot_u$assigned_annot, levels = unique(cto$cell_type))
	# finally, colors
	col_sc_v = cto_sc$color
	names(col_sc_v) = cto_sc$cell
	co_annot_u = sca_metacell_annotation_freq(mc@mc, col_sc_v [ names(mc@mc) ], collapse = FALSE)

	# save
	ctt = data.frame(
		metacell = colnames(mc@mc_fp),
		cell_type = ct_annot_u$assigned_annot,
		color = co_annot_u$assigned_annot,
		prev_mc = mc_annot$assigned_annot,
		prev_mc_freq = mc_annot$annot_freqs,
		prev_ct = ct_annot$assigned_annot,
		prev_ct_freq = ct_annot$annot_freqs
	)
	ctt = ctt [ order(ctt$cell_type), ]
	write.table(ctt, file = sprintf("%s/annotation_mc.%s.it4.peptidergic.tsv", out_fn, spi), sep = "\t", quote = FALSE, row.names = FALSE)

	# reorder metacells according to it4 annotations
	mc_order = ctt$metacell
	mc@mc_fp = mc@mc_fp[,mc_order]
	mc@colors = ctt$color
	names(mc@colors) = ctt$metacell
	
	# save reordered metacells
	metacell::scdb_add_mc(sprintf("%s_it4_pep_ord",run_name), mc = mc)
	mc = metacell::scdb_mc(sprintf("%s_it4_pep_ord",run_name))
	
	# # COUNTS
	# # metacell size and counts
	# mc_stat = scp_plot_mc_size_counts(
	# 	mc_object = mc, 
	# 	mc_color = mc@colors, 
	# 	mat_object = mat, 
	# 	output_file = sprintf("%s/%s.mc_size_counts.pdf", out_fn, run_name),
	# 	return_table = TRUE)
	# write.table(mc_stat, sprintf("%s/%s.mc_size_counts.csv", out_fn, run_name), quote = FALSE, sep = "\t", row.names = FALSE)
	
	# # check batch effects
	# scp_batch_effects_per_mc(
	# 	mc, 
	# 	mat, 
	# 	mc_color = mc@colors, 
	# 	output_file = sprintf("%s/%s.batch_effects_mc.pdf", out_fn, run_name))
	
	# 2D PROJECTION
	# get 2D projection of reclustered subset
	# override_params("data/amphi_metacell_params.yaml","metacell")
	metacell::mcell_mc2d_force_knn(
		mc2d_id = sprintf("%s_it4_pep",run_name),
		mc_id = sprintf("%s_it4_pep_ord",run_name), 
		feats_gset = sprintf("%s_it4_pep",run_name), 
		graph_id = sprintf("%s_it4_pep",run_name))
	mc2d = metacell::scdb_mc2d(sprintf("%s_it4_pep",run_name))
	mc2d_r = sca_reorder_mc2d(mc, mc2d)
	metacell::scdb_add_mc2d(sprintf("%s_it4_pep",run_name), mc2d = mc2d_r)
	mc2d = metacell::scdb_mc2d(sprintf("%s_it4_pep",run_name))

	# plot 2D projection of reclustered subset
	scp_plot_sc_2d(
		mc2d = mc2d,
		mc = mc,
		plot_mcs=TRUE,
		plot_mc_name=TRUE,
		plot_edges=TRUE,
		output_file=sprintf("%s/%s.2dproj.pdf", out_fn, run_name),
		width=8, height=8)

	# plot 2D projection of reclustered subset
	scp_plot_sc_2d(
		mc2d = mc2d,
		mc = mc,
		plot_mcs=TRUE,
		plot_mc_name=TRUE,
		plot_edges=TRUE,
		cex_sc = 0,
		output_file=sprintf("%s/%s.2dproj_mcs.pdf", out_fn, run_name),
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
		height = 12,
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
		height = 12,
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
		heatmap_colors = c("gray98","#91b9ba","#618e8f","#003358"),
		clust_col = mc@colors,
		width = 10,
		height = 12,
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
	mc_counts = sca_mc_gene_counts(mc,mat,0)
	mc_umifrac = sca_mc_gene_umifrac(mc, mc_counts)
	
	# plot tfs and other markers
	dir.create(sprintf("%s/markers/", out_fn))
	for (gene_subset in c("ecm","tfs","sig","sigother","sigpkin","siggpcr","flag","mus","chr","neu","rbp","myo","ion","sterol","ppsyn","carbanhy","neuropept","neugaba")) {
		
		if ( file.exists(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", gene_subset, spi)) ) {
			gene_subset_fn = sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", gene_subset, spi)
			gene_subset_v = read.table(gene_subset_fn, sep = "\t", row.names = 1)
			rownames(gene_subset_v) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = rownames(gene_subset_v))
		} else {
			gene_subset_fn = sprintf("../data/gene_annotations/%s.%s_genes.txt", gene_subset, spi)
			gene_subset_l = unique(read.table(gene_subset_fn, sep = "\t")[,1])
			gene_subset_l = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = gene_subset_l)
			gene_subset_v = data.frame(row.names = gene_subset_l, annot = paste(gene_annot[gene_subset_l,1], stringr::str_trunc(gene_annot[gene_subset_l,2], width = 60), sep = " | "))
			# if statement to avoid having 1-row lists: will add a random gene
			if (nrow(gene_subset_v) == 1)  {
				gene_subset_v = rbind(gene_subset_v, data.frame(row.names = rownames(mc_counts) [ ! rownames(mc_counts) %in% rownames(gene_subset_v) ] [1], annot = "PADDING"))
			}
		}
		
		# wrangle
		gene_subset_m = scp_barplot_heatmap_markers(
			mc_object = mc,
			mat_object = mat,
			mc_counts = mc_counts,
			markers_file = gene_subset_v,
			heatmap_colors =  c("white","orange","orangered2","#520c52"),
			output_file_heatmap = sprintf("%s/markers/%s.markers_%s.heatmap.pdf", out_fn, run_name, gene_subset),
			output_file_barplot = sprintf("%s/markers/%s.markers_%s.barplot.pdf", out_fn, run_name, gene_subset),
			T_totumi = 10,
			width = 16,
			height = NULL,
			use_raster = FALSE,
			min_gene_fc = 1.5,
			min_expression_fc = 1,
			max_expression_fc = 3,
			mc_color = mc@colors,
			print_barplots = TRUE
		)
	}	
	# gc
	gc()
	
	
	# save mc-level tables
	write.table(mc_umifrac, file = sprintf("%s/%s.matrix.mc_umifrac.csv", out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	write.table(mc_counts,  file = sprintf("%s/%s.matrix.mc_counts.csv",  out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	write.table(mc@mc_fp,   file = sprintf("%s/%s.matrix.mc_fp.csv",      out_fn, run_name), sep = "\t", quote = FALSE, row.names = TRUE)
	
	#### Save single cell tables ####

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

}

message("All done!")


