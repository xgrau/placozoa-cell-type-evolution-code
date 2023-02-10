# libraries
library("metacell")
source("../scripts/geneSetAnalysis.R")
source("../scripts/helper.R")
par(family  = "Arial")

# parameters
tgconfig::set_param("mc_cores", 1, "metacell") # this allows gstat to work without crashing due to excessive memory usage (and it's equally fast)

# data index
sps_list = c("Tadh","TrH2","Hhon","HoiH23")
gene_blacklist = list(
	"Tadh"   = c("../data/reference/Tadh_long.pep.gene_lists_ribosomal.csv",   "../data/reference/Tadh_long.pep.gene_lists_histones.csv"),
	"TrH2"   = c("../data/reference/TrH2_long.pep.gene_lists_ribosomal.csv",   "../data/reference/TrH2_long.pep.gene_lists_histones.csv"),
	"HoiH23" = c("../data/reference/HoiH23_long.pep.gene_lists_ribosomal.csv", "../data/reference/HoiH23_long.pep.gene_lists_histones.csv"),
	"Hhon"   = c("../data/reference/Hhon_long.pep.gene_lists_ribosomal.csv",   "../data/reference/Hhon_long.pep.gene_lists_histones.csv")
)

# output
out_fn = "results_metacell_it1/"

for (spi in sps_list) {


	#### Load data ####

	# init database
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	metacell::scfigs_init("results_figs/")

	# load umi table
	run_name = sprintf("scdr_%s", spi)
	mat = metacell::scdb_mat(sprintf("%s_it1",run_name))
	mat_cz = Matrix::colSums(mat@mat)
	
	
	#### Cell filtering ####
	
	# the `it1` matrix is already filtered after it0; all removed cells are lised in mat@ignore_cells
	# blacklisting criteria:
	# - small/large cells
	# - blacklisted based on clicktags
	# - blacklisted based on sick mcs
	
	# additional filtering criteria: excess of clicktags
	
	
	#### Gene set ####

	# gene stats
	metacell::mcell_add_gene_stat(
		gstat_id = sprintf("%s_it1", run_name),
		mat_id = sprintf("%s_it1", run_name),
		force = TRUE)
		
	# first, check if there are any genes to blacklist
	genes_blacklisted = c()
	if (!is.null(gene_blacklist[[spi]])) {
		for (bli in gene_blacklist[[spi]]) {
			genes_blacklisted = c(genes_blacklisted, read.table(bli, header = FALSE, sep = "\t")[,2])
		}
		genes_blacklisted = unique(genes_blacklisted)
		message(sprintf("%s | genes to blacklist = %i (from %i files)", spi, length(genes_blacklisted), length(gene_blacklist[[spi]]) ) )
	} else {
		message(sprintf("%s | genes to blacklist = %i (no file provided)", spi, length(genes_blacklisted)))
	}
	genes_blacklisted = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = genes_blacklisted)
	
	# further load batch-correlated genes identified in it0
	genes_blacklisted = union(genes_blacklisted, read.table(sprintf("results_figs/%s.batch_genes.txt", run_name))[,1])
	
	# define gene set
	metacell::mcell_gset_filter_multi(
		gstat_id = sprintf("%s_it1", run_name),
		gset_id = sprintf("%s_it1", run_name),
		T_tot = 100,
		T_top3 = 2,
		T_szcor = -0.05,
		T_niche = 0.01,
		force_new = TRUE,
		blacklist = genes_blacklisted)
	

	#### Metacells from coclustering ####

	# perform reclustering (to get a cgraph object)
	metacell::mcell_add_cgraph_from_mat_bknn(
		mat_id = sprintf("%s_it1", run_name),
		gset_id = sprintf("%s_it1", run_name),
		graph_id = sprintf("%s_it1", run_name),
		K = 100,
		dsamp = FALSE)
		
	# coclustering graph via resampling (to get a coclust object)
	metacell::mcell_coclust_from_graph_resamp(
		coc_id = sprintf("%s_it1", run_name),
		graph_id = sprintf("%s_it1", run_name),
		min_mc_size = 10,
		p_resamp = 0.75,
		n_resamp = 1000)
		
	# balanced coclustering (to get a mc object)
	metacell::mcell_mc_from_coclust_balanced(
		mc_id = sprintf("%s_it1", run_name),
		coc_id = sprintf("%s_it1", run_name),
		mat_id = sprintf("%s_it1", run_name),
		K = 30,
		min_mc_size = 12,
		alpha = 2)
		
	# # split metacells?
	# mcell_mc_split_filt_mod(
	# 	new_mc_id = sprintf("%s_it1_split", run_name),
	# 	mc_id = sprintf("%s_it1", run_name), 
	# 	mat_id = sprintf("%s_it1", run_name),
	# 	T_lfc = 8
	# )
		
	
	#### Dendrogram and order ####
	
	mc = metacell::scdb_mc(id = sprintf("%s_it1", run_name))
	mat = metacell::scdb_mat(id = sprintf("%s_it1", run_name))
	
	# obtain confusion matrix
	rec_confu_norm = mc_compute_norm_confu_matrix(
		mc_id = sprintf("%s_it1", run_name), 
		graph_id = sprintf("%s_it1", run_name), 
		max_deg = 100)

	# dendrogram
	rec_hc = mc_confusion_clustering(rec_confu_norm, clust_method = "ward.D2")
	cut_rec_hc = min(as.dist(rec_hc)) + (max(as.dist(rec_hc)) - min(as.dist(rec_hc))) / 2
	# cut_rec_hc = quantile(as.dist(rec_hc),0.25)
	mc_clusts = cutree(rec_hc, h=cut_rec_hc)

	# assign colors to each cluster of metacells
	color_palette = colorRampPalette(c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen","deepskyblue","cadetblue1","mediumblue","darkviolet","violet"))
	cluster_colors = color_palette(n=length(table(mc_clusts[rec_hc$order])))
	cluster_colors_per_mc = rep(cluster_colors, table(factor(mc_clusts[rec_hc$order],levels=unique(mc_clusts[rec_hc$order]))))
	names(cluster_colors_per_mc) = 1:length(cluster_colors_per_mc)
	
	# reorder metacells based on hierarchical clustering
	mc_reord = metacell::mc_reorder(mc,rec_hc$order)
	mc_reord@colors=cluster_colors_per_mc
	metacell::scdb_add_mc(id = sprintf("%s_it1",run_name), mc_reord)
	
	# DENDROGRAM
	# visualize tree and decide on a cut height for color assignment
	pdf(sprintf("%s/%s.mc_dendrogram.pdf", out_fn, run_name), height = 8, width = length(cluster_colors_per_mc) / 4)
	rec_phy = ape::as.phylo(rec_hc)
	ape::plot.phylo(rec_phy, direction="downwards", las=2, font=1, tip.color = "gray")
	if ( length(unique(cutree(rec_hc, h = cut_rec_hc))) > 1  ) {
		rect.hclust(rec_hc,h=cut_rec_hc, border="darkgray")
	}
	text(x=1:length(cluster_colors_per_mc),y=0, names(cluster_colors_per_mc), col = cluster_colors_per_mc, srt = 270)
	dev.off()
	
	# CONFUSION MATRIX
	# obtain confusion matrix (with reordered metacells)
	mc = metacell::scdb_mc(id = sprintf("%s_it1", run_name))
	rec_confu_norm = mc_compute_norm_confu_matrix(
		mc_id = sprintf("%s_it1", run_name), 
		graph_id = sprintf("%s_it1", run_name), 
		max_deg = 100)
	
	pdf(sprintf("%s/%s.mc_confumatrix.pdf", out_fn, run_name), width = 24, height = 24)
	print(plot_complex_heatmap(
		sqrt(rec_confu_norm),
		name = "sqrt(norm confu)",
		color_min = 0, color_max = max(sqrt(rec_confu_norm)), 
		cluster_row = FALSE, cluster_col = FALSE,
		colors_row = mc@colors, colors_col = mc@colors,
		use_raster = FALSE))
	dev.off()
	
	# clean memory
	gc()
	
}

message("All done!")
