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

# thresholds
cz_large_thr = 15000
cz_small_thr = 100

for (spi in sps_list) {

	#### Load data ####

	# init database
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	metacell::scfigs_init("results_figs/")

	# load umi table
	run_name = sprintf("scdr_%s", spi)
	mat = metacell::scdb_mat(run_name)
	# blacklist one sample
	# mat@mat = mat@mat [ , mat@cell_metadata$batch_set_id != "H13_5_ACMEMeOH_10x_10kc" ]
	mat_cz = Matrix::colSums(mat@mat)
	
	#### Cell filtering ####
	
	# select cells	
	cells_large = names(which(mat_cz > cz_large_thr))
	cells_small = names(which(mat_cz < cz_small_thr))
	
	# are there cells to exclude based on clicktags, other data, etc?
	list_libraries = unique(mat@cell_metadata$batch_set_id)
	list_dulibpath = sprintf("data/dual_samples.%s.class.tsv", list_libraries)
	list_dulibpath = list_dulibpath [ file.exists(list_dulibpath) ]
	cells_blacklisted = c()
	if (length(list_dulibpath) > 0 ) {
		for (cbl_fn in list_dulibpath) {
			cells_dt = read.table(cbl_fn, header = TRUE)
			cells_blacklisted = c(cells_blacklisted, cells_dt$long_cell_name [ cells_dt$classification != spi ])
		}
	}
	
	# manual sample blacklisting
	if (spi == "Hhon") {
		cells_blacklisted = c(cells_blacklisted, rownames(mat@cell_metadata) [ mat@cell_metadata$batch_set_id == "H13_2_ACME_10x_10kc" ])
	}
	
	# keep cells
	cells_keep = names( which(mat_cz >= cz_small_thr & mat_cz <= cz_large_thr) )
	cells_keep = cells_keep [ ! cells_keep %in% cells_blacklisted ]
	
	message(sprintf("%s | n = %i initial cells", spi, ncol(mat@mat)))
	message(sprintf("%s | n = %i small cells", spi, length(cells_small)))
	message(sprintf("%s | n = %i large cells", spi, length(cells_large)))
	message(sprintf("%s | n = %i blacklisted cells", spi, length(cells_blacklisted)))
	message(sprintf("%s | n = %i cells to keep", spi, length(cells_keep)))
	
	# visualise cell size distribution
	pdf(sprintf("results_figs/%s.cs_distribution.pdf", run_name), height = 6, width = 5)
	# plot histogram
	hh = hist(
		log10(mat_cz + 1),
		breaks = 60,
		main = sprintf("UMI/cell %s", run_name),
		col = "gray", 
		border = NA,
		xlab = "log10 UMI+1",
		xlim = c(0,5),
		las = 1)
	hist(log10(mat_cz + 1) [ names(mat_cz) %in% cells_keep ], col = alpha("blue", 0.7), add = TRUE, breaks = 60, border = NA)
	hist(log10(mat_cz + 1) [ names(mat_cz) %in% cells_blacklisted ], col = alpha("magenta3",0.5), add = TRUE, breaks = 60, border = NA)
	legend(
		"topleft", 
		legend = c(
			sprintf("all cells n = %i", length(mat_cz)), 
			sprintf("kept cells n = %i", length(cells_keep)),
			 sprintf("blacklisted cells n = %i", length(cells_blacklisted))), 
		fill = c("gray", "blue", "magenta3"), bty = "n", cex = 0.7, border = NA)
	# title(sub = sprintf("n = %i / %i / %i cells in each interval (out of %i)", length(cells_small), length(cells_keep), length(cells_large), length(mat_cz)), cex.sub = 0.9)
	abline(v = log10(c(cz_small_thr, cz_large_thr)), lty = 2, col = "red")
	text(x = log10(c(cz_small_thr, cz_large_thr)), y = max(max(hh$counts)), c(cz_small_thr, cz_large_thr), col = "darkred")
	
	# plot histogram per sample
	for (sai in unique(mat@cell_metadata$batch_set_id)) {
		bcs_sai = rownames(mat@cell_metadata) [ mat@cell_metadata$batch_set_id == sai ]
		mat_sai = mat@mat[, bcs_sai ]
		mat_sai_cz = Matrix::colSums(mat_sai)
		hh = hist(
			log10(mat_sai_cz + 1),
			breaks = 60,
			main = sprintf("UMI/cell %s | batch: %s", run_name, sai),
			col = "gray", 
			border = NA,
			xlab = "log10 UMI+1",
			xlim = c(0,5),
			las = 1)
		abline(v = log10(c(cz_small_thr, cz_large_thr)), lty = 2, col = "red")
		
	}
	
	# plot full distribution
	plot(
		sort(log10(mat_cz + 1)), 
		main = sprintf("UMI/cell %s", run_name),
		ylim = c(0,5),
		col = "blue", ylab = "log10 UMI+1", cex = 0.7)
	title(sub = sprintf("n = %i / %i / %i cells in each interval (out of %i)", length(cells_small), length(cells_keep), length(cells_large), length(mat_cz)), cex.sub = 0.9)
	abline(h = log10(c(cz_small_thr, cz_large_thr)), lty = 2, col = "red")
	text(y = log10(c(cz_small_thr, cz_large_thr)), x = 0, c(cz_small_thr, cz_large_thr), col = "darkred", pos = 4)
	dev.off()
	
	# remove large and small cells
	# also remove blacklisted cells
	metacell::mcell_mat_ignore_cells(
		new_mat_id = sprintf("%s_it0", run_name),
		mat_id = run_name,
		ig_cells=c(cells_small, cells_large, cells_blacklisted))
	# reload matrix
	mat = scdb_mat(sprintf("%s_it0", run_name))
	
	#### Metacell solution: initial run ####
	
	# gene stats
	metacell::mcell_add_gene_stat(
		gstat_id = sprintf("%s_it0", run_name),
		mat_id = sprintf("%s_it0", run_name),
		force = TRUE)
	
	# find feature gene set for this subset of data
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
	
	# perform reclustering (to get a cgraph object)
	metacell::mcell_add_cgraph_from_mat_bknn(
		mat_id = sprintf("%s_it0", run_name),
		gset_id = sprintf("%s_it0", run_name),
		graph_id = sprintf("%s_it0", run_name),
		K = 100,
		dsamp = FALSE)
		
	# coclustering graph via resampling (to get a coclust object)
	metacell::mcell_coclust_from_graph_resamp(
		coc_id = sprintf("%s_it0", run_name),
		graph_id = sprintf("%s_it0", run_name),
		min_mc_size = 10,
		p_resamp = 0.75,
		n_resamp = 1000)
		
	# balanced coclustering (to get a mc object)
	metacell::mcell_mc_from_coclust_balanced(
		mc_id = sprintf("%s_it0", run_name),
		coc_id = sprintf("%s_it0", run_name),
		mat_id = sprintf("%s_it0", run_name),
		K = 30,
		min_mc_size = 12,
		alpha = 2)
		
	
	#### Metacell solution: no batch genes ####
	
	# reload metacells and cell data
	mc = metacell::scdb_mc(id = sprintf("%s_it0", run_name))
	mat = metacell::scdb_mat(id = sprintf("%s_it0", run_name))
	
	# find batch-correlated genes at the metacell level
	pdf(sprintf("results_figs/%s.batch_genes.pdf", run_name), height = 9, width = 10)
	batch_data_mc = sca_batch_correlated_genes_mc(mc_object = mc, mat_object = mat, cor_thr = 0.5, method = "spearman", do_plots = TRUE)
	write.table(batch_data_mc$genes, sprintf("results_figs/%s.batch_genes.txt", run_name), row.names = FALSE, col.names = FALSE, quote = FALSE)
	dev.off()
	
	if (length(batch_data_mc$genes) > 0) {
		
		# define gene set
		metacell::mcell_gset_iter_multi(
			gstat_id = sprintf("%s_it0", run_name),
			gset_id = sprintf("%s_it0", run_name),
			T_tot = 100,
			T_top3 = 2,
			T_szcor = -0.05,
			T_niche = 0.01,
			force_new = TRUE,
			blacklist = union(genes_blacklisted, batch_data_mc$genes))
		
		# some diagnostic plots
		pdf(sprintf("results_figs/%s.cs_diagnostics.pdf", run_name), height = 4.5, width = 5)
		mcell_plot_gstats_mod(
			gstat_id = sprintf("%s_it0", run_name),
			gset_id = sprintf("%s_it0", run_name))
		dev.off()
		
		# perform reclustering (to get a cgraph object)
		metacell::mcell_add_cgraph_from_mat_bknn(
			mat_id = sprintf("%s_it0", run_name),
			gset_id = sprintf("%s_it0", run_name),
			graph_id = sprintf("%s_it0", run_name),
			K = 100,
			dsamp = FALSE)
		
		# coclustering graph via resampling (to get a coclust object)
		metacell::mcell_coclust_from_graph_resamp(
			coc_id = sprintf("%s_it0", run_name),
			graph_id = sprintf("%s_it0", run_name),
			min_mc_size = 10,
			p_resamp = 0.75,
			n_resamp = 1000)
		
		# balanced coclustering (to get a mc object)
		metacell::mcell_mc_from_coclust_balanced(
			mc_id = sprintf("%s_it0", run_name),
			coc_id = sprintf("%s_it0", run_name),
			mat_id = sprintf("%s_it0", run_name),
			K = 30,
			min_mc_size = 12,
			alpha = 2)
		
	}
	
	#### Dendrogram and order ####
	
	mc = metacell::scdb_mc(id = sprintf("%s_it0", run_name))
	mat = metacell::scdb_mat(id = sprintf("%s_it0", run_name))
	
	# obtain confusion matrix
	rec_confu_norm = mc_compute_norm_confu_matrix(
		mc_id = sprintf("%s_it0", run_name), 
		graph_id = sprintf("%s_it0", run_name), 
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
	metacell::scdb_add_mc(id = sprintf("%s_it0",run_name), mc_reord)
	
	# DENDROGRAM
	# visualize tree and decide on a cut height for color assignment
	pdf(sprintf("results_figs/%s.mc_dendrogram.pdf", run_name), height = 8, width = length(cluster_colors_per_mc) / 4)
	rec_phy = ape::as.phylo(rec_hc)
	ape::plot.phylo(rec_phy, direction="downwards", las=2, font=1, tip.color = "gray")
	if ( length(unique(cutree(rec_hc, h = cut_rec_hc))) > 1  ) {
		rect.hclust(rec_hc,h=cut_rec_hc, border="darkgray")
	}
	text(x=1:length(cluster_colors_per_mc),y=0, names(cluster_colors_per_mc), col = cluster_colors_per_mc, srt = 270)
	dev.off()
	
	# plot confusion matrix
	# obtain confusion matrix (with reordered metacells)
	mc = metacell::scdb_mc(id = sprintf("%s_it0", run_name))
	rec_confu_norm = mc_compute_norm_confu_matrix(
		mc_id = sprintf("%s_it0", run_name), 
		graph_id = sprintf("%s_it0", run_name), 
		max_deg = 100)
	
	pdf(sprintf("results_figs/%s.mc_confumatrix.pdf", run_name), width = 24, height = 24)
	print(plot_complex_heatmap(
		sqrt(rec_confu_norm),
		name = "sqrt(norm confu)",
		color_min = 0, color_max = max(sqrt(rec_confu_norm)), 
		cluster_row = FALSE, cluster_col = FALSE,
		colors_row = mc@colors, colors_col = mc@colors,
		use_raster = FALSE))
	dev.off()
	
	
}

message("All done!")
