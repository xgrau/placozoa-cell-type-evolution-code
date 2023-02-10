#### Input ####

# libraries
source("../scripts/helper.R")

# input
heatmap_colors = c("gray95","orange","orangered2","#520c52")
icc_fn = "results_alignment_icc/"
ann_fn = "../results_scatlas/results_metacell_it4/"
out_fn = "results_metacell_it4/"
dir.create(out_fn)

# list of comparisons
list_comp = list(
	c("Tadh","TrH2"),
	c("Tadh","Hhon"),
	c("Tadh","HoiH23"),
	c("TrH2","Hhon"),
	c("TrH2","HoiH23"),
	c("Hhon","HoiH23")
)

# list of focus levels
focus_list = list(
	"cts" = "cell_type",
	"bct" = "broad_cell_type"
	"mcs" = "metacell"
)

# list_methods = c("wpearson","wspearman","jaccard","shaferindex","jsd","kld","onls")
list_methods = c("wpearson","jaccard")
list_fcs = c(1.25, 1.5, 2.0)

for (com in list_comp) {
	
	# define input
	sp1 = com[1]
	sp2 = com[2]

	for (i in 1:length(focus_list)) {
		
		focid = names(focus_list)[[i]]
		focus = focus_list[[i]]

		# load cross-species marker data (from `s01`)
		message(sprintf("csps icc | load %s %s-%s csps data...",focus, sp1, sp2))
		csps_fn = sprintf("%s/csps_icc.%s.%s-%s.rds", icc_fn, focid, sp1, sp2)
		csps = readRDS(csps_fn)
		
		# load cell type data
		ctt_sp1_fn = sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp1)
		ctt_sp2_fn = sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp2)
		ctt_sp1 = read.table(ctt_sp1_fn, header = TRUE, comment.char = "", sep = "\t")
		ctt_sp2 = read.table(ctt_sp2_fn, header = TRUE, comment.char = "", sep = "\t")
		
		# get color annotations for each species
		ann_sp1 = unique(ctt_sp1[,c(focus,"color")])
		ann_sp2 = unique(ctt_sp2[,c(focus,"color")])
		ann_sp1 = ann_sp1 [ !duplicated(ann_sp1[,focus]), ]
		ann_sp2 = ann_sp2 [ !duplicated(ann_sp2[,focus]), ]
		ann_sp1 = ann_sp1 [ !grepl("trans",ann_sp1[,focus]), ]
		ann_sp2 = ann_sp2 [ !grepl("trans",ann_sp2[,focus]), ]

			
		# fi loop: fc thresholds to define covariable genes
		for (fi in list_fcs) {

			# get covariable genes
			csps_i = csps
			top_covariable_genes = csps_select_covariable_genes(
				sp1_fp = csps_i$sp1, 
				sp2_fp = csps_i$sp2,
				merged = csps_i$merged,
				cross_fc_thrs = fi,
				cross_n = 1)
			csps_i$top_cross_sp1 = top_covariable_genes[[1]]
			csps_i$top_cross_sp2 = top_covariable_genes[[2]]

			# cor_method = distance methods
			# fi = various FC thresholds to define covariable genes
			for (cor_method in list_methods) {
				
				# matrix for all icc genes
				message(sprintf("csps icc | %s-%s %s matrix (fc>=%.2f)", sp1, sp2, cor_method, fi))
				csps_com_all = csps_correlation_matrix(
					csps = csps_i,
					cor_method = cor_method,
					gene_weights = csps_i$og_pairs_ec_value [ csps_i$og_pairs[,1] %in% csps_i$top_cross_sp1 ],
					use_var_genes = TRUE, 
					fc_thrs = fi,
					add_sps_prefix = TRUE,
					prefix_sp1 = sp1,
					prefix_sp2 = sp2
				)

				# matrix for o2o genes
				message(sprintf("csps icc | %s-%s %s matrix (fc>=%.2f), o2o genes only", sp1, sp2, cor_method, fi))
				markers_o2o = csps_i$og_pairs[,1] [ csps_i$og_pairs_is_o2o & csps_i$og_pairs[,1] %in% csps_i$top_cross_sp1 ]
				csps_com_o2o = csps_correlation_matrix(
					csps = csps_i,
					cor_method = cor_method,
					gene_weights = csps_i$og_pairs_ec_value [ csps_i$og_pairs[,1] %in% markers_o2o ],
					use_var_genes = markers_o2o, 
					fc_thrs = fi,
					add_sps_prefix = TRUE,
					prefix_sp1 = sp1,
					prefix_sp2 = sp2
				)

				# matrix for TFs only
				message(sprintf("csps icc | %s-%s %s matrix (fc>=%.2f), TFs only", sp1, sp2, cor_method, fi))
				# load TF annots
				gene_subset_fn = sprintf("../data/gene_annotations/tfs.%s_genes.curated.csv", sp1)
				gene_subset_v = read.table(gene_subset_fn, sep = "\t", row.names = 1)
				rownames(gene_subset_v) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", sp1), vector_to_fix = rownames(gene_subset_v))
				markers_tfs = csps_i$og_pairs[,1] [ csps_i$og_pairs[,1] %in% rownames(gene_subset_v) & csps_i$og_pairs[,1] %in% csps_i$top_cross_sp1 ]
				csps_com_tfs = csps_correlation_matrix(
					csps = csps_i,
					cor_method = cor_method,
					gene_weights = csps_i$og_pairs_ec_value [ csps_i$og_pairs[,1] %in% markers_tfs ],
					use_var_genes = markers_tfs, 
					fc_thrs = fi,
					add_sps_prefix = TRUE,
					prefix_sp1 = sp1,
					prefix_sp2 = sp2,
					report_overlaps = TRUE
				)

				# write comparison objects
				saveRDS(csps_com_all, sprintf("%s/csps_icc.%s.all.%s-%s.%s.fc%.2f.rds", out_fn, focid, sp1, sp2, cor_method, fi))
				saveRDS(csps_com_o2o, sprintf("%s/csps_icc.%s.o2o.%s-%s.%s.fc%.2f.rds", out_fn, focid, sp1, sp2, cor_method, fi))
				saveRDS(csps_com_tfs, sprintf("%s/csps_icc.%s.tfs.%s-%s.%s.fc%.2f.rds", out_fn, focid, sp1, sp2, cor_method, fi))

				# plots
				list_csps_com = list(
					"all" = csps_com_all, 
					"o2o" = csps_com_o2o, 
					"tfs" = csps_com_tfs
				)
				for (i in 1:length(list_csps_com)) {
				
					# get info for this comparison
					csps_idf = names(list_csps_com)[i]
					csps_com = list_csps_com[[i]]
				
					# heatmap with col/row clustering
					min_val = 0
					max_val = 1
					cor_method_name = cor_method
					if (cor_method %in% c("jsd","kld")) { 
						min_val = quantile(csps_com$cor_matrix, 0.01)
						max_val = quantile(csps_com$cor_matrix, 0.99)
						cor_method_name = sprintf("1-sqrt(%s)", cor_method)
					} else if (cor_method %in% "jsdnp") {
						max_val = quantile(csps_com$cor_matrix, 0.01)
						min_val = quantile(csps_com$cor_matrix, 0.99)
					} else if (cor_method %in% "shaferindex") {
						min_val = 0
						max_val = quantile(csps_com$cor_matrix, 0.95)
					} else if (cor_method %in% "onls") {
						min_val = 0
						max_val = quantile(csps_com$cor_matrix, 0.99)
					}
					
					# exclude trans
					csps_com_com = csps_com$cor_matrix
					csps_com_com = csps_com_com [ !grepl("trans", rownames(csps_com_com)) , !grepl("trans", colnames(csps_com_com)) ]
					
					# plot
					csps_hm = csps_plot_annotated_matrix(
						mat = csps_com_com,
						name = cor_method_name,
						heatmap_colors = heatmap_colors,
						use_raster = FALSE,
						min_val = min_val,
						max_val = max_val,
						row_title = sp1,
						col_title = sprintf("%s (%i markers, FC>=%.2f)", sp2, length(csps_com$var_genes), fi),
						fontsize = 5,
						row_annot = list(ann_sp1),
						col_annot = list(ann_sp2),
						row_annot_legend = FALSE, 
						col_annot_legend = FALSE,
						row_cluster = FALSE,
						col_cluster = FALSE,
						do_dotplot = FALSE)
						
					pdf(sprintf("%s/csps_icc.%s.%s.%s-%s.%s.fc%.2f.pdf", out_fn, focid, csps_idf, sp1, sp2, cor_method, fi), 
						width =  5 + round(ncol(csps$sp2) / 20), 
						height = 4 + round(ncol(csps$sp1) / 20))
					tryCatch(print(csps_hm), error = function(e) { message("Caught an error!") })
					dev.off()
					
				}
				
			}
			
		}
		
	
	}
}

message("All comparisons done!")

