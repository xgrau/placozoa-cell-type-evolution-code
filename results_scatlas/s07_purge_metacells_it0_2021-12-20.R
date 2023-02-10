# libraries
library("metacell")
source("../scripts/Downstream_functions.R")
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
	mc_counts = sca_mc_gene_counts(mc, mat, 0)

	# load cell type annotations
	ctt_fn = sprintf("data/annotation_mc.%s.it0.tsv", spi)
	ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")

	# drop due to low umis
	mcs_to_drop_by_umi = sca_drop_mcs_by_umis(
		mc_object = mc, 
		mc_color = mc@colors,
		mat_object = mat, 
		mc_clusters = ctt$cell_type,
		do_global_filter = TRUE,
		mc_counts = mc_counts,
		min_median_umis = 0,
		min_total_umis = 0,
		output_file = sprintf("%s/%s.drop_mcs_by_umi.pdf", out_fn, run_name),
		do_plots = TRUE)
		
	# drop due to low fcs in general
	mcs_to_drop_by_fc  = sca_drop_mcs_by_fc(mc, min_fc = 1.5, min_markers = 10)
	
	# drop due to low fcs amongst tfs
	list_tfs = read.table(sprintf("../data/gene_annotations/tfs.%s_genes.txt", spi), sep = "\t")[,1]
	list_tfs = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = list_tfs)
	list_tfs = list_tfs [ list_tfs %in% rownames(mc@mc_fp) ]
	mcs_to_drop_by_tf = sca_drop_mcs_by_fc(mc, min_fc = 1.5, min_markers = 1, key_markers = list_tfs)

	# get metacells and cells to drop
	mcs_to_drop = as.character(sort(as.numeric(unique(c(mcs_to_drop_by_umi, mcs_to_drop_by_fc, mcs_to_drop_by_tf)))))
	cells_to_drop = names(mc@mc) [ mc@mc %in% as.numeric(mcs_to_drop) ]
	cells_to_drop_mc = mc@mc [ mc@mc %in% as.numeric(mcs_to_drop) ]
	
	# save table with dropped cells and reason
	drt = data.frame(
		cell = cells_to_drop,
		mc   = cells_to_drop_mc,
		fail_umi_filter = cells_to_drop_mc %in% mcs_to_drop_by_umi,
		fail_fc_filter = cells_to_drop_mc %in% mcs_to_drop_by_fc,
		fail_tf_filter = cells_to_drop_mc %in% mcs_to_drop_by_tf
	)
	write.table(drt, sprintf("%s/%s.drop_mcs_cell_list.csv", out_fn, run_name), row.names = FALSE, sep = "\t", quote = FALSE)
	
	# save new mat object for it1
	message(sprintf("save mat.%s_it1...", run_name))
	mcell_mat_ignore_cells(
		new_mat_id = sprintf("%s_it1",run_name), 
		mat_id = sprintf("%s_it0",run_name),
		ig_cells = union(cells_to_drop, mat@ignore_cells))
		
	# # log
	message(sprintf("%s drop: %i cells", run_name, length(cells_to_drop)))
	message(sprintf("%s drop: %i metacells: %s", run_name, length(mcs_to_drop), paste(mcs_to_drop, collapse = ",")))
	mat0 = scdb_mat(sprintf("%s_it0", run_name))
	mat1 = scdb_mat(sprintf("%s_it1", run_name))
	message(sprintf("%s it0 size: %i cells", run_name, ncol(mat0@mat)))
	message(sprintf("%s it1 size: %i cells", run_name, ncol(mat1@mat)))

}

message("all done!")