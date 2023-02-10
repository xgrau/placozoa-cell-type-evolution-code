#### Input ####

# libraries
library("metacell")
library("tgconfig")
source("../scripts/helper.R")
library("igraph")
library("WGCNA")

# paths to input data
mdb_fn = "../results_scatlas/data/scdb/"
out_fn = "results_alignment_icc/"
dir.create(out_fn, showWarnings = FALSE)

#### ICC markers per species pair ####

# Find gene markers for all species pairs
# run icc for each query species

list_com = list(
	c("Tadh","TrH2"),
	c("Tadh","Hhon"),
	c("Tadh","HoiH23"),
	c("TrH2","Hhon"),
	c("TrH2","HoiH23"),
	c("Hhon","HoiH23")
)

# init
metacell::scdb_init(mdb_fn, force_reinit = TRUE)

for (com in list_com) { 
	
	# log
	spr = com[1]
	spi = com[2]
	message(sprintf("cross-species markers | %s-%s", spr, spi))

	# define input
	spr_mcfp = metacell::scdb_mc(sprintf("scdr_%s_it4", spr))@mc_fp
	spi_mcfp = metacell::scdb_mc(sprintf("scdr_%s_it4", spi))@mc_fp
	og_pairs_fn = sprintf("../results_annotation/results_broccoli_ml/dir_step4/orthologous_pairs.txt")

	# load orthology pairs
	og_pairs = clean_og_pairs(og_pairs_fn, sp1 = spr, sp2 = spi, t2g = TRUE, t2g_sp1 = sprintf("../data/reference/%s_long.annot.gtf", spr), sprintf("../data/reference/%s_long.annot.gtf", spi))
	
	# get expression conservation values from ICC
	csps_genes_icc_fp = csps_markers_icc(
		mat_sp1 = spr_mcfp, 
		mat_sp2 = spi_mcfp, 
		og_pairs = og_pairs,
		niter = 100,
		icc_thr = 0.05,
		method = "pearson",
		do_duplicates = TRUE,
		num_o2o_ref_genes = 1000,
		use_variable_o2o_genes = FALSE,
		nthreads_icc = 10,
		nthreads_dup = 2)
	
	# save cross-species object
	saveRDS(
		csps_genes_icc_fp$csps,
		file = sprintf("%s/csps_icc.mcs.%s-%s.rds", out_fn, spr, spi))
	
	# save markers
	write.table(
		csps_genes_icc_fp$ec_markers,
		file = sprintf("%s/csps_icc.%s-%s.markers_ec.csv", out_fn, spr, spi),
		quote = FALSE, sep = "\t", row.names = FALSE)
		
	# save ec scores for all duplicate pairs
	write.table(
		csps_genes_icc_fp$ec_duplicates,
		file = sprintf("%s/csps_icc.%s-%s.paralogs_ec.csv", out_fn, spr, spi),
		quote = FALSE, sep = "\t", row.names = FALSE)
		
	gc()
		
}
