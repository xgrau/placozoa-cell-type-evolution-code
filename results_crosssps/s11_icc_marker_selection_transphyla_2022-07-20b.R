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
	c("Tadh","Spis"),
	c("Tadh","Nvec"),
	c("Tadh","Hvul"),
	c("Tadh","Mmus"),
	c("Tadh","Dmel"),
	c("Tadh","Spolac"),
	c("Tadh","Aque"),
	c("Tadh","Mlei"),
	c("TrH2","Spis"),
	c("TrH2","Nvec"),
	c("TrH2","Hvul"),
	c("TrH2","Mmus"),
	c("TrH2","Dmel"),
	c("TrH2","Spolac"),
	c("TrH2","Aque"),
	c("TrH2","Mlei"),
	c("Hhon","Spis"),
	c("Hhon","Nvec"),
	c("Hhon","Hvul"),
	c("Hhon","Mmus"),
	c("Hhon","Dmel"),
	c("Hhon","Spolac"),
	c("Hhon","Aque"),
	c("Hhon","Mlei"),
	c("HoiH23","Nvec"),
	c("HoiH23","Spis"),
	c("HoiH23","Hvul"),
	c("HoiH23","Mmus"),
	c("HoiH23","Dmel"),
	c("HoiH23","Spolac"),
	c("HoiH23","Aque"),
	c("HoiH23","Mlei")
)


for (com in list_com) { 
	
	# log
	spr = com[1]
	spi = com[2]
	message(sprintf("cross-species markers | %s-%s", spr, spi))

	# orthologs
	og_pairs_fn = sprintf("../data/results_broccoli_ml_extended_gmod/dir_step4/orthologous_pairs.txt")


	# define input
	metacell::scdb_init(mdb_fn, force_reinit = TRUE)
	spr_mcfp = metacell::scdb_mc(sprintf("scdr_%s_it4", spr))@mc_fp
	if (spi %in% c("Nvec","Hvul","Spis","Mmus","Dmel","Spolac","Mlei","Aque")) {
		metacell::scdb_init("../results_scatlas/data/scdb_outgroups", force_reinit = TRUE)
		spi_mcfp = metacell::scdb_mc(spi)@mc_fp
	}

	# load orthology pairs
	if (spi %in% c("Nvec","Hvul","Spis")) {
		og_pairs = clean_og_pairs(og_pairs_fn, sp1 = spr, sp2 = spi, t2g = TRUE, t2g_sp1 = sprintf("../data/reference/%s_long.annot.gtf", spr))
	} else if (spi %in% c("Dmel")) {
		og_pairs = clean_og_pairs(og_pairs_fn, sp1 = spr, sp2 = spi, t2g = TRUE, t2g_sp1 = sprintf("../data/reference/%s_long.annot.gtf", spr), t2g_sp2 = sprintf("../data/reference/%s_long.annot.gtf", spi), t2g_sp2_field = "gene_name")
	} else if (spi %in% c("Spolac")) {
		og_pairs = clean_og_pairs(og_pairs_fn, sp1 = spr, sp2 = spi, t2g = TRUE, t2g_sp1 = sprintf("../data/reference/%s_long.annot.gtf", spr))
		og_pairs[,2] = gsub("_i\\d+_.*","", og_pairs[,2])
	} else {
		og_pairs = clean_og_pairs(og_pairs_fn, sp1 = spr, sp2 = spi, t2g = TRUE, t2g_sp1 = sprintf("../data/reference/%s_long.annot.gtf", spr), t2g_sp2 = sprintf("../data/reference/%s_long.annot.gtf", spi))
	}
	og_pairs = unique(og_pairs)
	
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
