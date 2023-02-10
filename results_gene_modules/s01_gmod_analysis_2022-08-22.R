#### Input ####

# libraries
library("metacell")
library("scales")
source("../scripts/Gene_module_functions.R")
source("../scripts/Downstream_functions.R")
source("../scripts/helper.R")

## Placozoans ##

# paths to input data
# spslist
sps_list = c("Tadh","TrH2","Hhon","HoiH23")
sps_list = sps_list[4]
set_list = list(
	global =      list(inp_fn = "../results_scatlas/results_metacell_it4/",             out_fn = "results_gmod_it4c/",             mc_sprintf_string = "%s_it4",         mat_sprintf_string = "%s_it2",     ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv"),
	peptidergic = list(inp_fn = "../results_scatlas/results_metacell_it4_peptidergic/", out_fn = "results_gmod_it4c_peptidergic/", mc_sprintf_string = "%s_it4_pep_ord", mat_sprintf_string = "%s_it4_pep", ctt_sprtinf_string = "%s/annotation_mc.%s.it4.peptidergic.tsv")
)

# loop
for (set in set_list) {

	out_fn = set$out_fn
	inp_fn = set$inp_fn
	mc_sprintf_string = set$mc_sprintf_string
	mat_sprintf_string = set$mat_sprintf_string
	ctt_sprtinf_string = set$ctt_sprtinf_string

	dir.create(out_fn, showWarnings = FALSE)

	for (spi in sps_list) {
		
		# info
		run_name = sprintf("scdr_%s", spi)
		
		# load data
		message(sprintf("%s | load %s", out_fn, spi))
		metacell::scdb_init("../results_scatlas/data/scdb/",force_reinit=TRUE)
		mc  = metacell::scdb_mc( sprintf(mc_sprintf_string, run_name))
		# mat = metacell::scdb_mat(sprintf(mat_sprintf_string,run_name))
		
		# load cell type annotations for it4
		ctt_fn = sprintf(ctt_sprtinf_string, inp_fn, spi)
		ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")
		
		# # counts
		# ct = sca_cell_type_fp(ctt, mc, mat)
		# ct_counts  = sca_mc_gene_counts(mc_object = ct, mat_object = mat, T_totumi = 5)
		# ct_umifrac = sca_mc_gene_umifrac(mc_object = ct, mc_counts = ct_counts)
		
		# define variable genes
		fc_thr = 1.25
		var_genes = names(which(apply(mc@mc_fp, 1, max) > fc_thr))
		# var_genes = var_genes [ !grepl("orphan", var_genes) ]
		write.table(var_genes, sprintf("%s/gmod_%s.wgcna_object.var_genes.txt", out_fn, spi), quote = FALSE, row.names = FALSE, col.names = FALSE)
		
		# determine soft power
		gmod_determineSoftPowerWGCNA(
			data = mc@mc_fp[var_genes,], 
			output_file = sprintf("%s/gmod_%s.softpower.pdf", out_fn, spi)
		)
		
		# run wgcna with soft power value determined above (first value above line)
		mc_wgcna = gmod_runWGCNA(
			data = mc@mc_fp[var_genes,],
			propGenes = 1,
			softPower = 7,
			cor_method = "pearson",
			signedNetwork = TRUE)
		
		# plot dendrogram
		gmod_plotModulesCut(
			mc_wgcna, 
			output_file = sprintf("%s/gmod_%s.dendrogram.pdf", out_fn, spi)
		)
		
		# calculate module eigengenes
		mc_wgcna_me = gmod_calculateModuleEigengenes(
			mc_wgcna,
			split = 4, 
			minClusterSize = 10, 
			cutHeight = 0.99)
		
		mc_wgcna_gmods = gmod_moduleNonoverlapingMembership(
			mc_wgcna, 
			mc_wgcna_me, 
			kME_threshold = 0.5)
		
		# save objects for later use
		saveRDS(mc_wgcna,       sprintf("%s/gmod_%s.wgcna_object.wgcna.rds", out_fn, spi))
		saveRDS(mc_wgcna_gmods, sprintf("%s/gmod_%s.wgcna_object.gmods.rds", out_fn, spi))
		saveRDS(mc_wgcna_me,    sprintf("%s/gmod_%s.wgcna_object.ME.rds", out_fn, spi))
		
	}

}

