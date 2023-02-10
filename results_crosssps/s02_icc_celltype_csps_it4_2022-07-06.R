#### Input ####

# libraries
source("../scripts/helper.R")

# input
mdb_fn = "../results_scatlas/data/scdb/"
heatmap_colors = c("white","orange","orangered2","#520c52")

# list of comparisons
list_comp = list(
	c("Tadh","TrH2"),
	c("Tadh","Hhon"),
	c("Tadh","HoiH23"),
	c("TrH2","Hhon"),
	c("TrH2","HoiH23"),
	c("Hhon","HoiH23")
)

# list of sets to analyse (global and peptidergic)
set_list = list(
	global =      list(ann_fn = "../results_scatlas/results_metacell_it4/"            , icc_fn = "results_alignment_icc/", out_fn = "results_alignment_icc/",             mc_sprintf_string = "scdr_%s_it4",         mat_sprintf_string = "scdr_%s_it2",     ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv",   focus_list = list("cts" = "cell_type", "mcs" = "metacell", "bct" = "broad_cell_type")),
	peptidergic = list(ann_fn = "../results_scatlas/results_metacell_it4_peptidergic/", icc_fn = "results_alignment_icc/", out_fn = "results_alignment_icc_peptidergic/", mc_sprintf_string = "scdr_%s_it4_pep_ord", mat_sprintf_string = "scdr_%s_it4_pep", ctt_sprtinf_string = "%s/annotation_mc.%s.it4.peptidergic.tsv", focus_list = list("cts" = "cell_type", "mcs" = "metacell"))
)


for (set in set_list) {
	
	ann_fn = set$ann_fn
	icc_fn = set$icc_fn
	out_fn = set$out_fn
	focus_list = set$focus_list
	mc_sprintf_string  = set$mc_sprintf_string
	mat_sprintf_string = set$mat_sprintf_string
	ctt_sprtinf_string = set$ctt_sprtinf_string

	for (com in list_comp) {
		
		# define input
		sp1 = com[1]
		sp2 = com[2]
		
		# load cross-species marker data (from `s01`)
		csps_fn = sprintf("%s/csps_icc.mcs.%s-%s.rds", icc_fn, sp1, sp2)
		
		# load
		message(sprintf("csps icc | load %s-%s csps data...", sp1, sp2))
		csps = readRDS(csps_fn)
		
		# load cell type annotations
		message(sprintf("csps icc | load %s cell types...", sp1))
		ctt_sp1_fn = sprintf(ctt_sprtinf_string, ann_fn, sp1)
		ctt_sp1 = read.table(ctt_sp1_fn, header = TRUE, comment.char = "", sep = "\t")
		
		message(sprintf("csps icc | load %s cell types...", sp2))
		ctt_sp2_fn = sprintf(ctt_sprtinf_string, ann_fn, sp2)
		ctt_sp2 = read.table(ctt_sp2_fn, header = TRUE, comment.char = "", sep = "\t")

		# modify csps object to use cell types
		metacell::scdb_init(mdb_fn, force_reinit = TRUE)

		
		for (i in 1:length(focus_list)) {
			
			focid = names(focus_list)[[i]]
			focus = focus_list[[i]]

			# factors		
			ctt_sp1[,focus] = factor(as.character(ctt_sp1[,focus]), levels = unique(ctt_sp1[,focus]))
			ctt_sp2[,focus] = factor(as.character(ctt_sp2[,focus]), levels = unique(ctt_sp2[,focus]))


			# recalculate footprints
			message(sprintf("csps icc | %s footprints %s...", focus, sp1))
			csps_t_sp1 = sca_cell_type_fp(
				ctt_sp1[,c("metacell",focus,"color")], 
				mat_object = metacell::scdb_mat(sprintf(mat_sprintf_string, sp1)), 
				mc_object = metacell::scdb_mc(sprintf(mc_sprintf_string, sp1)))
			
			message(sprintf("csps icc | %s footprints %s...", focus, sp2))
			csps_t_sp2 = sca_cell_type_fp(
				ctt_sp2[,c("metacell",focus,"color")], 
				mat_object = metacell::scdb_mat(sprintf(mat_sprintf_string, sp2)), 
				mc_object = metacell::scdb_mc(sprintf(mc_sprintf_string, sp2)))
				
			# create new csps object, for metacells
			csps_t = csps
			csps_t$og_pairs = csps_t$og_pairs [ csps_t$og_pairs[,1] %in% rownames(csps_t_sp1@mc_fp) & csps_t$og_pairs[,2] %in% rownames(csps_t_sp2@mc_fp), ]
			csps_t$sp1 = csps_t_sp1@mc_fp[ csps_t$og_pairs[,1] , levels(ctt_sp1[,focus]) ]
			csps_t$sp2 = csps_t_sp2@mc_fp[ csps_t$og_pairs[,2] , levels(ctt_sp2[,focus]) ]
			csps_t$merged = quantile_normalisation(cbind(csps_t$sp1, csps_t$sp2))

			# save cell type csps
			saveRDS(csps_t, file = sprintf("%s/csps_icc.%s.%s-%s.rds", out_fn, focid, sp1, sp2))
			
		}
		
	}

}

message("All comparisons done!")

