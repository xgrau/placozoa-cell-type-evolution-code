#### Input ####

# libraries
source("../scripts/helper.R")

# input
mdb_fn = "../results_scatlas/data/scdb/"
heatmap_colors = c("white","orange","orangered2","#520c52")

# list of comparisons
list_comp = list(
	c("TrH2","Nvec"),
	c("TrH2","Hvul"),
	c("TrH2","Spis"),
	c("TrH2","Spolac"),
	c("TrH2","Mlei"),
	c("TrH2","Dmel"),
	c("TrH2","Mmus"),
	c("Tadh","Nvec"),
	c("Tadh","Hvul"),
	c("Tadh","Spis"),
	c("Tadh","Spolac"),
	c("Tadh","Mlei"),
	c("Tadh","Dmel"),
	c("Tadh","Mmus"),
	c("Hhon","Nvec"),
	c("Hhon","Hvul"),
	c("Hhon","Spis"),
	c("Hhon","Spolac"),
	c("Hhon","Mlei"),
	c("Hhon","Dmel"),
	c("Hhon","Mmus"),
	c("HoiH23","Nvec"),
	c("HoiH23","Hvul"),
	c("HoiH23","Spis"),
	c("HoiH23","Spolac"),
	c("HoiH23","Mlei"),
	c("HoiH23","Dmel"),
	c("HoiH23","Mmus")
)

# list of sets to analyse (global and peptidergic)
set_list = list(
	global = list(ann_fn = "../results_scatlas/results_metacell_it4/", icc_fn = "results_alignment_icc/", out_fn = "results_alignment_icc/", mc_sprintf_string = "scdr_%s_it4", mat_sprintf_string = "scdr_%s_it2", ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv",   focus_list = list("bct" = "broad_cell_type")),
	global = list(ann_fn = "../results_scatlas/results_metacell_it4/", icc_fn = "results_alignment_icc/", out_fn = "results_alignment_icc/", mc_sprintf_string = "scdr_%s_it4", mat_sprintf_string = "scdr_%s_it2", ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv",   focus_list = list("cts" = "cell_type"))
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
		ctt_sp2_fn = sprintf("../results_scatlas/data/scdb_outgroups/annot.%s.tsv", sp2)
		ctt_sp2 = read.table(ctt_sp2_fn, header = TRUE, comment.char = "", sep = "\t")

		
		for (i in 1:length(focus_list)) {
			
			focid = names(focus_list)[[i]]
			focus = focus_list[[i]]

			# factors		
			ctt_sp1[,focus] = factor(as.character(ctt_sp1[,focus]), levels = unique(ctt_sp1[,focus]))
			ctt_sp2[,focus] = factor(as.character(ctt_sp2[,focus]), levels = unique(ctt_sp2[,focus]))


			# recalculate footprints
			message(sprintf("csps icc | %s footprints %s...", focus, sp1))
			metacell::scdb_init(mdb_fn, force_reinit = TRUE)
			mat1 = metacell::scdb_mat(sprintf(mat_sprintf_string, sp1))
			mc1  = metacell::scdb_mc(sprintf(mc_sprintf_string, sp1))
			csps_t_sp1 = sca_cell_type_fp(
				ctt_sp1[,c("metacell",focus,"color")], 
				mat_object = mat1, 
				mc_object = mc1)
			
			message(sprintf("csps icc | %s footprints %s...", focus, sp2))
			metacell::scdb_init("../results_scatlas/data/scdb_outgroups", force_reinit = TRUE)
			mat2 = metacell::scdb_mat(sp2)
			mc2  = metacell::scdb_mc(sp2)
			csps_t_sp2 = sca_cell_type_fp(
				ctt_sp2[,c("metacell",focus,"color")], 
				mat_object = mat2, 
				mc_object = mc2)
				
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

