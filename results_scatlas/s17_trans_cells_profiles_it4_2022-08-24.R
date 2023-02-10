# libraries
library("scales")
library("metacell")
source("../scripts/Downstream_functions.R")
source("../scripts/helper.R")
source("../scripts/geneSetAnalysis.R")
graphics.off()

# files
out_fn = "results_trans_analysis"
dir.create(out_fn)

# list of species
sps_list = c("Tadh","TrH2","Hhon","HoiH23")
fp_thr = 2

# orthogroups
message(sprintf("data | load orthology"))
og_group_fn = "../data/results_broccoli_ml/orthogroup_conservation.csv"
ogg = read.table(og_group_fn, sep = "\t", header = TRUE)
ogg_v = ogg$orthogroup
names(ogg_v) = ogg$gene
ogv_v = ogg$gene
names(ogv_v) = ogg$orthogroup

### Trans profiles with AUCell ###

# trans definition
focus_peptidergic = c("peptidergic_alpha","peptidergic_beta","peptidergic_gamma","peptidergic_delta","peptidergic_epsilon","peptidergic_zeta","peptidergic_eta","peptidergic_theta","peptidergic_iota","peptidergic_kappa","peptidergic_lambda","peptidergic_mu","peptidergic_nu","peptidergic_omicron")

trans_pairs_list = list(
	Tadh = list(
		trans_lipophil_1_gland = list("trans_lipophil_1_gland", "lipophil", "gland_1"),
		trans_lipophil_1_epithelia = list("trans_lipophil_1_epithelia", "lipophil", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like")),
		trans_fibre_gland = list("trans_fibre_gland", "fibre", "gland_1"),
		trans_gland_epithelia = list("trans_gland_epithelia", "gland_1", c("epithelia_dorsal","epithelia_dorsal_like")),
		trans_epithelia_dorsal_ventral = list("trans_epithelia_dorsal_ventral", "epithelia_ventral", "epithelia_dorsal"),
		peptidergic_progenitors = list("peptidergic_progenitors", "epithelia_ventral", focus_peptidergic),
		epithelia_gland_like = list("epithelia_gland_like", "gland_1", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like"))
	),
	TrH2 = list(
		trans_lipophil_1_fibre = list("trans_lipophil_1_fibre", "lipophil", "fibre"),
		trans_lipophil_1_gland = list("trans_lipophil_1_gland", "lipophil", "gland_1"),
		trans_fibre_epithelia = list("trans_fibre_epithelia", "fibre", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like")),
		trans_epithelia_ventral_peptidergic = list("trans_epithelia_ventral_peptidergic", "epithelia_ventral", "peptidergic_iota"),
		trans_epithelia_dorsal_ventral = list("trans_epithelia_dorsal_ventral", "epithelia_ventral", c("epithelia_dorsal","epithelia_dorsal_like")),
		peptidergic_progenitors = list("peptidergic_progenitors", "epithelia_ventral", focus_peptidergic),
		epithelia_gland_like = list("epithelia_gland_like", "gland_1", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like"))
	),
	Hhon = list(
		trans_lipophil_1_gland = list("trans_lipophil_1_gland", "lipophil_1", "gland_1"),
		trans_lipophil_1_epithelia = list("trans_lipophil_1_epithelia", "lipophil_1", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like")),
		trans_lipophil_2_epithelia = list("trans_lipophil_2_epithelia", "lipophil_2", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like")),
		trans_fibre_gland = list("trans_fibre_gland", "fibre", "gland_1"),
		trans_gland_epithelia_ventral = list("trans_gland_epithelia_ventral", "gland_1", "epithelia_ventral"),
		trans_gland_epithelia_dorsal = list("trans_gland_epithelia_dorsal",  "gland_1", "epithelia_dorsal"),
		trans_epithelia_dorsal_ventral = list("trans_epithelia_dorsal_ventral",  "epithelia_ventral", "epithelia_dorsal"),
		peptidergic_progenitors = list("peptidergic_progenitors", "epithelia_ventral", focus_peptidergic),
		epithelia_gland_like = list("epithelia_gland_like", "gland_1", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like"))
	),
	HoiH23 = list(
		trans_lipophil_1_fibre = list("trans_lipophil_1_fibre", "lipophil_1", "fibre_1"),
		trans_lipophil_1_fibre = list("trans_lipophil_1_gland", "lipophil_1", "gland_1"),
		trans_lipophil_1_epithelia = list("trans_lipophil_1_epithelia", "lipophil_1", c("epithelia_ventral")),
		trans_lipophil_1_lipophil_2 = list("trans_lipophil_1_lipophil_2", "lipophil_1", "lipophil_2"),
		trans_lipophil_2_epithelia = list("trans_lipophil_2_epithelia", "lipophil_2", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like")),
		trans_lipophil_2_gland = list("trans_lipophil_2_gland", "lipophil_2", "gland_1"),
		trans_fibre_gland = list("trans_fibre_gland", c("fibre_1","fibre_2"), "gland_1"),
		trans_gland_epithelia_ventral = list("trans_gland_epithelia_ventral", "gland_1", "epithelia_ventral"),
		trans_gland_epithelia_dorsal = list("trans_gland_epithelia_dorsal", "gland_1", "epithelia_dorsal"),
		peptidergic_progenitors = list("peptidergic_progenitors", "epithelia_ventral", focus_peptidergic),
		epithelia_gland_like = list("epithelia_gland_like", "gland_1", c("epithelia_ventral","epithelia_dorsal","epithelia_dorsal_like"))
	)
)

# empty list for markers
trans_pairs_list_ogmarkers_i = list()
trans_pairs_list_ogmarkers_j = list()
trans_pairs_list_ogmarkers_t = list()

# loop species
for (spi in names(trans_pairs_list)) {
	
	# init database
	message(sprintf("trans profiles AUCell %s | AU ranks...", spi))
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	
	# first load reordered recluster mc solution
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it4",run_name))
	ct = metacell::scdb_mc(sprintf("%s_it4_cts",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it2",run_name))
	
    # metacell and single cell annots
    ctt_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
    ctt = read.table(ctt_fn, header = TRUE, sep = "\t", comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$broad_cell_type = factor(ctt$broad_cell_type, levels = unique(ctt$broad_cell_type))

	# color dictionaries
	dict_ctmc = as.character(ctt$cell_type)
	names(dict_ctmc) = ctt$metacell
	dict_ctco = as.character(ctt$color)
	names(dict_ctco) = ctt$cell_type
	dict_mcco = ctt$color
	names(dict_mcco) = ctt$metacell
	dict_scco = dict_mcco [ mc@mc ]
	names(dict_scco) = names(mc@mc)
	
	# reorder
	ct@mc_fp = ct@mc_fp [ , unique(ctt$cell_type) ]
	
	# load TFs
	set_fn = sprintf("../data/gene_annotations/tfs.%s_genes.curated.csv", spi)
	set = read.table(set_fn, col.names = c("transcript","annot"))
	set$gene = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), set$transcript)
	

	# for AUC based on umifrac, mc level (test)
	mc_counts = sca_mc_gene_counts(mc_object = mc, mat_object = mat)
	mc_umifra = sca_mc_gene_umifrac(mc, mc_counts)
	# for AUC based on umifrac, ct level (train)
	ctt_v = ctt$cell_type
	names(ctt_v) = ctt$metacell
	scs_v = as.character(ctt_v) [ mc@mc ]
	names(scs_v) = names(mc@mc)
	ct_counts = sca_mc_gene_counts(mc_object = mc, mat_object = mat, grouping_vector = scs_v)
	ct_umifra = sca_mc_gene_umifrac_noobj(mc_counts = ct_counts)

	# AUCell ranks
	auc_fp_cellrank = AUCell::AUCell_buildRankings(mc@mc_fp, nCores=4, plotStats=FALSE, verbose = FALSE)
	# auc_uf_cellrank = AUCell::AUCell_buildRankings(mc_umifra, nCores=4, plotStats=FALSE, verbose = FALSE)

	# trans pairs to evaluate
	trans_pairs_i = trans_pairs_list[[spi]]
	trans_pairs_markers = list()
	
	### Trans profiles based on umifracs, footprints, umifracsum ###
	
	pdf(sprintf("%s/trans_profile_AU.%s.pdf", out_fn, spi), width = 8, height = 30)
	layout(matrix(1:48, byrow = TRUE, ncol = 4))
	
	for (pai in trans_pairs_i) {
		
		tra = pai[[1]]
		rli = pai[[2]]
		rlj = pai[[3]]
		message(sprintf("trans profiles AUCell %s | AU score umifrac, %s", spi, tra))

		# mcs in trans query
		mcs_i = ctt [ ctt$cell_type %in% rli , "metacell" ]
		mcs_j = ctt [ ctt$cell_type %in% rlj , "metacell" ]
		mcs_t = ctt [ ctt$cell_type %in% tra , "metacell" ]

		
		## Footprint based markers ##
		# markers in cell type i and j
		mks_i = unique(names(which(apply(as.data.frame(ct@mc_fp [ , intersect(rli, colnames(ct@mc_fp)) ]), 1, function(r) max(r) >= fp_thr))))
		mks_j = unique(names(which(apply(as.data.frame(ct@mc_fp [ , intersect(rlj, colnames(ct@mc_fp)) ]), 1, function(r) max(r) >= fp_thr))))
		mks_t = unique(names(which(apply(as.data.frame(ct@mc_fp [ , intersect(tra, colnames(ct@mc_fp)) ]), 1, function(r) max(r) >= fp_thr))))
		# mks_i = mks_i [ !mks_i %in% mks_j ]
		# mks_j = mks_j [ !mks_j %in% mks_i ]
		
		# # remove antimarkers
		# mka_i = unique(names(which(apply(as.data.frame(ct@mc_fp [ , colnames(ct@mc_fp) [ !colnames(ct@mc_fp) %in% c(rli, colnames(ct@mc_fp) [ grep("trans_",colnames(ct@mc_fp)) ]) ] ]), 1, function(r) max(r) >= fp_thr))))
		# mka_j = unique(names(which(apply(as.data.frame(ct@mc_fp [ , colnames(ct@mc_fp) [ !colnames(ct@mc_fp) %in% c(rlj, colnames(ct@mc_fp) [ grep("trans_",colnames(ct@mc_fp)) ]) ] ]), 1, function(r) max(r) >= fp_thr))))
		# mks_i = mks_i [ ! mks_i %in% mka_i ]
		# mks_j = mks_j [ ! mks_j %in% mka_j ]
		
		# list of markers for AUCell		
		mks_l = list(mks_i, mks_j)
		names(mks_l) = c(rli[1], rlj[1])

		# plot AU based on footprint
		auc_fp_auc = AUCell::AUCell_calcAUC(geneSets = mks_l, rankings = auc_fp_cellrank, normAUC = TRUE, aucMaxRank = ceiling(1 * nrow(auc_fp_cellrank)), verbose = FALSE)
		plot(auc_fp_auc@assays@data$AUC[1,], auc_fp_auc@assays@data$AUC[2,], col = ctt$color, pch = 19, xlim = c(0.15, 1), ylim = c(0.15, 1), xlab = sprintf("AU %s markers", rli[1]), ylab = sprintf("AU %s markers", rlj[1]), main = sprintf("AUCell fp markers\n%s", tra), las = 1, cex.main = 0.7, cex = 0.7, cex.axis = 0.7, cex.lab = 0.7)
		text(auc_fp_auc@assays@data$AUC[1,mcs_t], auc_fp_auc@assays@data$AUC[2,mcs_t], tra, cex = 0.7, col = scales::alpha("thistle4", 0.8))
		abline(a = 0, b = 1, h= 0.5, v= 0.5, lty = 2)

		# plot UMIfrac sum, fp-based markers
		ufs_i = colSums(mc_umifra[mks_i,])
		ufs_j = colSums(mc_umifra[mks_j,])
		plot(ufs_i, ufs_j, col = ctt$color, pch = 19, xlab = sprintf("Sum UMIfrac %s markers", rli[1]), ylab = sprintf("Sum UMIfrac %s markers", rlj[1]), main = sprintf("AUCell fp markers, sumUMIfrac\n%s", tra), las = 1, cex.main = 0.7, cex = 0.7, cex.axis = 0.7, cex.lab = 0.7, xlim = c(0,600), ylim = c(0,600))
		text(ufs_i[mcs_t], ufs_j[mcs_t], tra, cex = 0.7, col = scales::alpha("thistle4", 0.8))
		abline(a = 0, b = 1, lty = 2)

		# plot UMIfrac sum, fp-based markers
		ufs_i = apply(mc@mc_fp[mks_i,], 2, mean)
		ufs_j = apply(mc@mc_fp[mks_j,], 2, mean)
		plot(ufs_i, ufs_j, col = ctt$color, pch = 19, xlab = sprintf("Mean fp %s markers", rli[1]), ylab = sprintf("Mean fp %s markers", rlj[1]), main = sprintf("AUCell fp markers, mean fp\n%s", tra), las = 1, cex.main = 0.7, cex = 0.7, cex.axis = 0.7, cex.lab = 0.7, xlim = c(0.5,12), ylim = c(0.5,12), log = "xy")
		text(ufs_i[mcs_t], ufs_j[mcs_t], tra, cex = 0.7, col = scales::alpha("thistle4", 0.8))
		abline(a = 0, b = 1, lty = 2)
		
		# store markers 
		trans_pairs_markers[[tra]] = list()
		trans_pairs_markers[[tra]][[tra]] = mks_t
		trans_pairs_markers[[tra]][[rli[1]]] = mks_i
		trans_pairs_markers[[tra]][[rlj[1]]] = mks_j
		
		# save list of markers
		mks_i_og = ogg_v [ mks_i ]
		mks_j_og = ogg_v [ mks_j ]
		mks_t_og = ogg_v [ mks_t ]
		mks_i_og = mks_i_og [ !is.na(mks_i_og) ]
		mks_j_og = mks_j_og [ !is.na(mks_j_og) ]
		mks_t_og = mks_t_og [ !is.na(mks_t_og) ]
		trans_pairs_list_ogmarkers_i[[paste(tra,spi, sep = "|")]] = mks_i_og
		trans_pairs_list_ogmarkers_j[[paste(tra,spi, sep = "|")]] = mks_j_og
		trans_pairs_list_ogmarkers_t[[paste(tra,spi, sep = "|")]] = mks_t_og
		
	}
	
	dev.off()
	
	pdf(sprintf("%s/trans_profile_Venn.%s.pdf", out_fn, spi), width = 4, height = 4)
	for (pai in trans_pairs_i) {
	
		tra = pai[[1]]
		rli = pai[[2]]
		rlj = pai[[3]]
		message(sprintf("trans profiles Venn %s | overlaps %s", spi, tra))
		mks_t = trans_pairs_markers[[tra]][[1]]
		mks_i = trans_pairs_markers[[tra]][[2]]
		mks_j = trans_pairs_markers[[tra]][[3]]
		vv = venn.three(mks_t, mks_i, mks_j, catname1 = tra, catname2 = rli[1], catname3 = rlj[1], col1 = "blue", col2 = "gray60", col3 = "gray80", eulerbool = TRUE, main = paste(spi, tra))
		
	}
	dev.off()

	pdf(sprintf("%s/trans_profile_Expr.%s.pdf", out_fn, spi), width = 8, height = 30)
	layout(matrix(1:48, byrow = TRUE, ncol = 4))
	for (pai in trans_pairs_i) {
	
		tra = pai[[1]]
		rli = pai[[2]]
		rlj = pai[[3]]
		message(sprintf("trans profiles Venn %s | overlaps %s", spi, tra))
		mks_t = trans_pairs_markers[[tra]][[1]]
		mks_i = trans_pairs_markers[[tra]][[2]]
		mks_j = trans_pairs_markers[[tra]][[3]]
		mks_ij_not = base::setdiff(union(mks_i, mks_j), mks_t)
		mks_ij_yet = mks_t [ mks_t %in% mks_i | mks_t %in% mks_j ]
		mks_ij_not = mks_ij_not [ mks_ij_not %in% rownames(mat@mat) ]
		mks_ij_yet = mks_ij_yet [ mks_ij_yet %in% rownames(mat@mat) ]
		
		oo = list(log10(rowSums(as.matrix(mat@mat [ mks_ij_yet, ]))), log10(rowSums(as.matrix(mat@mat [ mks_ij_not, ])))  )		
		tt = wilcox.test(oo[[1]], oo[[2]])
		boxplot(oo, las = 2, names = c("shared", "terminal"), col = "lightblue3", ylab = "log10(tUMI)", ylim = c(0,6), outline = FALSE)
		title(main = sprintf("markers %s+%s\nshared with %s", rli[1], rlj[1], tra), cex.main = 0.7, sub = sprintf("wilcox %s\np=%.1E", tt$alternative, tt$p.value), col.sub = "gray10", col.main = "gray20", cex.sub = 0.7)
		
		
	}
	dev.off()
	


}

message("AUCell done")



# comparable trans
trans_compare_csps = list(
	trans_lipo1_gland = c("trans_lipophil_1_gland|Tadh","trans_lipophil_1_gland|TrH2","trans_lipophil_1_gland|Hhon","trans_lipophil_1_gland|HoiH23"),
	trans_lipo1_epith = c("trans_lipophil_1_epithelia|Tadh","trans_lipophil_1_epithelia|Hhon","trans_lipophil_1_epithelia|HoiH23"),
	trans_lipo1_fibre = c("trans_lipophil_1_fibre|TrH2","trans_lipophil_1_fibre|HoiH23"),
	trans_fibre_gland = c("trans_fibre_gland|Tadh","trans_fibre_gland|Hhon","trans_fibre_gland|HoiH23"),
	trans_lipo2_epith = c("trans_lipophil_2_epithelia|Hhon","trans_lipophil_2_epithelia|HoiH23"),
	trans_gland_epive = c("trans_gland_epithelia_ventral|Hhon","trans_gland_epithelia_ventral|HoiH23","trans_gland_epithelia|Tadh"),
	trans_gland_epido = c("trans_gland_epithelia_dorsal|Hhon","trans_gland_epithelia_dorsal|HoiH23","trans_gland_epithelia|Tadh"),
	trans_epive_epido = c("trans_epithelia_dorsal_ventral|Hhon","trans_epithelia_dorsal_ventral|Tadh","trans_epithelia_dorsal_ventral|TrH2")
)


for (n in 1:length(trans_compare_csps)) {
	
	tcn = names(trans_compare_csps)[n]
	tcl = trans_compare_csps[[tcn]]
	
	# names for terminal i and j
	tcn_i = stringr::str_split(tcn,"_")[[1]][2]
	tcn_j = stringr::str_split(tcn,"_")[[1]][3]
	tcn_sps = stringr::str_split(trans_compare_csps[[tcn]],"\\|", simplify = TRUE)[,2]
	tcn_tra = stringr::str_split(trans_compare_csps[[tcn]],"\\|", simplify = TRUE)[,1]
	tcn_sps_i_v = paste(tcn_sps, tcn_i, sep = "|")
	tcn_sps_j_v = paste(tcn_sps, tcn_j, sep = "|")

	pdf(sprintf("%s/trans_profile_Venn.csps.trans.%s.pdf", out_fn,tcn), width = 4, height = 4)
	if (length(tcl) == 2) {
		
		vt = venn.two(
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][2]]],
			catname1 = trans_compare_csps[[tcn]][1],
			catname2 = trans_compare_csps[[tcn]][2],
			eulerbool = TRUE, main = tcn)
		vi = venn.two(
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][2]]],
			catname1 = tcn_sps_i_v[1],
			catname2 = tcn_sps_i_v[2],
			eulerbool = TRUE, main = tcn)
		vj = venn.two(
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][2]]],
			catname1 = tcn_sps_j_v[1],
			catname2 = tcn_sps_j_v[2],
			eulerbool = TRUE, main = tcn)
			
		# binomial test
		v_sha_t = vt$list_intersect
		n_sha_t = length(vt$list_intersect)
		n_uni_t = length(vt$union)
		f_sha_t = length(vt$list_intersect) / length(vt$union)
		f_sha_i = length(vi$list_intersect) / length(vi$union)
		f_sha_j = length(vj$list_intersect) / length(vj$union)
		# f_sha_a = mean(c(f_sha_i, f_sha_j))
		f_sha_a = (length(vi$list_intersect) + length(vj$list_intersect)) / (length(vi$union) + length(vj$union))
		n_exp_t_i = n_uni_t * f_sha_i
		n_exp_t_j = n_uni_t * f_sha_j
		n_exp_t_a = n_uni_t * f_sha_a
		btm = binom.test(x = n_sha_t, n = n_uni_t, p = mean(c(f_sha_i, f_sha_j)), alternative = "less")
		bti = binom.test(x = n_sha_t, n = n_uni_t, p = f_sha_i, alternative = "less")
		btj = binom.test(x = n_sha_t, n = n_uni_t, p = f_sha_j, alternative = "less")
		
	} else if (length(tcl) == 3) {
		
		vt = venn.three(
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][2]]],
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][3]]],
			catname1 = trans_compare_csps[[tcn]][1],
			catname2 = trans_compare_csps[[tcn]][2],
			catname3 = trans_compare_csps[[tcn]][3],
			eulerbool = TRUE, main = tcn)
		vi = venn.three(
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][2]]],
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][3]]],
			catname1 = tcn_sps_i_v[1],
			catname2 = tcn_sps_i_v[2],
			catname3 = tcn_sps_i_v[3],
			eulerbool = TRUE, main = tcn)
		vj = venn.three(
			trans_pairs_list_ogmarkers_j[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_j[[trans_compare_csps[[tcn]][2]]],
			trans_pairs_list_ogmarkers_j[[trans_compare_csps[[tcn]][3]]],
			catname1 = tcn_sps_j_v[1],
			catname2 = tcn_sps_j_v[2],
			catname3 = tcn_sps_j_v[3],
			eulerbool = TRUE, main = tcn)
			
		# binomial test
		v_sha_t = vt$intersect_123
		n_sha_t = length(vt$intersect_123)
		n_uni_t = length(vt$union)
		f_sha_t = length(vt$intersect_123) / length(vt$union)
		f_sha_i = length(vi$intersect_123) / length(vi$union)
		f_sha_j = length(vj$intersect_123) / length(vj$union)
		# f_sha_a = mean(c(f_sha_i, f_sha_j))
		f_sha_a = (length(vi$intersect_123) + length(vj$intersect_123)) / (length(vi$union) + length(vj$union))
		n_exp_t_i = n_uni_t * f_sha_i
		n_exp_t_j = n_uni_t * f_sha_j
		n_exp_t_a = n_uni_t * f_sha_a
		btm = binom.test(x = n_sha_t, n = n_uni_t, p = mean(c(f_sha_i, f_sha_j)), alternative = "less")
		bti = binom.test(x = n_sha_t, n = n_uni_t, p = f_sha_i, alternative = "less")
		btj = binom.test(x = n_sha_t, n = n_uni_t, p = f_sha_j, alternative = "less")
		
	} else if (length(tcl) == 4) {
		
		vt = venn.four(
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][2]]],
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][3]]],
			trans_pairs_list_ogmarkers_t[[trans_compare_csps[[tcn]][4]]],
			catname1 = trans_compare_csps[[tcn]][1],
			catname2 = trans_compare_csps[[tcn]][2],
			catname3 = trans_compare_csps[[tcn]][3],
			catname4 = trans_compare_csps[[tcn]][4],
			eulerbool = TRUE, main = tcn)
		vi = venn.four(
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][2]]],
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][3]]],
			trans_pairs_list_ogmarkers_i[[trans_compare_csps[[tcn]][4]]],
			catname1 = tcn_sps_i_v[1],
			catname2 = tcn_sps_i_v[2],
			catname3 = tcn_sps_i_v[3],
			catname4 = tcn_sps_i_v[4],
			eulerbool = TRUE, main = tcn)
		vj = venn.four(
			trans_pairs_list_ogmarkers_j[[trans_compare_csps[[tcn]][1]]],
			trans_pairs_list_ogmarkers_j[[trans_compare_csps[[tcn]][2]]],
			trans_pairs_list_ogmarkers_j[[trans_compare_csps[[tcn]][3]]],
			trans_pairs_list_ogmarkers_j[[trans_compare_csps[[tcn]][4]]],
			catname1 = tcn_sps_j_v[1],
			catname2 = tcn_sps_j_v[2],
			catname3 = tcn_sps_j_v[3],
			catname4 = tcn_sps_j_v[4],
			eulerbool = TRUE, main = tcn)
		
		# binomial test
		v_sha_t = vt$intersect_1234
		n_sha_t = length(vt$intersect_1234)
		n_uni_t = length(vt$union)
		f_sha_t = length(vt$intersect_1234) / length(vt$union)
		f_sha_i = length(vi$intersect_1234) / length(vi$union)
		f_sha_j = length(vj$intersect_1234) / length(vj$union)
		# f_sha_a = mean(c(f_sha_i, f_sha_j))
		f_sha_a = (length(vi$intersect_1234) + length(vj$intersect_1234)) / (length(vi$union) + length(vj$union))
		n_exp_t_i = n_uni_t * f_sha_i
		n_exp_t_j = n_uni_t * f_sha_j
		n_exp_t_a = n_uni_t * f_sha_a
		btm = binom.test(x = n_sha_t, n = n_uni_t, p = mean(c(f_sha_i, f_sha_j)), alternative = "less")
		bti = binom.test(x = n_sha_t, n = n_uni_t, p = f_sha_i, alternative = "less")
		btj = binom.test(x = n_sha_t, n = n_uni_t, p = f_sha_j, alternative = "less")
		
	}
	
	# plot observed/expected fractions
	par(mar = c(5.1, 12, 4.1, 2.1))
	b = barplot(
		c(f_sha_t, f_sha_a, f_sha_i, f_sha_j),
		names.arg = c("shared trans", "expected trans (all)", sprintf("expected trans (%s)", tcn_i), sprintf("expected trans (%s)", tcn_j)),
		horiz = TRUE, las = 2, col = c("darkolivegreen3","snow4","snow3","snow3"),
		beside = FALSE, ylim = c(10,0),
		border = NA, xlab = "# genes", xlim = c(0,0.5)
	)
	title(main = sprintf("trans O/E in %s", tcn), cex.main = 1, sub = sprintf("trans O/E binomial test\nn=%i genes in trans", length(n_sha_t)), cex.sub = 0.5)
	v_pval = c(
		"",
		sprintf("p=%.1E | o/e=%i/%.1f", btm$p.value, n_sha_t, n_exp_t_a),
		sprintf("p=%.1E | o/e=%i/%.1f", bti$p.value, n_sha_t, n_exp_t_i),
		sprintf("p=%.1E | o/e=%i/%.1f", btj$p.value, n_sha_t, n_exp_t_j)
	)
	text(0, b, v_pval, adj = 0, col = "gray10", cex = 0.7)


	# fract
	par(mar = c(5.1, 4.1, 4.1, 2.1))
	suppressMessages(metacell::scdb_init("data/scdb/",force_reinit=TRUE))
	
	# first load reordered recluster mc solution
	run_name = sprintf("scdr_%s", tcn_sps[1])
	mc = metacell::scdb_mc(sprintf("%s_it4",run_name))
	ct = metacell::scdb_mc(sprintf("%s_it4_cts",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it2",run_name))
	# for AUC based on umifrac, mc level (test)
	mc_counts = sca_mc_gene_counts(mc_object = mc, mat_object = mat)
	mc_umifra = sca_mc_gene_umifrac(mc, mc_counts)

	
    # metacell and single cell annots
    ctt_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", tcn_sps[1])
    ctt = read.table(ctt_fn, header = TRUE, sep = "\t", comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$broad_cell_type = factor(ctt$broad_cell_type, levels = unique(ctt$broad_cell_type))
	dict_ctmc = as.character(ctt$cell_type)
	names(dict_ctmc) = ctt$metacell
	dict_ctco = as.character(ctt$color)
	names(dict_ctco) = ctt$cell_type
	dict_mcco = ctt$color
	names(dict_mcco) = ctt$metacell
	dict_scco = dict_mcco [ mc@mc ]
	names(dict_scco) = names(mc@mc)
	
	# reorder
	ct@mc_fp = ct@mc_fp [ , unique(ctt$cell_type) ]

	# plot umifracs in terminal, highlighting shared OGs
	rli = trans_pairs_list[[tcn_sps[1]]][[tcn_tra[1]]][[2]]
	rlj = trans_pairs_list[[tcn_sps[1]]][[tcn_tra[1]]][[3]]
	mcs_i = ctt [ ctt$cell_type %in% rli, "metacell" ] 
	mcs_j = ctt [ ctt$cell_type %in% rlj , "metacell" ] 
	mks_t = unique(names(which(apply(as.data.frame(ct@mc_fp [ , intersect(tcn_tra[1], colnames(ct@mc_fp)) ]), 1, function(r) max(r) >= fp_thr))))
	mks_t_shared = mks_t [ mks_t %in% ogv_v [ v_sha_t ] ]
	ufs_t_i = sort(rowSums(mc_umifra[mks_t,mcs_i]))
	ufs_t_j = sort(rowSums(mc_umifra[mks_t,mcs_j]))
	col_ufs_t_i = rep(alpha(dict_ctco [ rli[1] ], 0.6), length(ufs_t_i))
	col_ufs_t_j = rep(alpha(dict_ctco [ rlj[1] ], 0.6), length(ufs_t_j))
	col_ufs_t_i [ names(ufs_t_i) %in%  mks_t_shared] = "gray10"
	col_ufs_t_j [ names(ufs_t_j) %in%  mks_t_shared] = "gray10"
	plot(ufs_t_i, log = "y", las = 2, pch = 19,  col = col_ufs_t_i, ylim = c(0.001,1e4), ylab = "sum UMIfrac", xlab = "")
	points(ufs_t_j, pch = 19, col = col_ufs_t_j)
	# boxplot
	wti = wilcox.test(ufs_t_i [ !names(ufs_t_i) %in%  mks_t_shared ], ufs_t_i [ names(ufs_t_i) %in%  mks_t_shared ])
	wtj = wilcox.test(ufs_t_j [ !names(ufs_t_j) %in%  mks_t_shared ], ufs_t_j [ names(ufs_t_j) %in%  mks_t_shared ])
	boxplot(
		ufs_t_i [ !names(ufs_t_i) %in%  mks_t_shared ], 
		ufs_t_i [ names(ufs_t_i) %in%  mks_t_shared ], 
		ufs_t_j [ !names(ufs_t_j) %in%  mks_t_shared ],
		ufs_t_j [ names(ufs_t_j) %in%  mks_t_shared ],
		outline = FALSE,
		las = 1, ylab = "sum UMIfrac",
		names = c(rli[1], "shared", rlj[1], "shared"),
		col = c(dict_ctco [ rli[1] ], "gray", dict_ctco [ rlj[1] ], "gray"),
		log = "y", ylim = c(0.001,1e4))
	title(sub = sprintf("wilcox-2s trans shared v %s: p=%.1E | v %s shared p=%.1E\n%i trans genes, %i shared", rli[1], wti$p.value, rlj[1], wtj$p.value, length(mks_t), length(mks_t_shared)), col = "gray20", cex.sub = 0.5)
	dev.off()
	
	# save shared OGs
	ogg_i = ogg [ ogg$orthogroup %in% v_sha_t, ]
	write.table(ogg_i, sprintf("%s/trans_profile_Venn.csps.trans.%s.csv", out_fn,tcn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


}




# cross-species overlap in trans markers
jac_csps_m_i  = matrix(nrow = length(trans_pairs_list_ogmarkers_i), ncol = length(trans_pairs_list_ogmarkers_j))
jac_csps_m_j  = matrix(nrow = length(trans_pairs_list_ogmarkers_j), ncol = length(trans_pairs_list_ogmarkers_j))
jac_csps_m_t  = matrix(nrow = length(trans_pairs_list_ogmarkers_t), ncol = length(trans_pairs_list_ogmarkers_t))
jac_csps_m_ij = matrix(nrow = length(trans_pairs_list_ogmarkers_t), ncol = length(trans_pairs_list_ogmarkers_t))
rownames(jac_csps_m_i) = rownames(jac_csps_m_j) = rownames(jac_csps_m_t) = rownames(jac_csps_m_ij) = names(trans_pairs_list_ogmarkers_t)
colnames(jac_csps_m_i) = colnames(jac_csps_m_j) = colnames(jac_csps_m_t) = colnames(jac_csps_m_ij) = names(trans_pairs_list_ogmarkers_t)
for (i in 1:length(trans_pairs_list_ogmarkers_t)) {
	for (j in 1:length(trans_pairs_list_ogmarkers_t)) {
		ni = names(trans_pairs_list_ogmarkers_t)[i]
		nj = names(trans_pairs_list_ogmarkers_t)[j]
		jac_csps_m_i [ni,nj] = jaccard_index(trans_pairs_list_ogmarkers_i[[ni]], trans_pairs_list_ogmarkers_i[[nj]])
		jac_csps_m_j [ni,nj] = jaccard_index(trans_pairs_list_ogmarkers_j[[ni]], trans_pairs_list_ogmarkers_j[[nj]])
		jac_csps_m_t [ni,nj] = jaccard_index(trans_pairs_list_ogmarkers_t[[ni]], trans_pairs_list_ogmarkers_t[[nj]])
		jac_csps_m_ij[ni,nj] = jaccard_index(trans_pairs_list_ogmarkers_i[[ni]], trans_pairs_list_ogmarkers_j[[nj]])
	}
}

rownames(jac_csps_m_ij) = names(trans_pairs_list_ogmarkers_i)
colnames(jac_csps_m_ij) = names(trans_pairs_list_ogmarkers_j)

jac_csps_m_d_it = jac_csps_m_t - jac_csps_m_i

# general heatmap?
hm = plot_complex_heatmap(
	jac_csps_m_t [ !grepl("like|progenitors", rownames(jac_csps_m_t)), !grepl("like|progenitors", colnames(jac_csps_m_t)) ], 
	color_mat = c("darkred","red","orange","white","#d6e72e","#6fb600","#003f4d"),
	color_min = -1, color_max = 1, fontsize = 5,
	cluster_col = TRUE, cluster_row = TRUE, name = "Jaccard")
pdf(sprintf("%s/trans_profile_AU.csps.pdf", out_fn, spi), width = 6, height = 6)
print(hm)
dev.off()



# vectors of top markers per trans cell type per species (only among other species)
for (spi in sps_list) {
	
	pdf(sprintf("%s/trans_profile_Jaccard.%s.pdf", out_fn, spi), width = 12, height = 30)
	layout(matrix(1:48, byrow = TRUE, ncol = 4))
	tra_spi = rownames(jac_csps_m_i) [ grepl(spi, rownames(jac_csps_m_i)) ]
	tra_spo = rownames(jac_csps_m_i) [ !grepl(spi, rownames(jac_csps_m_i)) ]
	
	for (tri in tra_spi) {
		
		vo = head(sort(jac_csps_m_t[tri, tra_spo ], decreasing = TRUE), 20)
		barplot(vo, las = 2, pch = 19, col = "lightblue3", ylim = c(0,1), cex.names = 0.5)
		title(main = sprintf("%s top pairs",tri), cex.main = 1)
		
	}
	
	dev.off()

}




### Observed/expected trans cells ###

# loop species
pdf(sprintf("%s/trans_obs_exp.pdf", out_fn, spi), width = 4, height = 16)
layout(matrix(1:4, byrow = FALSE, ncol = 1))
par(mar = c(5.1, 12, 4.1, 2.1))
for (spi in names(trans_pairs_list)) {

    # metacell and single cell annots
    ctt_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
    sct_fn = sprintf("results_metacell_it4/scdr_%s.matrix.sc_annot.csv", spi)

    # load
    ctt = read.table(ctt_fn, header = TRUE, sep = "\t", comment.char = "")
    sct = read.table(sct_fn, header = TRUE, sep = "\t", comment.char = "")
	
	# factors
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$broad_cell_type = factor(ctt$broad_cell_type, levels = unique(ctt$broad_cell_type))
	sct$cell_type = factor(sct$cell_type, levels = levels(ctt$cell_type))
	
	# colors per ct
	cts_col = ctt [ , c("cell_type","color") ]
	cts_col = cts_col [ !duplicated(cts_col$cell_type) , ] 
	bct_col = ctt [ , c("broad_cell_type","color") ]
	bct_col = bct_col [ !duplicated(bct_col$broad_cell_type) , ] 

	# merge
	sct = merge(sct, ctt[,c("metacell","broad_cell_type")], by = "metacell", all.x  = TRUE, all.y = TRUE)

	trans_pairs_list_o = data.frame()
	for (pai in trans_pairs_list[[spi]]) {
	
		tra = pai[[1]]
		rli = pai[[2]]
		rlj = pai[[3]]
		
		scn_i = length(sct [ sct$cell_type %in% rli, "cell" ])
		scn_j = length(sct [ sct$cell_type %in% rlj, "cell" ])
		scn_t_obs = length(sct [ sct$cell_type %in% tra, "cell" ])
		scf_t_obs = scn_t_obs / nrow(sct)
		scf_i = scn_i / nrow(sct)
		scf_j = scn_j / nrow(sct)
		scf_t_exp = scf_i * scf_j
		scn_t_exp = scf_i * scf_j * nrow(sct)
		
		
		dd = data.frame(scn_t_obs, scn_t_exp, scn_i, scn_j, nrow(sct))
		colnames(dd) = c("trans_obs","trans_exp","terminal_i_obs","terminal_j_obs", "terminal_tot_obs")
		rownames(dd) = tra
		
		tt = binom.test(x = dd$trans_obs, n = dd$terminal_tot_obs, p = dd$terminal_i_obs/dd$terminal_tot_obs * dd$terminal_j_obs/dd$terminal_tot_obs )
		dd$pval = tt$p.value

		trans_pairs_list_o = rbind(trans_pairs_list_o, dd)

	}
	
	b = barplot(t(as.matrix(trans_pairs_list_o[,1:2])), names.arg = rownames(trans_pairs_list_o), horiz = TRUE, las = 2, col = c("darkolivegreen3","snow3"), beside = TRUE, ylim = c(50,0), border = NA, xlim = c(0,1500), xlab = "# cells")
	title(main = sprintf("trans O/E in %s", spi), cex.main = 1, sub = "trans O/E binomial test")
	v_pval = c()
	for(v in 1:length(trans_pairs_list_o[,"pval"])) {
		v_pval = c(v_pval, c("",sprintf("p=%.1E | o/e=%i/%.1f", trans_pairs_list_o[,"pval"][v], trans_pairs_list_o[,"trans_obs"][v], trans_pairs_list_o[,"trans_exp"][v])))
	}
	
	text(0, b, v_pval, adj = 0, col = "gray10", cex = 0.7)
	legend("bottomright", c("Observed","Expected"), fill = c("darkolivegreen3","snow3"), border = NA, bty = "n")

}
dev.off()



message("All done!")
