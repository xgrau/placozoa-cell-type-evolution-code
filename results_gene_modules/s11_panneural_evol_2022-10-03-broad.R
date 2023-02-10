#### Input ####

# libraries
library("stringr")
library("scales")
suppressMessages(source("../scripts/helper.R"))
library("ape")

# output folder
out_fn = "results_panneural_markers_broad/"

# reference species (for focused gene lists)
sps_ref  = c("Tadh","TrH2","Hhon","HoiH23")
sps_list = c("Tadh","TrH2","Hhon","HoiH23","Nvec","Spis","Hvul","Mmus","Dmel")

# list of genes (long format)
mgs_fn = "results_panneural_markers_broad/matrix.all.long.ogexpression.csv"
mgs_d  = read.table(mgs_fn, header = TRUE)

# orthogroup names (for interpretation only)
ngs_fn = "results_panneural_markers/matrix.all.ognames.csv"
ngs    = read.table(ngs_fn, header = TRUE)
ngs_v  = ngs$orthogroup_name
names(ngs_v) = ngs$orthogroup
ngs_v [ is.na(ngs_v) ] = names(ngs_v) [ is.na(ngs_v) ]

## load orthogroups ##

# orthogroups
message(sprintf("%s | prepare orthogroups", out_fn))
og_group_fn = "../data/results_broccoli_ml_extended_gmod/dir_step3/orthologous_groups.txt"
ogg = read.table(og_group_fn, sep = "\t", col.names = c("OG","genes"))
ogg_l = data.frame(
	OG   = unlist(sapply(1:nrow(ogg), function(i) { v = ogg[i,"genes"] ; vs = unlist(strsplit(v, " ")) ; lv = length(vs) ; vo = rep(ogg[i,"OG"], lv) })),
	gene = unlist(sapply(1:nrow(ogg), function(i) { v = ogg[i,"genes"] ; vs = unlist(strsplit(v, " ")) }))
)
ogg_l$species = gsub("_.*", "", ogg_l$gene)
ogg_l$species = factor(ogg_l$species, levels = sps_list)
ogg_l = ogg_l [ !is.na(ogg_l$species), ]

# clean spongilla
ogg_l [ ogg_l$species == "Spolac", "gene"] = gsub("_i\\d+_.*","", ogg_l [ ogg_l$species == "Spolac", "gene"])

# vector dictionary
ogg_v = ogg_l$OG
names(ogg_v) = ogg_l$gene


# load architectures
pfam_d = data.frame()
for (spi in c("Tadh","TrH2","Hhon","HoiH23")) {
	message(sprintf("data | pfam %s", spi))
	pfam_i = read.table(sprintf("../data/reference/%s_long.pep.annotations.csv", spi), sep = "\t", col.names = c("transcript","name","architecture"))
	pfam_i$species = spi
	pfam_d = rbind(pfam_d, pfam_i)
}
pfam_d_dict = pfam_d$architecture
names(pfam_d_dict) = pfam_d$transcript


#### Load expression data ####

mc_list     = list()
ct_list     = list()
ct_list_col = list()
for (spi in sps_list) {
	message(sprintf("load | %s", spi))
	if (spi %in% sps_ref) {
		mc_sprintf_string = "%s_it4"
		ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv"
		run_name = sprintf("scdr_%s", spi)
		focus_ct = "peptidergic"
		suppressMessages(metacell::scdb_init("../results_scatlas/data/scdb/",force_reinit=TRUE))
		mc_list[[spi]] = metacell::scdb_mc( sprintf(mc_sprintf_string, run_name))@mc_fp
		rownames(mc_list[[spi]]) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), rownames(mc_list[[spi]]), t2g = FALSE)
		ct_list[[spi]] = read.table(sprintf(ctt_sprtinf_string, "../results_scatlas/results_metacell_it4/", spi), header = TRUE, comment.char = "", sep = "\t")
		ct_list_col[[spi]] = ifelse(ct_list[[spi]]$broad_cell_type == focus_ct, "royalblue3", "snow2")
	} else {
		mc_sprintf_string = "%s"
		ctt_sprtinf_string = "annot.%s.csv"
		run_name = spi
		focus_ct = ifelse(spi == "Spolac", "neuroid", "neuron")
		suppressMessages(metacell::scdb_init("../results_scatlas/data/scdb_outgroups/",force_reinit=TRUE))
		mc_list[[spi]] = metacell::scdb_mc( sprintf(mc_sprintf_string, run_name))@mc_fp
		ct_list[[spi]] = read.table(sprintf("../results_scatlas/data/scdb_outgroups/annot.%s.tsv", spi), header = TRUE, comment.char = "", sep = "\t")
		ct_list_col[[spi]] = ifelse(ct_list[[spi]]$broad_cell_type == focus_ct, "blue", "snow2")
		if (spi %in% c("Mmus","Mlei")) {
			rownames(mc_list[[spi]]) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), rownames(mc_list[[spi]]), t2g = FALSE)
		} else if (spi == "Dmel") {
			rownames(mc_list[[spi]]) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), rownames(mc_list[[spi]]), t2g = FALSE, gene_id = "gene_name")
		}
	}
}


#### Loop over topologies ####

maxval = 400

# loop over species tree topologies
for (topi in c("CB","CP")) {

	# species tree
	message(sprintf("topology %s | load data", topi))
	phy_fn = sprintf("species_tree_%s.newick", topi)

	# input files, for expression
	msu_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogexpression.%s.dollo_summary.csv", topi)
	mpr_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogexpression.%s.dollo_pres.csv", topi)
	mga_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogexpression.%s.dollo_gain.csv", topi)
	mlo_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogexpression.%s.dollo_loss.csv", topi)

	# input files, for conservation
	csu_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogpresence.%s.dollo_summary.csv", topi)
	cpr_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogpresence.%s.dollo_pres.csv", topi)
	cga_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogpresence.%s.dollo_gain.csv", topi)
	clo_fn = sprintf("results_panneural_markers_broad/anc.all.long.ogpresence.%s.dollo_loss.csv", topi)

	# load data
	msu = read.table(msu_fn, sep = "\t", header = TRUE, row.names = 1)
	mpr = read.table(mpr_fn, sep = "\t", header = TRUE, row.names = 1)
	mga = read.table(mga_fn, sep = "\t", header = TRUE, row.names = 1)
	mlo = read.table(mlo_fn, sep = "\t", header = TRUE, row.names = 1)

	# load data
	csu = read.table(csu_fn, sep = "\t", header = TRUE, row.names = 1)
	cpr = read.table(cpr_fn, sep = "\t", header = TRUE, row.names = 1)
	cga = read.table(cga_fn, sep = "\t", header = TRUE, row.names = 1)
	clo = read.table(clo_fn, sep = "\t", header = TRUE, row.names = 1)

	# presence to binary
	mpr [ mpr > 1 ] = 1
	cpr [ cpr > 1 ] = 1

	#### Prepare species tree ####

	# read species tree
	phy = ape::read.tree(phy_fn)
	sps_list = phy$tip.label
	nod_list = c(phy$tip.label, phy$node.label)

	# prepare species dataframe
	# dataframe of edges
	phy_edge           = as.data.frame(phy$edge)
	colnames(phy_edge) = c("edge_start","edge_end")
	phy_edge$ix_edges = as.numeric(rownames(phy_edge))
	phy_edge$ends_in_tip = phy_edge$edge_end <= length(phy$tip.label)
	# dataframe of nodes
	phy_nods = data.frame(taxa = c(phy$tip.label, phy$node.label))
	phy_nods$edge_end = as.numeric(rownames(phy_nods))
	phy_nods$is_tip   = phy_nods$edge_end <= length(phy$tip.label)
	# merge them
	phy_edge = merge(phy_edge, phy_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")

	# which is the root node?
	root_label = phy$node.label[1]


	##### Prepare Dollo matrices #####

	# number of species in which each og is present
	msu$n_presences_sps = stringr::str_count(msu$presence, ",") + 1

	# filter out OGs present in less than three placozoans?
	# sps_plc = c("Tadh","TrH2","Hhon","HoiH23")
	# ogs_f = names(which(apply(mpr[ ,sps_plc ], 1, function(r) sum(r > 0) >= 2)))
	# filter out OGs gained in spongilla?
	# ogs_f = rownames(mga) [ which(sapply(mga[,"Spolac"], function(r) r == 0)) ]
	
	# filter out OGs present in only one species? (currently disabled)
	ogs_f = rownames(msu) [ msu$n_presences_sps >= 0 ]

	# apply to expression matrices
	msu_f = msu [ ogs_f , ]
	mpr_f = mpr [ ogs_f , ]
	mga_f = mga [ ogs_f , ]
	mlo_f = mlo [ ogs_f , ]

	# apply to conservation matrices
	csu_f = csu [ ogs_f , ]
	cpr_f = cpr [ ogs_f , ]
	cga_f = cga [ ogs_f , ]
	clo_f = clo [ ogs_f , ]

	# reorder nodes according to phylogeny, expression matrices
	mpr_f = mpr_f [ , nod_list ]
	mga_f = mga_f [ , nod_list ]
	mlo_f = mlo_f [ , nod_list ]
	
	# reorder nodes according to phylogeny, conservation matrices
	cpr_f = cpr_f [ , nod_list ]
	cga_f = cga_f [ , nod_list ]
	clo_f = clo_f [ , nod_list ]

	# summarise gains, losses and presences per node
	mss_f = data.frame(row.names = nod_list)
	mss_f$taxa = rownames(mss_f)
	mss_f$gain_exp = colSums(mga_f > 0, na.rm = TRUE)
	mss_f$loss_exp = colSums(mlo_f > 0, na.rm = TRUE)
	mss_f$pres_exp = colSums(mpr_f > 0, na.rm = TRUE)
	mss_f$gain_con = colSums(cga_f > 0, na.rm = TRUE)
	mss_f$loss_con = colSums(clo_f > 0, na.rm = TRUE)
	mss_f$pres_con = colSums(cpr_f > 0, na.rm = TRUE)

	# how many genes among those newly expressed are novel, or not?
	mss_f$gain_exp_and_gain_con     = colSums(mga_f==1 & cga_f==1, na.rm = TRUE)
	mss_f$gain_exp_and_nogain_con   = colSums(mga_f==1 & cga_f==0, na.rm = TRUE)
	# mss_f$nogain_exp_and_gain_con   = colSums(mga_f==0 & cga_f==1, na.rm = TRUE)
	# mss_f$nogain_exp_and_nogain_con = colSums(mga_f==0 & cga_f==0 & mpr_f==1, na.rm = TRUE)

	# how many genes among those not with lost expression are actually absent?
	mss_f$loss_exp_and_abst_con = colSums(mlo_f==1 & cpr_f==0, na.rm = TRUE)

	##### Plots, all #####

	# info to plot
	phy_data = merge(phy_edge, mss_f, by.x = "taxa", by.y = "taxa", all.x = TRUE)
	rownames(phy_data) = phy_data$taxa
	phy_data = phy_data [ order(phy_data$ix_edges), ]
	phy_data = phy_data [ !is.na(phy_data$ix_edges), ]

	# plot aggregated cladogram
	message(sprintf("topology %s | summary tree", topi))
	pdf(sprintf("%s/matrix.%s.tree.pdf", out_fn, topi), height = 10, width = 10)
	layout(mat = matrix(c(1,1), nrow = 2))
	ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1)

	# plot gains and losses as pies
	gain_df = data.frame(
		gain_exp_and_gain_con = phy_data$gain_exp_and_gain_con / phy_data$pres_exp,
		gain_exp_and_nogain_con = phy_data$gain_exp_and_nogain_con / phy_data$pres_exp
		# nogain_exp_and_gain_con = phy_data$nogain_exp_and_gain_con / phy_data$pres_exp,
		# nogain_exp_and_nogain_con = phy_data$nogain_exp_and_nogain_con / phy_data$pres_exp
	)
	gain_df$present_no_gain_exp = 1 - rowSums(gain_df)
	ape::edgelabels(pie = gain_df, piecol = c("chartreuse4","olivedrab2","lightcyan2"), cex = sqrt(phy_data$pres_exp / maxval) * 2)
	ape::edgelabels(pie = c(phy_data$loss_exp_and_abst_con / phy_data$loss_exp), piecol = c("indianred4","indianred1"), cex = sqrt(phy_data$loss_exp / maxval) * 2, adj = -0.3)
	ape::edgelabels(
		text = sprintf("%s\np=%i\n+%i | -%i", phy_data$taxa, phy_data$pres_exp, phy_data$gain_exp, phy_data$loss_exp),
		col = alpha("gray20", 0.9), frame="none", cex = 0.7)
	title(main = sprintf("%s toplogy, neural marker usage", topi))
	title(sub  = sprintf("root expressed p=%i (out of %i conserved)\ntotal genes n=%i", mss_f[root_label,"pres_exp"], mss_f[root_label,"pres_con"], nrow(mga_f)))

	# plot legend
	legend(
		"bottomleft",
		c("new expression & new gene", "new expression & old gene", "already expressed", "lost expression","lost expression & absent gene"),
		fill = c("chartreuse4","olivedrab2","lightcyan2","indianred1","indianred4"),
		bty = "n", cex = 1)
		
	# plot node size legend
	ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1, edge.color = FALSE, show.tip.label = FALSE)
	legend_vec = c(10,20,50,100,200,500,1000)
	ape::tiplabels(pie = c(rep(1, length(legend_vec)), rep(0, length(phy$tip.label) - length(legend_vec))), piecol = c("violet"), cex = sqrt(legend_vec / maxval) * 2, offset = -2)
	ape::tiplabels(c(legend_vec, rep(0, length(phy$tip.label) - length(legend_vec))), col = c("darkorchid"), frame = "none", offset = -3)

	## summary table ##
	
	lss_f = data.frame(
		orthogroup = rownames(msu_f),
		orthogroup_name_Mus = gsub("OG_\\d+:","",ngs_v [ rownames(msu_f) ]),
		presence_expression = msu_f$presence,
		gain_expression = msu_f$gain,
		loss_expression = msu_f$loss,
		presence_conservation = csu_f$presence,
		gain_conservation = csu_f$gain,
		loss_conservation = csu_f$loss,
		n_presences_sps = msu_f$n_presences_sps
	)

	##### Plots, subsets #####
	
	for (seti in c("tfs","ion","sig","neu","siggpcr","myo")) {
		
		message(sprintf("topology %s | summary tree for %s", topi, seti))
		set_t = data.frame()
		for (spi in sps_list) {
			set_fn = sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", seti, spi)
			if (file.exists(set_fn)) {
				set_i = read.table(set_fn, col.names = c("gene","annot"))
				set_t = rbind(set_t, set_i)
			}
		}
		ogs_i = unique(mgs_d [ mgs_d$gene %in% set_t$gene, "orthogroup" ])
		
		# add to summary
		lss_f [ , sprintf("is_%s", seti) ] = lss_f$orthogroup %in% ogs_i
		
		# apply to expression matrices
		msu_i = msu [ ogs_i , ]
		mpr_i = mpr [ ogs_i , ]
		mga_i = mga [ ogs_i , ]
		mlo_i = mlo [ ogs_i , ]

		# apply to conservation matrices
		csu_i = csu [ ogs_i , ]
		cpr_i = cpr [ ogs_i , ]
		cga_i = cga [ ogs_i , ]
		clo_i = clo [ ogs_i , ]

		# reorder nodes according to phylogeny, expression matrices
		mpr_i = mpr_i [ , nod_list ]
		mga_i = mga_i [ , nod_list ]
		mlo_i = mlo_i [ , nod_list ]
		
		# reorder nodes according to phylogeny, conservation matrices
		cpr_i = cpr_i [ , nod_list ]
		cga_i = cga_i [ , nod_list ]
		clo_i = clo_i [ , nod_list ]

		# summarise gains, losses and presences per node
		mss_i = data.frame(row.names = nod_list)
		mss_i$taxa = rownames(mss_i)
		mss_i$gain_exp = colSums(mga_i > 0, na.rm = TRUE)
		mss_i$loss_exp = colSums(mlo_i > 0, na.rm = TRUE)
		mss_i$pres_exp = colSums(mpr_i > 0, na.rm = TRUE)
		mss_i$gain_con = colSums(cga_i > 0, na.rm = TRUE)
		mss_i$loss_con = colSums(clo_i > 0, na.rm = TRUE)
		mss_i$pres_con = colSums(cpr_i > 0, na.rm = TRUE)

		# how many genes among those newly expressed are novel, or not?
		mss_i$gain_exp_and_gain_con     = colSums(mga_i==1 & cga_i==1, na.rm = TRUE)
		mss_i$gain_exp_and_nogain_con   = colSums(mga_i==1 & cga_i==0, na.rm = TRUE)
		# mss_i$nogain_exp_and_gain_con   = colSums(mga_i==0 & cga_i==1, na.rm = TRUE)
		# mss_i$nogain_exp_and_nogain_con = colSums(mga_i==0 & cga_i==0 & mpr_i==1, na.rm = TRUE)

		# how many genes among those not with lost expression are actually absent?
		mss_i$loss_exp_and_abst_con = colSums(mlo_i==1 & cpr_i==0, na.rm = TRUE)

		##### Plots, all #####

		# info to plot
		phy_data = merge(phy_edge, mss_i, by.x = "taxa", by.y = "taxa", all.x = TRUE)
		rownames(phy_data) = phy_data$taxa
		phy_data = phy_data [ order(phy_data$ix_edges), ]
		phy_data = phy_data [ !is.na(phy_data$ix_edges), ]

		# plot aggregated cladogram
		ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1)

		# plot gains and losses as pies
		gain_df = data.frame(
			gain_exp_and_gain_con = phy_data$gain_exp_and_gain_con / phy_data$pres_exp,
			gain_exp_and_nogain_con = phy_data$gain_exp_and_nogain_con / phy_data$pres_exp
			# nogain_exp_and_gain_con = phy_data$nogain_exp_and_gain_con / phy_data$pres_exp,
			# nogain_exp_and_nogain_con = phy_data$nogain_exp_and_nogain_con / phy_data$pres_exp
		)
		gain_df$present_no_gain_exp = 1 - rowSums(gain_df)
		ape::edgelabels(pie = gain_df, piecol = c("chartreuse4","olivedrab2","lightcyan2"), cex = sqrt(phy_data$pres_exp / maxval) * 2)
		ape::edgelabels(pie = c(phy_data$loss_exp_and_abst_con / phy_data$loss_exp), piecol = c("indianred4","indianred1"), cex = sqrt(phy_data$loss_exp / maxval) * 2, adj = -0.2)
		ape::edgelabels(
			text = sprintf("%s\np=%i\n+%i | -%i", phy_data$taxa, phy_data$pres_exp, phy_data$gain_exp, phy_data$loss_exp),
			col = alpha("gray20", 0.9), frame="none", cex = 0.7)
		title(main = sprintf("%s toplogy, neural marker usage, %s", topi, seti))
		title(sub  = sprintf("root expressed p=%i (out of %i conserved)\ntotal genes n=%i", mss_i[root_label,"pres_exp"], mss_i[root_label,"pres_con"], nrow(mga_i)))

		# plot legend
		legend(
			"bottomleft",
			c("new expression & new gene", "new expression & old gene", "already expressed", "lost expression","lost expression & absent gene"),
			fill = c("chartreuse4","olivedrab2","lightcyan2","indianred1","indianred4"),
			bty = "n", cex = 1)
			
	}
	
	dev.off()
	
	# add ortholog member genes
	message("summary | add individual marker genes...")
	for (spi in sps_list) { 

		if (spi %in% c("Dmel","Mlei","Nvec","Mmus","Spis","Hvul")) {
			mar_fn = sprintf("results_panneural_markers_broad/markers.neuron.%s.txt", spi)
		} else if (spi %in% c("Tadh","TrH2","Hhon","HoiH23")) {
			mar_fn = sprintf("results_panneural_markers_broad/markers.peptidergic.%s.txt", spi)
		} else if (spi %in% c("Aque","Spolac")) {
			mar_fn = sprintf("results_panneural_markers_broad/markers.neuroid.%s.txt", spi)
		}
		mar_i = read.table(mar_fn, header = TRUE)
		ogg_l_i = ogg_l [ ogg_l$species == spi, ]		
		ogg_l_i$species = NULL
		ogg_v_i = ogg_l_i$gene
		names(ogg_v_i) = ogg_l_i$OG
		
		lss_f[,sprintf("genes_%s", spi)] = sapply(lss_f$orthogroup, function(v) { 
			genes_in_og = ogg_l_i [ ogg_l_i$OG == v , "gene" ]
			genes_in_og = genes_in_og [ genes_in_og %in% mar_i$gene ]
			genes_in_og = paste(unique(sort(genes_in_og)), collapse = ",")
		})
		
	}

	for (seti in c("tfs","ion","sig","neu","siggpcr","myo")) {
		
		message(sprintf("topology %s | subset plots for %s", topi, seti))
		set_t = data.frame()
		for (spi in sps_list) {
			set_fn = sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", seti, spi)
			if (file.exists(set_fn)) {
				set_i = read.table(set_fn, col.names = c("gene","annot"))
				set_i$species = spi
				set_t = rbind(set_t, set_i)
			}
		}
		ogs_i = unique(mgs_d [ mgs_d$gene %in% set_t$gene, "orthogroup" ])
		
		
		### Plot evolutionary histories ###
		
		# matrix of presence in tips, for genes present in a given ancestor
		ancestor_focus = "Placozoa"
		mpr_f_i = mpr_f   [ ogs_i [ ogs_i %in% rownames(mpr_f) ] , ]
		mpr_f_i = mpr_f_i [ sapply(mpr_f_i[,ancestor_focus], function(v) v> 0) , sps_list ]
		mpr_f_i = mpr_f_i [ order(rowSums(mpr_f_i), decreasing = TRUE) , ]
		
		# matrix of conservation in tips
		cpr_f_i = cpr_f   [ ogs_i [ ogs_i %in% rownames(cpr_f) ] , ]
		cpr_f_i = cpr_f_i [ sapply(cpr_f_i[,ancestor_focus], function(v) v> 0) , sps_list ]
		cpr_f_i = cpr_f_i [ rownames(mpr_f_i) , ]
			
		message(sprintf("topology %s | individual evolutionary histories for %s", topi, seti))
		pdf(sprintf("%s/matrix.%s.histories_set_%s.pdf", out_fn, topi, seti), height = 3, width = 20)
		layout(mat = matrix(1:6, ncol = 6, byrow = TRUE))
		for (ogi in 1:nrow(mpr_f_i)) {
			ogi_name = rownames(mpr_f_i)[ogi]
			ogi_v = factor(mpr_f_i[ogi,], levels = c(0,1))
			cgi_v = factor(cpr_f_i[ogi,], levels = c(0,1))
			ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1)
			title(main = sprintf("%s %s\n%s", ancestor_focus, ogi_name, ngs_v [ ogi_name ]), cex.main = 0.7)
			ape::tiplabels(pch = c(1,19) [ ogi_v ], col = c("chartreuse3"), frame = "none", offset = +2.5, cex = 2)
			ape::tiplabels(pch = c(1,19) [ cgi_v ], col = c("lightcyan4"),  frame = "none", offset = +3.5, cex = 2)
		}
		dev.off()
		
		# add curated annot
		annot_in_og_l = sapply(1:nrow(lss_f), function(o) {
			genes_in_og = c()
			for (spi in c("Tadh","TrH2","Hhon","HoiH23")) {
				genes_in_og = c(genes_in_og, unlist(stringr::str_split(lss_f[o,sprintf("genes_%s", spi)], ",")))
			}
			genes_in_og = genes_in_og [ genes_in_og != ""]
			
			annot_in_og = set_t [ set_t$gene %in% genes_in_og, "annot" ]
			annot_in_og = paste(unique(sort(annot_in_og)), collapse = ",")
			return(annot_in_og)
		} )
		annot_in_og_l [ is.na(annot_in_og_l) ] = ""
		lss_f[sprintf("annot_%s", seti)] = annot_in_og_l


		### Gene expression ###
		
		message(sprintf("topology %s | metacell footprints for %s", topi, seti))
		pdf(sprintf("%s/matrix.%s.expression_set_%s.pdf", out_fn, topi, seti), height = 15, width = 12)
		for (ogi in 1:nrow(mpr_f_i)) {
			ogi_name = rownames(mpr_f_i)[ogi]
			ogi_v = mpr_f_i[ogi,]
			layout(mat = matrix(1:16, nrow = 8, byrow = FALSE))
			for (spi in sps_list) {
				genes_to_plot_d = ogg_l [ ogg_l$OG == ogi_name & ogg_l$species == spi, ]
				genes_to_plot_d = genes_to_plot_d [ genes_to_plot_d$gene %in% rownames(mc_list[[spi]]), ]
				# genes_to_plot_d = mgs_d [ mgs_d$orthogroup == ogi_name & mgs_d$species == spi, ]
				if (nrow(genes_to_plot_d) > 0 & spi %in% names(mc_list)) {
					for (gei in 1:nrow(genes_to_plot_d)) {
						gei_name = genes_to_plot_d$gene[gei]
						mc_v = mc_list[[spi]] [ gei_name, ]
						mc_v [ mc_v <= 1 ] = 1
						barplot(mc_v, col = ct_list[[spi]]$color, log = "y", border = ct_list_col[[spi]], ylim = c(1,10),  ylab = "fp", xlab = "metacells", main = sprintf("%s %s | %s", ngs_v [ ogi_name ], spi, gei_name), space = 0, las = 2)
					}
				} else {
					barplot(0, main = sprintf("%s %s | no data", ngs_v [ ogi_name ], spi), las = 2)
				}
			}
		}
		dev.off()
		
	}
	
	# sort summary table
	message("summary | sort...")
	lss_f$gain_expression = factor(lss_f$gain_expression, levels = c(phy$node.label, phy$tip.label))
	lss_f = lss_f [ order(-lss_f$is_tfs, -lss_f$n_presences_sps, lss_f$gain_expression), ]
	
	# add pfam
	pfam_in_og_l = sapply(1:nrow(lss_f), function(o) {
		genes_in_og = c()
		for (spi in c("Tadh","TrH2","Hhon","HoiH23")) {
			genes_in_og = c(genes_in_og, unlist(stringr::str_split(lss_f[o,sprintf("genes_%s", spi)], ",")))
		}
		genes_in_og = genes_in_og [ genes_in_og != ""]
		annot_in_og = pfam_d_dict [ genes_in_og ]
		annot_in_og = paste(unique(sort(annot_in_og)), collapse = ",")
		return(annot_in_og)
	} )
	pfam_in_og_l [ is.na(pfam_in_og_l) ] = ""
	lss_f[sprintf("pfam_architecture")] = pfam_in_og_l

	
	
	# save
	message("summary | save...")
	write.table(lss_f, sprintf("%s/matrix.%s.tree_sum.csv", out_fn, topi), row.names = FALSE, sep = "\t", quote = FALSE)

}