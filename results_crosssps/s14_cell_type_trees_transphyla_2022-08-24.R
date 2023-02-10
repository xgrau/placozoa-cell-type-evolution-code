#### Input ####

# libraries
source("../scripts/helper.R")
library("umap")
library("phangorn")
library("ape")

# input
heatmap_colors = c("white","orange","orangered2","#520c52")
out_fn = "results_ct_trees/"
out_fn = "results_trees_cell_types_transphyla/"
dir.create(out_fn, showWarnings = FALSE)

# species
list_species = c("Tadh","TrH2","Hhon","HoiH23","Nvec","Spis","Hvul")
# list_species = c("Tadh","TrH2","Hhon","HoiH23","Nvec","Spis")
list_outgroups = c("Nvec","Spis","Hvul")
spr = "TrH2"
spq = list_species [ !list_species %in% spr ]

# completeness
fraction_genes = 0.85

# fc threshold
fc_thr = 1.5

# fraction of variance required in PCA based dendrogram
fraction_required = 0.4

# bootstrap
num_bs = 1000
collapse_bs = num_bs * 0.3


# analyses
set_list = list(
	ann_fn = "../results_scatlas/results_metacell_it4/",
	icc_fn = "results_alignment_icc/",
	ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv",
	focus_list = list("bct" = "broad_cell_type"))

ignore_cell_types = c("unknown","calicoblast","cnidocyte","precursors","mitotic_host_cells","unannotated","Hemocyte","Female_Reproductive_System","Male_Reproductive_System","Gland")

set_id = "all"
ann_fn = set_list$ann_fn
icc_fn = set_list$icc_fn
ctt_sprtinf_string = set_list$ctt_sprtinf_string
focus_list = set_list$focus_list

for (i in 1:length(focus_list)) {
	
	focid = names(focus_list)[[i]]
	focus = focus_list[[i]]

	# load cell type data for all species
	ann_cts = c()
	for (spi in list_species) {
		if (spi %in% c(list_outgroups)) {
			ctt_spi_fn = sprintf("../results_scatlas/data/scdb_outgroups/annot.%s.tsv", spi)
		} else {
			ctt_spi_fn = sprintf(ctt_sprtinf_string, ann_fn, spi)
		}
		ctt_spi = read.table(ctt_spi_fn, header = TRUE, comment.char = "", sep = "\t")
		# get color annotations for each species
		ann_spi = unique(ctt_spi[,c(focus,"color")])
		ann_spi = ann_spi [ !duplicated(ann_spi[,focus]), ]
		ann_spi_v = ann_spi[,2]
		names(ann_spi_v) = paste(spi, ann_spi[,1], sep = "|")
		# concatenate annotations
		ann_cts = c(ann_cts, ann_spi_v)
	}
	
	# loop query species
	for (spi in spq) {

		# load cross-species marker data (from `s01`)
		message(sprintf("csps %s | load %s %s-%s csps data...", set_id, focid, spr, spi))
		if (file.exists(sprintf("%s/csps_icc.%s.%s-%s.rds", icc_fn, focid, spr, spi))) {
			csps_i = readRDS(sprintf("%s/csps_icc.%s.%s-%s.rds", icc_fn, focid, spr, spi))
		} else {
			csps_i = readRDS(sprintf("%s/csps_icc.%s.%s-%s.rds", icc_fn, focid, spi, spr))
			csps_i$sp1n = csps_i$sp2
			csps_i$sp2n = csps_i$sp1
			csps_i$sp1 = csps_i$sp1n
			csps_i$sp2 = csps_i$sp2n
		}
		
		if (spi %in% list_outgroups) {
			csps_i$sp2 = csps_i$sp2 [ , ! colnames(csps_i$sp2) %in% ignore_cell_types ]
		}
		
		if (spi == spq[1]) {
			csps_i_m = cbind(csps_i$sp1, csps_i$sp2)
			colnames(csps_i_m) = c(paste(spr, colnames(csps_i$sp1), sep = "|"), paste(spi, colnames(csps_i$sp2), sep = "|"))
			csps_m = csps_i_m
		} else {
			csps_i_m = csps_i$sp2
			rownames(csps_i_m) = rownames(csps_i$sp1)
			colnames(csps_i_m) = paste(spi, colnames(csps_i$sp2), sep = "|")
			# merge and add NA to missing
			csps_m = merge(csps_m, csps_i_m, by = "row.names", all = TRUE)
			rownames(csps_m) = csps_m$Row.names
			csps_m$Row.names = NULL
			csps_m = as.matrix(csps_m)
			# merge only unique matches
			# csps_m = csps_m [ intersect(rownames(csps_m), rownames(csps_i_m)), ]
			# csps_m = cbind(csps_m, csps_i_m [ intersect(rownames(csps_m), rownames(csps_i_m)), ])
		}

	}

	# drop trans
	csps_m = csps_m [ ,!grepl("trans",colnames(csps_m)) ]
	
	# remove markers with an excess of absences (keep genes present in 50% of tissues)
	csps_m = csps_m [ rowSums(!is.na(csps_m)) >= ncol(csps_m) * fraction_genes , ]
	
	# na to zero
	csps_m [ is.na(csps_m) ] = 0
	
	# quantile normalisation
	csps_m_n = quantile_normalisation(csps_m)
	colnames(csps_m_n) = colnames(csps_m)
	csps_m = csps_m_n

	# footprint discretisation
	message(sprintf("csps %s | %s discretise fps at fc>=%.2f...", set_id, focid, fc_thr))

	# binarise
	csps_m_b = (csps_m >= fc_thr) * 1

	message(sprintf("csps %s | %s n=%i markers in %i %ss...", set_id, focid, nrow(csps_m_b), ncol(csps_m_b), focus))

	
	# reduce
	pdf(sprintf("%s/csps.%s.dimred.%s.pdf", out_fn, set_id, focid))
	# umap
	message(sprintf("csps %s | %s reduce UMAP...", set_id, focid))
	csps_u = umap::umap(t(csps_m_b))
	plot(csps_u$layout, col = ann_cts [ rownames(csps_u$layout) ], pch = c(15,17,18,19) [ factor(gsub("\\|.*","", rownames(csps_u$layout))) ] )
	title(main = "umap")
	text(csps_u$layout, rownames(csps_u$layout), col = scales::alpha(ann_cts [ rownames(csps_u$layout) ], 0.2), cex = 0.6)

	# pca
	message(sprintf("csps %s | %s reduce PCA...", set_id, focid))
	if (ncol(csps_m_b) < nrow(csps_m_b)) {
	
		csps_p = princomp(csps_m_b)
		plot(csps_p$loadings[,c(1,2)], col = ann_cts [ rownames(csps_p$loadings) ], pch = c(15,17,18,19) [ factor(gsub("\\|.*","", rownames(csps_p$loadings))) ] )
		title(main = "pca")
		text(csps_p$loadings[,c(1,2)], rownames(csps_p$loadings), col = scales::alpha(ann_cts [ rownames(csps_p$loadings) ], 0.2), cex = 0.6)

		# top?
		pov = csps_p$sdev^2/sum(csps_p$sdev^2)
		n_pcs = which(cumsum(pov) >= fraction_required)[1]
	
		plot(pov, col = "blue", main = "fraction variance", sub = sprintf(">=%.1fpp variance explained at PC %i", fraction_required * 100, n_pcs))

	}
	
	dev.off()


	# plot height
	plot_height = ceiling(ncol(csps_m_b) / 6 + 4)

	# dendrogram based on first X PCs
	if (ncol(csps_m_b) < nrow(csps_m_b)) {
		
		message(sprintf("csps %s | %s reduce PCA, keep n=%i PCs...", set_id, focid, n_pcs))
		csps_p_f = csps_p$loadings[,1:n_pcs]

		message(sprintf("csps %s | %s plot PCA dendrogram...", set_id, focid))
		csps_p_f_h = hclust(dist(csps_p_f, method = "manhattan"), method = "average")
		csps_p_f_t = as.phylo(csps_p_f_h)

		pdf(sprintf("%s/csps.%s.dendrogram.%s.PCA.pdf", out_fn, set_id, focid), height = plot_height, width = 12)
		ape::plot.phylo(ladderize(csps_p_f_t), tip.color = ann_cts [ csps_p_f_t$tip.label ], font = 1, underscore = TRUE, type = "phylo", main = sprintf("upgma from %i PCs", n_pcs))
		ape::add.scale.bar()
		dev.off()
		ape::write.tree(ladderize(csps_p_f_t), sprintf("%s/csps.%s.dendrogram.%s.PCA.newick", out_fn, set_id, focid))
		
	}

	# binary matrix
	ali_f = as.phyDat(t(csps_m_b), type = "USER", levels = names(table(csps_m_b)), names = rownames(csps_m_b), return.index = TRUE)
	
	# distance and initial parsimony tree
	message(sprintf("csps %s | %s parsimony tree...", set_id, focid))
	parsimony_fun =  function(x) { phangorn::upgma(phangorn::dist.logDet(x)) }
	ali_f_phy = parsimony_fun(ali_f)

	# bootstrap in UPGMA (may fail, retry often)
	boo_fail = 0
	boo_success = FALSE
	while (boo_fail < num_bs && !boo_success) {
		ali_f_boo = tryCatch( phangorn::bootstrap.phyDat(ali_f, parsimony_fun, multicore = TRUE, mc.cores = 30, bs = num_bs), error = function(e) e)
		if (!inherits(ali_f_boo, "simpleError")) {
			boo_success = TRUE
			message(sprintf("bootstrap try: %ith iteration succeeds", boo_fail))
		} else {
			message(sprintf("bootstrap try: %ith iteration fails", boo_fail))
			boo_fail = boo_fail + 1
		}
	}
	
	
	# collapse poorly supported nodes
	collapse_nodes = which(prop.clades(ali_f_phy, ali_f_boo) < collapse_bs) + length(ali_f_phy$tip.label)
	rownames(ali_f_phy$edge) = as.character(1:nrow(ali_f_phy$edge))
	collapse_edges = sort(as.numeric(rownames(ali_f_phy$edge) [ ali_f_phy$edge[,2] %in% collapse_nodes ]))
	ali_f_phy_col = ali_f_phy
	ali_f_phy_col$edge.length[ collapse_edges ] = 0
	ali_f_phy_col = ape::di2multi(ali_f_phy_col)
	
	# plot
	pdf(sprintf("%s/csps.%s.dendrogram.%s.UPGMA.pdf", out_fn, set_id, focid), height = plot_height, width = 12)
	phangorn::plotBS(ali_f_phy, ali_f_boo, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP")
	ape::add.scale.bar()
	phangorn::plotBS(ali_f_phy_col, ali_f_boo, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP", use.edge.length = FALSE)
	phangorn::plotBS(ali_f_phy, ali_f_boo, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with TBE", fc_thr), root.edge = TRUE, method = "TBE", digits = 1)
	ape::add.scale.bar()
	ape::plot.phylo(ali_f_phy, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, underscore = TRUE, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f", fc_thr), root.edge = TRUE, use.edge.length = TRUE)
	# ape::edgelabels()
	# ape::nodelabels()
	ape::add.scale.bar()
	dev.off()
	ape::write.tree(ali_f_phy, sprintf("%s/csps.%s.dendrogram.%s.UPGMA.newick", out_fn, set_id, focid))
	
	# save matrix
	saveRDS(csps_m, sprintf("%s/csps.%s.cspsmatrix.%s.rds", out_fn, set_id, focid))
	
	message(sprintf("csps %s | %s trees done!", set_id, focid))
	
}
