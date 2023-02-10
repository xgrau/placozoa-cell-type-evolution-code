# libraries
source("../scripts/helper.R")
source("../scripts/Downstream_functions.R")
library("scales")
graphics.off()
range01 = function(x) { (x-min(x)) / (max(x)-min(x)) }


### Common heatmap panneural ###

tgo_fl = list.files("results_panneural_markers/enrichments/", pattern = "out.*.topgo.csv", full.names = TRUE)
pfa_fl = list.files("results_panneural_markers/enrichments/", pattern = "out.*.hypergeo.csv", full.names = TRUE)
tgo_ids_v  = gsub(".topgo.csv", "", tgo_fl)
pfa_ids_v  = gsub(".hypergeo.csv", "", pfa_fl)
ids_v = intersect(tgo_ids_v, pfa_ids_v)
ids_v = ids_v [ grepl("gain",ids_v) ]

tfu = data.frame()
for (pre in ids_v) {
	
	tgo_fn = sprintf("%s.topgo.csv", pre)
	pfa_fn = sprintf("%s.hypergeo.csv", pre)
	tgo = read.delim(tgo_fn, sep = "\t", header = TRUE)
	pfa = read.table(pfa_fn, sep = "\t", header = TRUE)
	
	tfi = data.frame(
		term = c( apply(tgo[,c("GO.ID","Term")], 1, function(v) paste(v, collapse = " | ")) , pfa$annot ),
		pval = c( tgo$pval_test, pfa$pval ),
		ontology = c(tgo$ontology, rep("Pfam",nrow(pfa))),
		fg = c(tgo$Significant, pfa$freq_in_fg),
		oe = c(tgo$Significant / tgo$Expected, pfa$freq_in_fg / (pfa$freq_in_bg / pfa$total_annot_in_bg) )
	)
	tfi$node = basename(pre)
	
	tfu = rbind(tfu, tfi)
	
}

# filter out pfam
tfu = tfu [ tfu$ontology != "Pfam",]

tfu_p = tfu [ , c("term","node","pval") ]
tfu_n = tfu [ , c("term","node","fg") ]
mfu_p = reshape(tfu_p, idvar = "term", timevar = "node", direction = "wide")
mfu_n = reshape(tfu_n, idvar = "term", timevar = "node", direction = "wide")
rownames(mfu_p) = mfu_p[,1]
rownames(mfu_n) = mfu_n[,1]
mfu_p[,1] = NULL
mfu_n[,1] = NULL
mfu_p = as.matrix(mfu_p)
mfu_n = as.matrix(mfu_n)
mfu_n [ is.na(mfu_n) ] = 0
mfu_p [ is.na(mfu_p) ] = 1


# select markers
# mfu_top = scp_select_top_markers(-log10(mfu_p), matrix_thr = 4, n_top_markers = 30, n_markers_rollmean = 1)

# load must ontologies
mus_top = read.table("../results_scatlas/data/list_ontologies_panneural2.txt", sep = "\t")[,1]
term_dict = rownames(mfu_p)
names(term_dict) = gsub(" \\| .*","", term_dict)
mus_top = term_dict [ mus_top ]
mus_top = mus_top [ mus_top %in% rownames(mfu_p) ]
mus_top = scp_select_top_markers(-log10(mfu_p[mus_top,]), matrix_thr = 1, n_top_markers = 200, n_markers_rollmean = 1)
# mus_top = scp_select_top_markers(-log10(mfu_p[union(mfu_top,mus_top),]), matrix_thr = 1, n_top_markers = 200, n_markers_rollmean = 1)

mfu_pf = mfu_p [ mus_top, ]
mfu_nf = mfu_n [ mus_top, ]

tfu_col = tfu$ontology
names(tfu_col) = tfu$term
tfu_col [ tfu_col == "Pfam" ] = "darkorchid3"
tfu_col [ tfu_col == "BP" ] = "blue3"
tfu_col [ tfu_col == "MF" ] = "darkslategray3"
tfu_col [ tfu_col == "CC" ] = "darkslategray4"


hm = plot_complex_heatmap(
	-log10(mfu_pf),
	color_min = 1,
	color_max = 6,
	color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
	cluster_row = FALSE,
	cluster_col = FALSE,
	colors_row = tfu_col,
	do_dotplot = TRUE,
	dot_size_mat = sqrt(mfu_nf),
	dot_size_min = 1,
	dot_size_max = sqrt(20),
	cex_dotplot = 0.08,
	col_dotplot_border = "black"
)
pdf("results_panneural_markers/enrichments/summary_enrichments_joint.pdf", width = 6, height = 15)
print(hm)
dev.off()





### Common heatmap modules ###

tgo_fl = list.files("results_gmod_it4_evol_communities//enrichment/", pattern = "out.all.*.topgo.csv", full.names = TRUE)
pfa_fl = list.files("results_gmod_it4_evol_communities//enrichment/", pattern = "out.all.*.hypergeo.csv", full.names = TRUE)
tgo_ids_v  = gsub(".topgo.csv", "", tgo_fl)
pfa_ids_v  = gsub(".hypergeo.csv", "", pfa_fl)
ids_v = intersect(tgo_ids_v, pfa_ids_v)
ids_v = ids_v [ !grepl("loss",ids_v) ]

# ids_v = ids_v [ grepl("gc002a", ids_v) | grepl("gc022", ids_v) | grepl("gc024", ids_v) | grepl("gc007b", ids_v) | grepl("gc027", ids_v) ]

# colors per community, ordered by cell type label
gcm = read.table("../results_gene_modules/results_gmod_it4_ct_annotations//module_communities.curated.csv", sep = "\t", header = TRUE)
gcm$classification = factor(gcm$classification, levels = c("lipophil","unknown_1","unknown_2","fibre","gland","epithelia","peptidergic","none"))
ctt_fn = sprintf("../results_scatlas/results_metacell_it4/annotation_mc.TrH2.it4.reordered.tsv")
ctt = read.table(ctt_fn, header = TRUE)
gcm$label = factor(gcm$label, levels = unique(ctt$cell_type))
gcm = gcm [ order(gcm$classification, gcm$label, gcm$community), ]
gcm$community = factor(gcm$community, levels = unique(gcm$community))

# colors
com_com_u = unique(gcm [ , c("community","color") ])
cols_cts = com_com_u$color
names(cols_cts) = com_com_u$community



tfu = data.frame()
for (pre in ids_v) {
	
	tgo_fn = sprintf("%s.topgo.csv", pre)
	pfa_fn = sprintf("%s.hypergeo.csv", pre)
	tgo = read.delim(tgo_fn, sep = "\t", header = TRUE)
	pfa = read.table(pfa_fn, sep = "\t", header = TRUE)
	
	tfi = data.frame(
		term = c( apply(tgo[,c("GO.ID","Term")], 1, function(v) paste(v, collapse = " | ")) , pfa$annot ),
		pval = c( tgo$pval_test, pfa$pval ),
		ontology = c(tgo$ontology, rep("Pfam",nrow(pfa))),
		fg = c(tgo$Significant, pfa$freq_in_fg),
		oe = c(tgo$Significant / tgo$Expected, pfa$freq_in_fg / (pfa$freq_in_bg / pfa$total_annot_in_bg) )
	)
	tfi$node = gsub("out\\.all\\.","",basename(pre))
	
	tfu = rbind(tfu, tfi)
	
}

# filter out pfam
tfu = tfu [ tfu$ontology != "Pfam",]

# matrices
tfu_p = tfu [ , c("term","node","pval") ]
tfu_n = tfu [ , c("term","node","fg") ]
mfu_p = reshape(tfu_p, idvar = "term", timevar = "node", direction = "wide")
mfu_n = reshape(tfu_n, idvar = "term", timevar = "node", direction = "wide")
rownames(mfu_p) = mfu_p[,1]
rownames(mfu_n) = mfu_n[,1]
mfu_p[,1] = NULL
mfu_n[,1] = NULL
mfu_p = as.matrix(mfu_p)
mfu_n = as.matrix(mfu_n)
mfu_n [ is.na(mfu_n) ] = 0
mfu_p [ is.na(mfu_p) ] = 1

# order based on community classification?
colnames(mfu_p) = gsub("^pval\\.", "", colnames(mfu_p))
colnames(mfu_n) = gsub("^fg\\.", "", colnames(mfu_n))
mfu_p = mfu_p [ , levels(gcm$community) ]
mfu_n = mfu_n [ , levels(gcm$community) ]

# top markers
mfu_top = scp_select_top_markers(-log10(mfu_p), matrix_thr = 4, n_top_markers = 10, n_markers_rollmean = 1)

# load must ontologies
mus_top = read.table("../results_scatlas/data/list_ontologies_modules.txt", sep = "\t")[,1]
term_dict = rownames(mfu_p)
names(term_dict) = gsub(" \\| .*","", term_dict)
mus_top = term_dict [ mus_top ]
mus_top = mus_top [ mus_top %in% rownames(mfu_p) ]
# mus_top = scp_select_top_markers(-log10(mfu_p[mus_top,]), matrix_thr = 1, n_top_markers = 200, n_markers_rollmean = 1)
mus_top = scp_select_top_markers(-log10(mfu_p[union(mfu_top,mus_top),]), matrix_thr = 1, n_top_markers = 200, n_markers_rollmean = 1)

# filter 
mfu_pf = mfu_p [ mus_top, ]
mfu_nf = mfu_n [ mus_top, ]


tfu_col = tfu$ontology
names(tfu_col) = tfu$term
tfu_col [ tfu_col == "Pfam" ] = "darkorchid3"
tfu_col [ tfu_col == "BP" ] = "blue3"
tfu_col [ tfu_col == "MF" ] = "darkslategray3"
tfu_col [ tfu_col == "CC" ] = "darkslategray4"


hm = plot_complex_heatmap(
	-log10(mfu_pf),
	color_min = 1,
	color_max = 12,
	color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
	cluster_row = FALSE,
	cluster_col = FALSE,
	colors_row = tfu_col,
	colors_col = cols_cts,
	do_dotplot = TRUE,
	dot_size_mat = sqrt(mfu_nf),
	dot_size_min = sqrt(2),
	dot_size_max = sqrt(20),
	cex_dotplot = 0.007,
	col_dotplot_border = "black"
)
pdf("results_gmod_it4_evol_communities//enrichment/summary_enrichments_joint.pdf", width = 12, height = 40)
print(hm)
dev.off()





### Panneural per-node plots ###

min_pval = 1e-12

# individual plots
pdf("results_panneural_markers/enrichments/summary_enrichments.pdf", width = 6, height = 24)
# layout(matrix(1:5, byrow=TRUE,ncol=1))
for (pre in ids_v) {
	
	tgo_fn = sprintf("%s.topgo.csv", pre)
	pfa_fn = sprintf("%s.hypergeo.csv", pre)
	tgo = read.delim(tgo_fn, sep = "\t", header = TRUE)
	pfa = read.table(pfa_fn, sep = "\t", header = TRUE)
	
	tfu = data.frame(
		term = c( apply(tgo[,c("GO.ID","Term")], 1, function(v) paste(v, collapse = " | ")) , pfa$annot ),
		pval = c( tgo$pval_test, pfa$pval ),
		ontology = c(tgo$ontology, rep("Pfam",nrow(pfa))),
		fg = c(tgo$Significant, pfa$freq_in_fg),
		oe = c(tgo$Significant / tgo$Expected, pfa$freq_in_fg / (pfa$freq_in_bg / pfa$total_annot_in_bg) )
	)
	print(pre)
	
	# filter
	tfu$pval [ is.na(tfu$pval) ] = min(tfu$pval [ !is.na(tfu$pval) ])
	tfu$ontology = factor(tfu$ontology, levels = c("MF","BP","CC","Pfam"))
	tfu = tfu [ tfu$fg >= 1 & tfu$pval <= 1e-2, ]
	tfu = tfu [ order(tfu$pval), ]
	
	tfu$symbol = 20
	tfu$symbol [ tfu$pval < min_pval ] = 17
	tfu$pval [ tfu$pval < min_pval ] = min_pval

	# which to show?	
	ix_show = which(tfu$pval <= 1e-3 & tfu$oe >= 2)
	if (length(ix_show) == 0) {
		ix_show = 1:nrow(tfu)
	}
	
	# # # plot
	# plot(tfu$oe, tfu$pval, log = "xy", col = c("cadetblue3","dodgerblue3","cadetblue4","darkorchid3")[ tfu$ontology ], pch = tfu$symbol, ylim = c(1, min_pval), las = 2, ylab = "p-value", xlab = "Observed/expected")
	# text(tfu$term[ix_show], x = tfu$oe[ix_show], y = tfu$pval[ix_show], col = scales::alpha("thistle4",0.7), cex = 0.5, pos = 4)
	# title(main = basename(pre))

	tfu$size = log10(tfu$oe)
	tfu$size [ tfu$size > 5 ] = 5
	tfu$size = sqrt(tfu$size)
	# plot(tfu$fg, tfu$pval, log = "xy", col = c("cadetblue3","dodgerblue3","cadetblue4","darkorchid3")[ tfu$ontology ], pch = tfu$symbol, ylim = c(1e-2, min_pval), las = 2, ylab = "p-value", xlab = "Observed", cex = tfu$size)
	# text(tfu$term[ix_show], x = tfu$fg[ix_show], y = tfu$pval[ix_show], col = scales::alpha("thistle3",0.7), cex = 0.5, pos = 4)
	# title(main = basename(pre))
	write.table(tfu, sprintf("%s.enrichments_table.csv",pre), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
	
	tfu_f = tfu [ tfu$pval < 1e-3 , ]
	plot(y=1:nrow(tfu_f) , x = tfu_f$pval, xlim = c(1e-3, min_pval), cex = tfu_f$size, pch = 20, col = c("cadetblue3","dodgerblue3","cadetblue4","darkorchid3")[ tfu_f$ontology ], log = "x", las = 2, ylab = "", xlab = "p-value", yaxt = "n")
	text(tfu_f$term, x = tfu_f$pval, y = 1:nrow(tfu_f), col = scales::alpha(c("cadetblue3","dodgerblue3","cadetblue4","darkorchid3")[ tfu_f$ontology ],0.7), cex = 0.5, pos = 4)
	title(main = basename(pre))
	
}
dev.off()



### Per module plots ###

tgo_fl = list.files("results_gmod_it4_evol_communities//enrichment/", pattern = "out.all.*.topgo.csv", full.names = TRUE)
pfa_fl = list.files("results_gmod_it4_evol_communities//enrichment/", pattern = "out.all.*.hypergeo.csv", full.names = TRUE)
tgo_ids_v  = gsub(".topgo.csv", "", tgo_fl)
pfa_ids_v  = gsub(".hypergeo.csv", "", pfa_fl)
ids_v = intersect(tgo_ids_v, pfa_ids_v)
ids_v = ids_v [ !grepl("loss",ids_v) ]

min_pval = 1e-30
max_pval = 1e-4

pdf("results_gmod_it4_evol_communities//enrichment/summary_enrichments.pdf", width = 6, height = 24)
# layout(matrix(1:5, byrow=TRUE,ncol=1))
for (pre in ids_v) {
	
	tgo_fn = sprintf("%s.topgo.csv", pre)
	pfa_fn = sprintf("%s.hypergeo.csv", pre)
	tgo = read.delim(tgo_fn, sep = "\t", header = TRUE)
	pfa = read.table(pfa_fn, sep = "\t", header = TRUE)
	
	tfu = data.frame(
		term = c( apply(tgo[,c("GO.ID","Term")], 1, function(v) paste(v, collapse = " | ")) , pfa$annot ),
		pval = c( tgo$pval_test, pfa$pval ),
		ontology = c(tgo$ontology, rep("Pfam",nrow(pfa))),
		fg = c(tgo$Significant, pfa$freq_in_fg),
		oe = c(tgo$Significant / tgo$Expected, pfa$freq_in_fg / (pfa$freq_in_bg / pfa$total_annot_in_bg) )
	)
	print(pre)
	
	# filter
	tfu$pval [ is.na(tfu$pval) ] = min(tfu$pval [ !is.na(tfu$pval) ])
	tfu$ontology = factor(tfu$ontology, levels = c("MF","BP","CC","Pfam"))
	tfu = tfu [ tfu$fg >= 1 & tfu$pval <= 1e-2, ]
	tfu = tfu [ order(tfu$pval), ]
	
	tfu$symbol = 20
	tfu$symbol [ tfu$pval < min_pval ] = 17
	tfu$pval [ tfu$pval < min_pval ] = min_pval

	# which to show?	
	ix_show = which(tfu$pval <= 1e-3 & tfu$oe >= 2)
	if (length(ix_show) == 0) {
		ix_show = 1:nrow(tfu)
	}
	
	tfu$size = log10(tfu$oe)
	tfu$size [ tfu$size > 5 ] = 5
	tfu$size = sqrt(tfu$size)
	# plot(tfu$fg, tfu$pval, log = "xy", col = c("cadetblue3","dodgerblue3","cadetblue4","darkorchid3")[ tfu$ontology ], pch = tfu$symbol, ylim = c(1e-2, min_pval), las = 2, ylab = "p-value", xlab = "Observed", cex = tfu$size)
	# text(tfu$term[ix_show], x = tfu$fg[ix_show], y = tfu$pval[ix_show], col = scales::alpha("thistle3",0.7), cex = 0.5, pos = 4)
	# title(main = basename(pre))
	write.table(tfu, sprintf("%s.enrichments_table.csv",pre), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
	
	tfu_f = tfu [ tfu$pval < max_pval , ]
	plot(y=1:nrow(tfu_f) , x = tfu_f$pval, xlim = c(max_pval, min_pval), cex = tfu_f$size, pch = 20, col = c("cadetblue3","dodgerblue3","cadetblue4","darkorchid3")[ tfu_f$ontology ], log = "x", las = 2, ylab = "", xlab = "p-value",  yaxt = "n")
	text(tfu_f$term, x = tfu_f$pval, y = 1:nrow(tfu_f), col = scales::alpha(c("cadetblue3","dodgerblue3","cadetblue4","darkorchid3")[ tfu_f$ontology ],0.7), cex = 0.5, pos = 4)
	title(main = basename(pre))
	
}
dev.off()