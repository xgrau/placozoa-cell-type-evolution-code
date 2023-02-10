# libraries
library("igraph")
source("../scripts/helper.R")
source("../scripts/geneSetAnalysis.R")
graphics.off()

#### Input ####

# spslist
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

# paths to input data
inp_fn = "results_gmod_it4/"
ing_fn = "results_gmod_it4_ct_annotations/"
out_fn = "results_gmod_it4_evol_modules/"
dir.create(out_fn, showWarnings = FALSE)
graphics.off()

# pairwise jaccard values
jac_t = read.table(sprintf("%s/module_communities.pairwise_jaccard.csv", ing_fn), header = TRUE)

# load manual annotations
mod = read.table(sprintf("%s/module_communities.curated.csv", ing_fn), header = TRUE)

# apply manual annotations
jac_t_f = data.frame()
for (moi in unique(mod$community)) {
	mois = mod [ mod$community == moi , "module"]
	jac_t_f_i = jac_t [ jac_t$to %in% mois & jac_t$from %in% mois, ]
	jac_t_f = rbind(jac_t_f, jac_t_f_i)
}
jac_g = igraph::graph_from_data_frame(jac_t_f, directed = FALSE)


# annotations per module per species
sps_gmod_ct = data.frame()
pdf(sprintf("%s/pairwise_jaccard.pdf", out_fn), width = 8, height = 8)
for (ni in 1:length(sps_list)) {

	# load gene module cell type lists
	spi = sps_list[ni]
	spi_gmod_ct = read.table(sprintf("%s/top_ct_per_module.global.%s.bct.csv", ing_fn, spi), sep = "\t", header = TRUE)
	spi_gmod_ct$module = paste(spi, spi_gmod_ct$module, sep = ".")
	sps_gmod_ct = rbind(sps_gmod_ct, spi_gmod_ct)

}


# annotations for broad cell types per species
bct_list_t = c()
for (spi in sps_list) {

	# load cell type annotations for it4
	ctt_fn = sprintf("../results_scatlas/results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
	ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")
	bct_list = unique(ctt$broad_cell_type)
	bct_list_t = union(bct_list, bct_list_t)

}


# find communities of modules with a minimum overlap, and group by cell type to which they map
jac_g_components = igraph::components(jac_g)
jac_g_components_d = data.frame("module" = names(jac_g_components$membership), "community" = sprintf("gc%03d",jac_g_components$membership))

# add cell type
gmod_bct_t = data.frame()
gmod_cts_t = data.frame()
for (ni in 1:length(sps_list)) {

	# message
	spi = sps_list[ni]
	gmod_bct_i = read.table(sprintf("%s/top_ct_per_module.global.%s.bct.csv", ing_fn, spi), sep = "\t", header = TRUE)
	gmod_cts_i = read.table(sprintf("%s/top_ct_per_module.global.%s.cts.csv", ing_fn, spi), sep = "\t", header = TRUE)
	gmod_bct_i$module = paste(spi, gmod_bct_i$module, sep = ".")
	gmod_cts_i$module = paste(spi, gmod_cts_i$module, sep = ".")
	gmod_bct_t = rbind(gmod_bct_t, gmod_bct_i)
	gmod_cts_t = rbind(gmod_cts_t, gmod_cts_i)

}
colnames(gmod_bct_t) = c("module","bct")
colnames(gmod_cts_t) = c("module","cts")

# map each gene module to their top tissue
jac_g_components_d = merge(jac_g_components_d, gmod_bct_t, by.x = "module", by.y = "module", all.x = TRUE)
jac_g_components_d = merge(jac_g_components_d, gmod_cts_t, by.x = "module", by.y = "module", all.x = TRUE)

# add species
jac_g_components_d$species = gsub("\\..*","", jac_g_components_d$module)
jac_g_components_d$species = factor(jac_g_components_d$species, levels = sps_list)

# order
jac_g_components_d = jac_g_components_d [ order(jac_g_components_d$community, jac_g_components_d$species) , ]

# annotate comunities
jac_g_components_d_a_lab = aggregate(bct ~ community, data = jac_g_components_d,     function(v) { names(sort(table(v), decreasing = TRUE))[1] })
jac_g_components_d_a_bct = aggregate(bct ~ community, data = jac_g_components_d,     function(v) { paste(unique(v), collapse = ",") })
jac_g_components_d_a_cts = aggregate(cts ~ community, data = jac_g_components_d,     function(v) { paste(unique(v), collapse = ",") })
jac_g_components_d_a_sps = aggregate(species ~ community, data = jac_g_components_d, function(v) { paste(sort(unique(v)), collapse = ",") })
jac_g_components_d_a_nsp = aggregate(species ~ community, data = jac_g_components_d, function(v) { length(unique(v)) })
# aggregate
jac_g_components_d_a = data.frame(community = sort(unique(jac_g_components_d$community)))
rownames(jac_g_components_d_a) = jac_g_components_d_a[,1]
jac_g_components_d_a = merge(jac_g_components_d_a, jac_g_components_d_a_lab, by = "community", all.x = TRUE)
jac_g_components_d_a = merge(jac_g_components_d_a, jac_g_components_d_a_bct, by = "community", all.x = TRUE)
jac_g_components_d_a = merge(jac_g_components_d_a, jac_g_components_d_a_cts, by = "community", all.x = TRUE)
jac_g_components_d_a = merge(jac_g_components_d_a, jac_g_components_d_a_nsp, by = "community", all.x = TRUE)
jac_g_components_d_a = merge(jac_g_components_d_a, jac_g_components_d_a_sps, by = "community", all.x = TRUE)
colnames(jac_g_components_d_a) = c("community","label","bct","cts","num_species","species")

# plot graph
# colors for node borders
sps_gmod_ct$mc_top = NULL
sps_gmod_ct = unique(sps_gmod_ct)
rownames(sps_gmod_ct) = sps_gmod_ct$module
jac_g_b_vec = sps_gmod_ct [ names(V(jac_g)) , "ct_top"  ]
jac_g_b_vec = factor(jac_g_b_vec, levels = bct_list_t)
# colors
colors_for_cts = c("khaki2","gray","darkolivegreen2","chartreuse3","darkorange3","orange","orchid1","royalblue3","darkolivegreen4")
jac_g_b_color = colors_for_cts [ jac_g_b_vec ]
jac_g_b_color [ is.na(jac_g_b_color) ] = "gray50"


# colors for nodes
jac_g_v_color = rep("gray50", length(V(jac_g)))
jac_g_v_color [ names(V(jac_g)) %in% bct_list_t ] = "snow4"
jac_g_v_color [ grepl("^Tadh.",names(V(jac_g)))    ] = "palegreen1"
jac_g_v_color [ grepl("^TrH2.",names(V(jac_g)))    ] = "palegreen3"
jac_g_v_color [ grepl("^Hhon.",names(V(jac_g)))    ] = "hotpink1"
jac_g_v_color [ grepl("^HoiH23.",names(V(jac_g)))  ] = "hotpink3"

# colors for labels
jac_g_v_collb = jac_g_v_color
jac_g_v_collb [ jac_g_v_collb == "snow4" ] = "black"
jac_g_v_collb = scales::alpha(colorspace::darken(jac_g_v_collb, 0.1), 1)

set.seed(1)
jac_g_lay = igraph::layout_with_fr(jac_g, niter = 5000)

pdf(sprintf("%s/module_communities.pairwise_graph.pdf", ing_fn), heigh = 8, width = 8)
igraph::plot.igraph(
	jac_g,
	vertex.size  = 4,
	vertex.color = jac_g_v_collb,
	vertex.frame.color = jac_g_v_collb,
	vertex.label.family = "sans",
	vertex.label.color = jac_g_b_color,
	vertex.label.cex =  0.5,
	edge.arrow.size = 1,
	layout = jac_g_lay,
	edge.color = scales::alpha("snow3", 0.8),
	edge.width = ((E(jac_g)$weight)^(1/1.5)) * 6
)
legend("topleft",    c("Tadh","TrH2","Hhon","HoiH23"), fill = c("palegreen1","palegreen3","hotpink1","hotpink3"), bty = "n", cex = 0.5, title = "Species")
legend_weights = seq(0, 1, length.out = 11)
legend("bottomleft", legend = legend_weights, lwd = (legend_weights ^ (1/1.15)) * 6, bty = "n", cex = 0.5, title = "Jaccard")
legend("bottomright", legend = levels(jac_g_b_vec), col = colors_for_cts, pch = 16, bty = "n", cex = 0.5, title = "CT")
title(sub = "Jaccard overlap between modules (4 sps)\nEdge weight = jaccard")
dev.off()



### Descriptive plots ###


mod$species = gsub("\\..*","",mod$module)
mod$species = factor(mod$species, levels = sps_list)



# open
pdf(sprintf("%s/module_communities.distribution.pdf", ing_fn), height = 6, width = 4)

# num sps per community dist
mod_f = mod [ ,c("community","species") ]
mod_f = mod_f [ !duplicated(mod_f) ,]
xtab_nsps_per_comm = table(table(mod_f$community))
pie(xtab_nsps_per_comm, label = sprintf("%s sps n=%i", names(xtab_nsps_per_comm), xtab_nsps_per_comm), col = c("lightcyan1","lightcyan2","lightcyan3","lightcyan4"))
title(main = "num communities covered by num sps", cex.main = 1)

# num modules per community dist
xtab_nmod_per_comm = table(table(mod$community))
barplot(xtab_nmod_per_comm, horiz = TRUE, las = 1, xlab = "num modules")
title(main = "num modules per community", cex.main = 1)

# sps presence in each community
xtab_sps_per_comm = table(mod_f$species, mod_f$community)
hh = plot_complex_heatmap(
	xtab_sps_per_comm,
	cell_border = gpar(col = "white", lwd = 1, lty = 1),
	heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
	color_max = 1,
	use_raster = FALSE,
	cluster_row = FALSE,
	cluster_col = FALSE
)
draw(hh)

# num communities per broad cell type
mod_f = mod [ ,c("community","classification") ]
mod_f = mod_f [ !duplicated(mod_f) ,]
mod_f$classification = factor(mod_f$classification, levels = c("lipophil","fibre","unknown_1","unknown_2","gland","epithelia","peptidergic","none"))
xtab_ntype = table(mod_f$classification)

barplot(xtab_ntype, las = 1, horiz = TRUE, names.arg = sprintf("%s\nn=%i", names(xtab_ntype), xtab_ntype))
title(main = "num communities per ct specificty", cex.main = 1)

dev.off()






### Functional enrichments in each module and community ###

ing_fn = "results_gmod_it4_evol_communities/"
dir.create(sprintf("%s/enrichment", ing_fn))

# communties, manually defined
mod_dict = mod$community
names(mod_dict) = mod$module

# spslist
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

# for each species
# load annotations
tot_gmod_f = data.frame()
tot_gos_l  = c()
tot_pfm_l  = c()
for (spi in sps_list) {
			
	# message 
	message(sprintf("gmod functional enrichment | %s load", spi))

	# load GOs species i
	ref_gos = gsa_topgo_load_emapper(sprintf("../data/reference/%s_ensembl.GO.csv",spi), index_col_GOs = 2)
	names(ref_gos) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), names(ref_gos), t2g = TRUE)
	
	# load pfam species i
	ref_pfa = gsa_enrichment_load_pfam_list(sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi))
	names(ref_pfa) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), names(ref_pfa), t2g = TRUE)
	
	# concatenate species
	tot_gos_l = c(tot_gos_l, ref_gos)
	tot_pfm_l = c(tot_pfm_l, ref_pfa)
	
	# load gene module gene lists
	spi_gmod = read.table(sprintf("%s/gmod_%s.gmod_annotation.csv", inp_fn, spi), sep = "\t", header = TRUE)
	spi_gmod$gene_module = paste(spi, spi_gmod$gene_module, sep = ".")
	
	# add community
	spi_gmod$community   = mod_dict [ spi_gmod$gene_module ]
	spi_gmod_f = spi_gmod [ !is.na(spi_gmod$community), ]
	tot_gmod_f = rbind(tot_gmod_f, spi_gmod_f)
	
	community_list = sort(unique(spi_gmod$community))
	
	for (comi in community_list) {
		message(sprintf("gmod functional enrichment | %s enrichments %s", spi, comi))
		spi_gmod_i = spi_gmod_f [ spi_gmod_f$community == comi, ]
		enr_gos = gsa_topgo_enrichment(ref_gos, spi_gmod_i$gene, output_prefix = sprintf("%s/enrichment/out.%s", ing_fn, spi), name_fg = comi)
		enr_pfa = gsa_enrichment_hypergeometric(ref_pfa, spi_gmod_i$gene, output_prefix = sprintf("%s/enrichment/out.%s", ing_fn, spi), name_fg = comi)
		
		# joint plot
		if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
			enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
			colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
			enr_pfa_b$ontology = "Pfam"
			enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
			enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
			enr = enr [ enr$pval_test <= 0.1 , ]
			enr = enr [ enr$Significant > 0 , ]
			enr$Term = gsub(" ","_",enr$Term)
			enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
			pdf(sprintf("%s/enrichment/scatter.%s.%s.pdf", ing_fn, spi, comi), width = 16, height = 16)
			layout(matrix(1:4, nrow = 2, byrow = TRUE))
			gsa_enrichment_scatter_plot(enr, main = comi, annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
			dev.off()
		}

	}
			
}


# all species, concatenated
community_list = sort(unique(tot_gmod_f$community))
for (comi in community_list) {
	
	markers_community = tot_gmod_f [ tot_gmod_f$community == comi, "gene"]
	message(sprintf("gmod functional enrichment | %s enrichments %s", "all", comi))
	enr_gos = gsa_topgo_enrichment(tot_gos_l, markers_community, output_prefix = sprintf("%s/enrichment/out.%s", ing_fn, "all"), name_fg = comi)
	enr_pfa = gsa_enrichment_hypergeometric(tot_pfm_l, markers_community, output_prefix = sprintf("%s/enrichment/out.%s", ing_fn, "all"), name_fg = comi)
	
	# joint plot
	if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
		enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
		colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
		enr_pfa_b$ontology = "Pfam"
		enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
		enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
		enr = enr [ enr$pval_test <= 0.1 , ]
		enr = enr [ enr$Significant > 0 , ]
		enr$Term = gsub(" ","_",enr$Term)
		enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
		pdf(sprintf("%s/enrichment/scatter.%s.%s.pdf", ing_fn, "all", comi), width = 16, height = 16)
		layout(matrix(1:4, nrow = 2, byrow = TRUE))
		gsa_enrichment_scatter_plot(enr, main = comi, annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
		dev.off()
	}

}

