#### Input ####

# libraries
source("../scripts/helper.R")

# paths to input data
inp_fn = "results_gmod_it4/"
out_fn = "results_gmod_it4_overlaps/"
dir.create(out_fn, showWarnings = FALSE)

# spslist
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

# load communities
gmc = read.table("results_gmod_it4_ct_annotations/module_communities.curated.csv", header = TRUE)
list_communities = unique(gmc$community)

# output
out_fn = "results_gmod_it4_ct_annotations/"

fp_thr = 1.5

# load cell types
ctt_cts_t = data.frame()
for (spi in sps_list) {
	
	# load cell type info
	ctt_fn = sprintf("../results_scatlas/results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
	ctt = read.table(ctt_fn, header = TRUE, sep = "\t", comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$broad_cell_type = factor(ctt$broad_cell_type, levels = unique(ctt$broad_cell_type))
	ctt_cts = unique(ctt[,c("cell_type","color")])
	ctt_cts$species = spi
	ctt_cts_t = rbind(ctt_cts_t, ctt_cts)
	
}

# orthogroups
og_group_fn = "../data/results_broccoli_ml/dir_step3/orthologous_groups.txt"
ogg = read.table(og_group_fn, sep = "\t", col.names = c("OG","genes"))
ogg_l = data.frame(
	orthogroup = unlist(sapply(1:nrow(ogg), function(i) { v = ogg[i,"genes"] ; vs = unlist(strsplit(v, " ")) ; lv = length(vs) ; vo = rep(ogg[i,"OG"], lv) })),
	gene       = unlist(sapply(1:nrow(ogg), function(i) { v = ogg[i,"genes"] ; vs = unlist(strsplit(v, " ")) }))
)
ogg_l$species = gsub("_.*", "", ogg_l$gene)
ogg_l$species = factor(ogg_l$species, levels = sps_list)
ogg_v = ogg_l$orthogroup
names(ogg_v) = ogg_l$gene

# orthogroups with genes rather than transcripts
ogg_v_t = ogg_v
ogg_v_t_names = names(ogg_v_t)
for (spi in sps_list) {
	ogg_v_t_names [ grepl(sprintf("^%s_",spi), ogg_v_t_names) ] = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), ogg_v_t_names [ grepl(sprintf("^%s_",spi), ogg_v_t_names) ], t2g = TRUE)
}
names(ogg_v_t) = ogg_v_t_names

# all marker genes per community
gmo_l = list()
for (spi in sps_list) {

	# load gene modules
	gmo_i = readRDS(sprintf("results_gmod_it4/gmod_%s.wgcna_object.gmods.rds", spi))
	names(gmo_i) = paste(spi, names(gmo_i), sep = ".")
	gmo_l = c(gmo_l, gmo_i)

}

# get aggregated signal per cell type and community
pdf(sprintf("%s/community_activity.pdf", out_fn), width = 8, height = 16)
gma_d = data.frame()
gma_l = list()
for (spi in sps_list) {

	# load cell type info
	ctt_fn = sprintf("../results_scatlas/results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
	ctt = read.table(ctt_fn, header = TRUE, sep = "\t", comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$broad_cell_type = factor(ctt$broad_cell_type, levels = unique(ctt$broad_cell_type))
	# ctt = ctt [ ctt$broad_cell_type != "trans", ]
	ctt_cts = unique(ctt[,c("cell_type","color")])
	ctt_cts$species = spi
	
	# init database
	metacell::scdb_init("../results_scatlas/data/scdb/",force_reinit=TRUE)
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it4_cts",run_name))
	mc@mc_fp = mc@mc_fp [ , as.character(ctt_cts$cell_type) ]
	
	# mcfps as ogs
	mog_fp = mc@mc_fp
	rownames(mog_fp) = ogg_v_t [ rownames(mog_fp) ]
	mog_fp = mog_fp  [ !is.na(rownames(mog_fp)), ]

	# start species i	
	layout(matrix(1:12, byrow = TRUE, ncol = 2))
	gma_i = matrix(nrow = ncol(mc@mc_fp), ncol = length(list_communities))
	rownames(gma_i) = colnames(mc@mc_fp)
	colnames(gma_i) = list_communities
	for (cmi in list_communities) {
		cmi_name = gmc[ gmc$community == cmi, "label" ][1]
		gmo_in_cmi_v = gmc [ gmc$community == cmi , "module" ]
		gen_in_cmi_v = c()
		if (length(gmo_in_cmi_v) > 0) {
			for (gmi in gmo_in_cmi_v) {
				gen_in_cmi_v = c(gen_in_cmi_v, gmo_l[[gmi]])
			}
			# active genes
			# gen_in_cmi_v = gen_in_cmi_v [ gen_in_cmi_v %in% rownames(mc@mc_fp) ]
			# active OGs
			ogs_in_cmi_v = ogg_v_t [ gen_in_cmi_v ]
			ogs_in_cmi_v = unique(as.character(ogs_in_cmi_v [ !is.na(ogs_in_cmi_v) ]))
			ogs_in_cmi_v = ogs_in_cmi_v [ ogs_in_cmi_v %in% rownames(mog_fp) ]
			# activity with genes
			# act_in_cmi_i = apply(mc@mc_fp [ gen_in_cmi_v,  ], 2, function(c) sum(c >= fp_thr)) / length(gen_in_cmi_v)
			# activity with OGs
			act_in_cmi_i = apply(mog_fp [ ogs_in_cmi_v,  ], 2, function(c) sum(c >= fp_thr)) / length(ogs_in_cmi_v)
			
			# plot activity
			barplot(act_in_cmi_i, las = 2, ylab = "Fraction markers", col = ctt_cts$color, cex.names=  0.6, ylim = c(0,1))
			title(main = sprintf("%s %s\n%s", cmi, cmi_name, spi ), cex.main = 1)
			gma_i [ , cmi ] = act_in_cmi_i
		} else {
			barplot(0,0)
			title(main = sprintf("%s %s\n%s", cmi, cmi_name, spi ))
		}
		
	}
	
	# save (per species and concatenated)
	gma_l[[spi]] = gma_i
	rownames(gma_i) = paste(spi, rownames(gma_i), sep = "|")
	gma_d = rbind(gma_d, gma_i)
	
}

dev.off()


## load tree for cell type order: all ##
cts_phy = ape::read.tree("../results_crosssps/results_trees_cell_types/csps.global.dendrogram.cts.UPGMA.newick")
pep_phy = ape::read.tree("../results_crosssps/results_trees_cell_types/csps.peptidergic.dendrogram.cts.UPGMA.newick")

# reorder for heatmap
gma_d_t_1 = gma_d [ rev(cts_phy$tip.label [ !grepl("peptidergic",cts_phy$tip.label) ]), ]
gma_d_t_2 = gma_d [ rev(pep_phy$tip.label [  grepl("peptidergic",pep_phy$tip.label) ]), ]
gma_d_t   = rbind(gma_d_t_1, gma_d_t_2)
gma_d_t [ is.na(gma_d_t) ] = 0

# annotate columns
dict_communities = gmc [ !duplicated(gmc$community), "label" ]
names(dict_communities) = gmc [ !duplicated(gmc$community), "community" ]
colnames(gma_d_t) = paste(dict_communities [ colnames(gma_d) ], colnames(gma_d_t))

# column order
gmc_order = scp_select_top_markers(t(gma_d_t), matrix_thr = 0, n_top_markers = 1e4, n_markers_rollmean = 2)

# annotate rows
dict_ctt_cts_t = ctt_cts_t$color
names(dict_ctt_cts_t) = paste(ctt_cts_t$species, ctt_cts_t$cell_type, sep = "|")

hm1 = plot_complex_heatmap(
	gma_d_t [,gmc_order],
	name = "frac",
	title_col = "Orthogroup usage across cts and gene modules",
	color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
	cluster_row = FALSE,
	cluster_col = FALSE,
	colors_row = dict_ctt_cts_t,
	color_min = 0.1,
	color_max = 0.8,
	use_raster = FALSE,
	cell_border = gpar(col = "white", lwd = 1, lty = 1),
	heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
)

colnames(gma_d) = paste(dict_communities [ colnames(gma_d) ], colnames(gma_d))
hm2 = plot_complex_heatmap(
	gma_d [,gmc_order],
	name = "frac",
	title_col = "Orthogroup usage across cts and gene modules",
	color_mat = c("gray99","#accbcc","#508490","#004066","#000738"),
	cluster_row = TRUE,
	cluster_col = FALSE,
	colors_row = dict_ctt_cts_t,
	color_min = 0.0,
	color_max = 0.6,
	use_raster = FALSE,
	cell_border = gpar(col = "white", lwd = 1, lty = 1),
	heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
)

pdf(sprintf("%s/community_activity.cts.pdf", out_fn), width = 8, height = 16)
print(hm1)
print(hm2)
dev.off()

message("Done! Tak.")

