# libraries
library("alluvial")
library("igraph")
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/geneSetAnalysis.R"))

# paths to input data
inp_fn = "results_gmod_it4/"
ing_fn = "results_gmod_it4_ct_annotations/"
out_fn = "results_gmod_it4_ct_annotations/"
dir.create(out_fn, showWarnings = FALSE)

# spslist
sps_list = c("Tadh","TrH2","Hhon","HoiH23")
list_focus_cts = c("lipophil", "peptidergic", "fibre", "epithelia", "gland", "unknown_1", "unknown_2")

# orthogroups
og_groups_fn = "../data/results_broccoli_ml/dir_step3/orthologous_groups.txt"

# load orthogroups
ogg_w = read.table(og_groups_fn, sep = "\t", col.names = c("orthogroup","genes"))
ogg_d = data.frame(
	orthogroup = unlist( sapply( 1:length(ogg_w$orthogroup), function(n) rep(ogg_w$orthogroup[n], length(unlist(strsplit(ogg_w$genes[n], split = " ")))) ) ),
	transcript = unlist( sapply( ogg_w$genes,                function(v) strsplit(v, split = " ")) )
)
rownames(ogg_d) = NULL
ogg_d$species = gsub("_.*","",ogg_d$transcript)
ogg_d = ogg_d [ ogg_d$species %in% sps_list, ]

# orthogroups transcript to gene
ogg_d$gene = NA
for (spi in sps_list) {
	ogg_d$gene [ ogg_d$species == spi ] = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), ogg_d$transcript [ ogg_d$species == spi ] , t2g = TRUE)
}

# vector
ogg_v = ogg_d$orthogroup
names(ogg_v) = ogg_d$gene



#### Map gene modules to orthogroups ####

# loop
for (ni in 1:length(sps_list)) {

	# message
	spi = sps_list[ni]
	message(sprintf("gmod overlap | %s: load", spi))

	# load gene module gene lists
	spi_gmod = read.table(sprintf("%s/gmod_%s.gmod_annotation.csv", inp_fn, spi), sep = "\t", header = TRUE)
	# load gene module cell type lists
	spi_gmod_ct = read.table(sprintf("%s/top_ct_per_module.global.%s.bct.csv", ing_fn, spi), sep = "\t", header = TRUE)

	# map each gene module to their top tissue
	spi_gmod_m = merge(spi_gmod, spi_gmod_ct, by.x = "gene_module", by.y = "module")

	# gene to orthogroup mapping
	spi_gmod_m$orthogroup = ogg_v [ spi_gmod_m$gene ]

	# load tfs
	spi_tfs = read.table(sprintf("../data/gene_annotations/tfs.%s_genes.curated.csv", spi), col.names = c("gene","annotation"))
	spi_tfs$gene = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), spi_tfs$gene, t2g = TRUE)
	spi_tfs_v = spi_tfs$annotation
	names(spi_tfs_v) = spi_tfs$gene
	# add tf names to table
	spi_gmod_m$is_transcription_factor = spi_gmod_m$gene %in% spi_tfs$gene
	spi_gmod_m$gene_name [ spi_gmod_m$is_transcription_factor ] = spi_tfs_v [ spi_gmod_m$gene ] [ spi_gmod_m$is_transcription_factor ]

	# remove incomparable rows
	spi_gmod_m = spi_gmod_m [ !is.na(spi_gmod_m$orthogroup), ]

	# add sps name to gene module
	spi_gmod_m$gene_module = paste(spi, spi_gmod_m$gene_module, sep = ".")

	# save
	write.table(spi_gmod_m, sprintf("%s/gmod_%s.orthogroups.csv", out_fn, spi), sep = "\t", quote = FALSE, row.names = FALSE)

	if (ni == 1) {
		sps_gmod_m = spi_gmod_m [ !is.na(spi_gmod_m$orthogroup), c("orthogroup","gene_module","ct_top") ]
		colnames(sps_gmod_m)[2] = paste("module", spi, sep = "_")
		colnames(sps_gmod_m)[3] = paste("ct",     spi, sep = "_")
	} else {
		sps_gmod_m = merge(sps_gmod_m, spi_gmod_m [ !is.na(spi_gmod_m$orthogroup), c("orthogroup","gene_module","ct_top") ], by = "orthogroup", all.x = TRUE, all.y = TRUE)
		colnames(sps_gmod_m)[ ncol(sps_gmod_m) - 1 ] = paste("module", spi, sep = "_")
		colnames(sps_gmod_m)[ ncol(sps_gmod_m)     ] = paste("ct",     spi, sep = "_")
	}

}



# align gene modules across species based on shared orthogroups
sps_gmod_ct = data.frame()
pdf(sprintf("%s/pairwise_jaccard.pdf", out_fn), width = 8, height = 8)
for (ni in 1:length(sps_list)) {

	# load gene module cell type lists
	spi = sps_list[ni]
	spi_gmod_ct = read.table(sprintf("%s/top_ct_per_module.global.%s.bct.csv", ing_fn, spi), sep = "\t", header = TRUE)
	spi_gmod_ct$module = paste(spi, spi_gmod_ct$module, sep = ".")
	sps_gmod_ct = rbind(sps_gmod_ct, spi_gmod_ct)

	for (nj in 1:length(sps_list)) {
		if (ni < nj) {

			# load
			spi = sps_list[ni]
			spj = sps_list[nj]
			spi_gmod_m = read.table(sprintf("%s/gmod_%s.orthogroups.csv", out_fn, spi), sep = "\t", header = TRUE)
			spj_gmod_m = read.table(sprintf("%s/gmod_%s.orthogroups.csv", out_fn, spj), sep = "\t", header = TRUE)

			# compare module per module
			# create matrix of og activity per module
			spi_modog = spi_gmod_m [ , c("gene_module","orthogroup") ]
			spj_modog = spj_gmod_m [ , c("gene_module","orthogroup") ]
			spi_modog$present = 1
			spj_modog$present = 1
			spi_modog_m = reshape(spi_modog, idvar="orthogroup", timevar="gene_module", direction="wide")
			spj_modog_m = reshape(spj_modog, idvar="orthogroup", timevar="gene_module", direction="wide")
			spj_modog_m [ is.na(spj_modog_m) ] = 0
			spi_modog_m [ is.na(spi_modog_m) ] = 0
			rownames(spi_modog_m) = spi_modog_m$orthogroup
			rownames(spj_modog_m) = spj_modog_m$orthogroup
			spi_modog_m$orthogroup = NULL
			spj_modog_m$orthogroup = NULL
			colnames(spi_modog_m) = gsub("^present\\.","",colnames(spi_modog_m))
			colnames(spj_modog_m) = gsub("^present\\.","",colnames(spj_modog_m))

			# shared ogs
			shared_ogs = intersect(rownames(spi_modog_m), rownames(spj_modog_m))

			# jaccard similarity
			jac_ij = jaccard(spi_modog_m[shared_ogs,], spj_modog_m[shared_ogs,])

			# plot pairwise
			row_order = unique(apply(jac_ij, 2, function(r) names(which.max(r)) ))
			row_order = c(row_order, setdiff(rownames(jac_ij), row_order))
			jac_ij_o = jac_ij [ row_order , ]
			hm = plot_complex_heatmap(
				jac_ij_o,
				color_mat = c("gray98","#d6e72e","#6fb600","#003f4d"),
				color_min = 0,
				color_max = 0.6,
				use_raster = FALSE,
				cluster_row = FALSE,
				cluster_col = FALSE,
				cell_border = gpar(col = "white", lwd = 1, lty = 1),
				heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
			)
			print(hm)

			thr_jac = 0.1
			jac_v_ij = apply(jac_ij, 1, function(r) names(which(r>thr_jac)))
			jac_d_ij = data.frame(
				to = unlist(sapply(1:length(jac_v_ij), function(e) rep(names(jac_v_ij)[e], lengths(jac_v_ij)[e]))),
				from = unlist(jac_v_ij))
			rownames(jac_d_ij) = NULL

			# add edge weights
			jac_d_ij$weight = sapply(1:nrow(jac_d_ij), function(i) {
				r = jac_d_ij[i,1]
				c = jac_d_ij[i,2]
				jac_ij[r,c]
			} )

			if (ni == 1 && nj == 2) {
				jac_t = jac_d_ij
			} else {
				jac_t = rbind(jac_t, jac_d_ij)
			}

		}
	}

}
dev.off()

# save
write.table(jac_t, sprintf("%s/module_communities.pairwise_jaccard.csv", ing_fn), sep = "\t", quote = FALSE, rownames = FALSE)


#### Find communities of modules based on jaccard overlap ####

# dataframe to network
jac_t = read.table(sprintf("%s/module_communities.pairwise_jaccard.csv", ing_fn), header = TRUE)
jac_g = igraph::graph_from_data_frame(jac_t [ jac_t$weight >= 0.125, ], directed = FALSE)

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

# save communities
write.table(jac_g_components_d,   sprintf("%s/module_communities.raw.csv", ing_fn), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(jac_g_components_d_a, sprintf("%s/module_communities.raw.annot.csv", ing_fn), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

