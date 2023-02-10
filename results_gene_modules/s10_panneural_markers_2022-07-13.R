#### Input ####

# libraries
library("metacell")
suppressMessages(source("../scripts/helper.R"))
library("scales")

# paths
out_fn = "results_panneural_markers"
dir.create(out_fn, showWarnings = FALSE)

# sps lists
sps_list_plc = c("Tadh","TrH2","Hhon","HoiH23")
sps_list_oth = sps_list = c("Mmus","Dmel", "Nvec","Hvul","Spis", "Mlei", "Spolac")
sps_list_all = c(sps_list_plc, sps_list_oth)

# parameters to select markers:
# - genes expressed with fc>2 in at least 10% of the neural metacells for that species
fc_thr = 2.0
mc_fra = 0.1


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
ogg_l$species = factor(ogg_l$species, levels = sps_list_all)
ogg_l = ogg_l [ !is.na(ogg_l$species), ]

# clean spongilla
ogg_l [ ogg_l$species == "Spolac", "gene"] = gsub("_i\\d+_.*","", ogg_l [ ogg_l$species == "Spolac", "gene"])

# vector dictionary
ogg_v = ogg_l$OG
names(ogg_v) = ogg_l$gene



## Placozoans ##

# paths input data
inp_fn = "../results_scatlas/results_metacell_it4/"
mc_sprintf_string = "%s_it4"
ctt_sprtinf_string = "%s/annotation_mc.%s.it4.reordered.tsv"

# loop
for (spi in sps_list_plc) {
	
	# load data
	run_name = sprintf("scdr_%s", spi)
	message(sprintf("%s | %s load", out_fn, spi))
	suppressMessages(metacell::scdb_init("../results_scatlas/data/scdb/",force_reinit=TRUE))
	mc = metacell::scdb_mc( sprintf(mc_sprintf_string, run_name))
	
	# save as transcripts
	rownames(mc@mc_fp) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), rownames(mc@mc_fp), t2g = FALSE)
	
	# subset to genes with orthologs
	mc_f = mc@mc_fp [ rownames(mc@mc_fp) %in% names(ogg_v), ]
	
	# load cell type annotations
	ctt_fn = sprintf(ctt_sprtinf_string, inp_fn, spi)
	ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")
	
	# neural genes
	mc_neu_mcs = ctt [ ctt$broad_cell_type == "peptidergic", "metacell" ]
	mc_neu_fp  = mc_f [ , mc_neu_mcs ]
	var_genes_neu = names(which(apply(as.data.frame(mc_neu_fp), 1, function(r) sum(r > fc_thr) >= max(1, ceiling(mc_fra * length(r))) )))
	var_genes_neu = data.frame("gene" = var_genes_neu, "orthogroup" = ogg_v [ var_genes_neu ])
	write.table(var_genes_neu, sprintf("%s/markers.peptidergic.%s.txt", out_fn, spi), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
	message(sprintf("%s | %s n=%i markers in neural", out_fn, spi, nrow(var_genes_neu)))

}



## Outgroup species ##

# paths to input data
inp_fn = "../results_scatlas/data/scdb_outgroups/"
mc_sprintf_string = "%s"
ctt_sprtinf_string = "%s/annot.%s.tsv"

# loop over species
for (spi in sps_list_oth) {
	
	# load data
	run_name = spi
	message(sprintf("%s | %s load", out_fn, spi))
	suppressMessages(metacell::scdb_init("../results_scatlas/data/scdb_outgroups/",force_reinit=TRUE))
	mc = metacell::scdb_mc( sprintf(mc_sprintf_string, run_name))
	
	# save as transcripts
	if (spi %in% c("Mmus","Mlei")) {
		rownames(mc@mc_fp) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), rownames(mc@mc_fp), t2g = FALSE)
	} else if (spi == "Dmel") {
		rownames(mc@mc_fp) = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", spi), rownames(mc@mc_fp), t2g = FALSE, gene_id = "gene_name")
	}

	# subset to genes with orthologs
	mc_f = mc@mc_fp [ rownames(mc@mc_fp) %in% names(ogg_v), ]
	
	# load cell type annotations
	ctt_fn = sprintf(ctt_sprtinf_string, inp_fn, spi)
	ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")
	
	
	# select neural metacells
	if (spi %in% c("Spolac")) {

		# neural genes
		mc_neu_mcs = ctt [ ctt$broad_cell_type == "neuroid", "metacell" ]
		mc_neu_fp  = mc_f [ , mc_neu_mcs ]
		var_genes_neu = names(which(apply(as.data.frame(mc_neu_fp), 1, function(r) sum(r > fc_thr) >= max(1, ceiling(mc_fra * length(r))) )))
		var_genes_neu = data.frame("gene" = var_genes_neu, "orthogroup" = ogg_v [ var_genes_neu ])
		write.table(var_genes_neu, sprintf("%s/markers.neuroid.%s.txt", out_fn, spi), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		message(sprintf("%s | %s n=%i markers in neuroid", out_fn, spi, nrow(var_genes_neu)))

	} else {
		
		# neural genes
		mc_neu_mcs = ctt [ ctt$broad_cell_type == "neuron", "metacell" ]
		mc_neu_fp  = mc_f [ , mc_neu_mcs ]
		var_genes_neu = names(which(apply(as.data.frame(mc_neu_fp), 1, function(r) sum(r > fc_thr) >= max(1, ceiling(mc_fra * length(r))) )))
		var_genes_neu = data.frame("gene" = var_genes_neu, "orthogroup" = ogg_v [ var_genes_neu ])
		write.table(var_genes_neu, sprintf("%s/markers.neuron.%s.txt", out_fn, spi), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		message(sprintf("%s | %s n=%i markers in neural", out_fn, spi, nrow(var_genes_neu)))

	}
	
}



## Panneural orthogroup matrix ##

# get vector of OGs to check
message(sprintf("%s | create OG matrix", out_fn))
ogs_v   = c()
for (spi in sps_list_all) {
	
	if (spi %in% sps_list_plc) {
		mar_fn = sprintf("%s/markers.peptidergic.%s.txt", out_fn, spi)
	} else if (spi %in% c("Spolac","Aque")) {
		mar_fn = sprintf("%s/markers.neuroid.%s.txt", out_fn, spi)
	} else {
		mar_fn = sprintf("%s/markers.neuron.%s.txt", out_fn, spi)
	}
	
	mar_v = read.table(mar_fn, sep = "\t", header = TRUE)
	ogs_v = unique(sort(c(ogs_v, mar_v[,2])))
	
}

ogs_order = unique(ogg_l$OG)
ogs_order = ogs_order [ ogs_order %in% ogs_v ]

# fill matrix
ogs_m = matrix(nrow = length(ogs_v), ncol = length(sps_list_all))
rownames(ogs_m) = ogs_order
colnames(ogs_m) = sps_list_all

mar_t = data.frame()
for (spi in sps_list_all) {
	
	if (spi %in% sps_list_plc) {
		mar_fn = sprintf("%s/markers.peptidergic.%s.txt", out_fn, spi)
	} else if (spi %in% c("Spolac","Aque")) {
		mar_fn = sprintf("%s/markers.neuroid.%s.txt", out_fn, spi)
	} else {
		mar_fn = sprintf("%s/markers.neuron.%s.txt", out_fn, spi)
	}
	
	mar_v = read.table(mar_fn, sep = "\t", header = TRUE)
	mar_v$orthogroup = factor(mar_v$orthogroup, levels = ogs_order)
	ogs_m[,spi] = table(mar_v$orthogroup)
	
	# save
	mar_t = rbind(mar_t, mar_v[,c("gene","orthogroup")])
	
}

# save as presence matrix
ogs_m = data.frame(ogs_m)
ogs_m$orthogroup = rownames(ogs_m)
ogs_m = ogs_m[,c(ncol(ogs_m),1:(ncol(ogs_m)-1))]
write.table(ogs_m, sprintf("%s/matrix.all.counts.ogexpression.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# save as long format
mar_t$species = gsub("_.*","", mar_t$gene)
mar_t$species = factor(mar_t$species, levels = sps_list_all)
mar_t = mar_t [ !is.na(mar_t$orthogroup), ]
mar_t = mar_t [ order(mar_t$orthogroup, mar_t$species), ]
write.table(mar_t, sprintf("%s/matrix.all.long.ogexpression.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# # filtered version
# ogs_m_f = ogs_m
# # ogs_m_f = ogs_m_f [ apply(ogs_m_f[,sps_list_plc],     1, function(r) sum(r > 0) >= 3 ), ]
# ogs_m_f = ogs_m_f [ apply(ogs_m_f[,c("Dmel","Mmus")], 1, function(r) sum(r > 0) >= 1 ), ]
# ogs_m_f_v = rownames(ogs_m_f)
# write.table(ogs_m_f, sprintf("%s/matrix.fil.counts.ogexpression.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# # tfs only
# tfs_v = c()
# for (spi in sps_list_plc) {
# 	tfs_v = unique(c(tfs_v, read.table(sprintf("../data/gene_annotations/tfs.%s_genes.txt", spi))[,1]))
# }
# ogs_is_tf = unique(ogg_v [ names(ogg_v) %in% tfs_v ])
# ogs_is_tf = intersect(ogs_is_tf, rownames(ogs_m))
# write.table(ogs_m[ogs_is_tf,], sprintf("%s/matrix.tfs.counts.ogexpression.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# annotate orthogroups according to placozoans or Mus (if no placozoans)
# placozoans
ref_d = data.frame()
for (spi in sps_list_plc) {
	message(sprintf("%s | annotate OGs to %s", out_fn, spi))
	ref_txs = read.table(sprintf("%s/markers.peptidergic.%s.txt", out_fn, spi), sep = "\t", header = TRUE)
	
	# load gene annotations (based on transcripts, need to be changed to genes)		
	gene_annot = read.table(sprintf("../data/reference/%s_long.pep.annotations.csv", spi), sep = "\t", row.names = 1)
	# load TF-specific annotations
	dic_genetf = read.table(sprintf("../data/gene_annotations/tfs.%s_genes.curated.csv", spi), sep = "\t", row.names = 1, col.names = c("gene","OG"))
	dic_genetf_gene_name = gsub("like:","like_",dic_genetf$OG)
	dic_genetf_gene_name = stringr::str_split(dic_genetf_gene_name, pattern = ":", simplify = TRUE)[,2]
	names(dic_genetf_gene_name) = rownames(dic_genetf)
	gene_annot_tfs = merge(gene_annot, dic_genetf_gene_name, by.x = 0, by.y = 0, all.x = TRUE, all.y = FALSE)[,"y"]
	gene_annot [ !is.na(gene_annot_tfs) , "V2" ] = gene_annot_tfs [ !is.na(gene_annot_tfs) ]
	ref_di = data.frame("gene" = rownames(gene_annot), "gene_name" = gene_annot[,1])
	ref_d = rbind(ref_d, merge(ref_di, ogg_l, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE))
	
}

# mus
message(sprintf("%s | annotate OGs to Mmus", out_fn))
ref_txs = read.table(sprintf("%s/markers.neuron.%s.txt", out_fn, "Mmus"), sep = "\t", header = TRUE)
ref_gff = rtracklayer::readGFF("../data/reference/Mmus_long.annot.gtf")
ref_di = data.frame("gene" = ref_gff$transcript, "gene_name" = ref_gff$gene_name)
ref_di = ref_di [ ref_di$gene %in% ref_txs$gene,  ]
ref_di = unique(ref_di)
ref_d = rbind(ref_d, merge(ref_di, ogg_l, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE))
ref_d = ref_d [ order(ref_d$species), ]
# merge
ref_d = unique(ref_d [ , c("gene_name","OG") ])
ref_d = aggregate(gene_name ~ OG, data = ref_d, function(v) paste( sort(unique(c(v))), collapse = "/") )
ref_d$gene_name = gsub("^/","",ref_d$gene_name)
ref_d$gene_name = gsub("/$","",ref_d$gene_name)
rownames(ref_d) = ref_d$OG
ref_d$orthogroup_name = paste(ref_d$OG, ref_d$gene_name, sep = ":")
ref_d_v = ref_d$orthogroup_name
names(ref_d_v) = ref_d$OG

# add annot and sort
ogs_annot = data.frame("orthogroup" = ogs_order, "orthogroup_name" = ref_d_v [ ogs_order ])
rownames(ogs_annot) = ogs_order
ogs_annot = ogs_annot [ ogs_order, ]

# save
write.table(ogs_annot,             sprintf("%s/matrix.all.ognames.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


## Conservation orthogroup matrix ##

# fill ortholog presence matrix
ogp_m = matrix(nrow = length(ogs_v), ncol = length(sps_list_all))
rownames(ogp_m) = ogs_order
colnames(ogp_m) = sps_list_all

map_t = data.frame()
for (spi in sps_list_all) {
	
	map_v = data.frame(
		"gene" = ogg_l [ ogg_l$species == spi, "gene" ],
		"orthogroup" = factor(ogg_l [ ogg_l$species == spi, "OG" ], levels = ogs_order)
	)
	ogp_m[,spi] = table(map_v$orthogroup)
	
	map_t = rbind(map_v, map_t)
	
}

# long format
map_t$species = gsub("_.*","", map_t$gene)
map_t$species = factor(map_t$species, levels = sps_list_all)
map_t = map_t [ !is.na(map_t$orthogroup), ]
map_t = map_t [ order(map_t$orthogroup, map_t$species), ]

# save
ogp_m = data.frame(ogp_m)
ogp_m$orthogroup = rownames(ogp_m)
ogp_m = ogp_m[,c(ncol(ogp_m),1:(ncol(ogp_m)-1))]
write.table(ogp_m,             sprintf("%s/matrix.all.counts.ogpresence.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(map_t,             sprintf("%s/matrix.all.long.ogpresence.csv",   out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(ogp_m[ogs_m_f_v,], sprintf("%s/matrix.fil.counts.ogpresence.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(ogp_m[ogs_is_tf,], sprintf("%s/matrix.tfs.counts.ogpresence.csv", out_fn), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

