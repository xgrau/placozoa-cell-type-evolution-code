#### Input ####

# libraries
source("../scripts/helper.R")
library("igraph")
range01 = function(x) { (x-min(x)) / (max(x)-min(x)) }
graphics.off()

# input
heatmap_colors = c("white","orange","orangered2","#520c52")
icc_fn = "results_alignment_icc/"
ann_fn = "../results_scatlas/results_metacell_it4/"
inp_fn = "results_metacell_it4_transphyla/"
out_fn = "results_metacell_it4_ctpairwise/"
dir.create(out_fn)

# list of comparisons
list_comp = list(
	c("Tadh","Nvec"),
	c("Tadh","Hvul"),
	c("Tadh","Spis"),
	c("Tadh","Mmus"),
	c("Tadh","Dmel"),
	c("Tadh","Mlei"),
	c("Tadh","Spolac"),
	c("TrH2","Nvec"),
	c("TrH2","Hvul"),
	c("TrH2","Spis"),
	c("TrH2","Mmus"),
	c("TrH2","Dmel"),
	c("TrH2","Mlei"),
	c("TrH2","Spolac"),
	c("Hhon","Nvec"),
	c("Hhon","Hvul"),
	c("Hhon","Spis"),
	c("Hhon","Mmus"),
	c("Hhon","Dmel"),
	c("Hhon","Mlei"),
	c("Hhon","Spolac"),
	c("HoiH23","Nvec"),
	c("HoiH23","Hvul"),
	c("HoiH23","Spis"),
	c("HoiH23","Mmus"),
	c("HoiH23","Dmel"),
	c("HoiH23","Mlei"),
	c("HoiH23","Spolac")
)

sps_list = c("Tadh","TrH2","Hhon","HoiH23")

list_bil = c("Dmel","Mmus")
list_cni = c("Nvec","Spis","Hvul")
list_oth = c("Mlei","Spolac")

# list of focus levels
focus_list = list(
	"bct" = "broad_cell_type"
)

focid = "bct"
focus = "broad_cell_type"
list_focus_cts = c("neuron")

cor_method = "wpearson"
fp_thr = c(1.5)
fraction_shared_ogs = 0.3


for (ctf in list_focus_cts[1]) {
	
	ovs_d = data.frame()
	
	for (com in list_comp) {
		
		# define input
		sp1 = com[1]
		sp2 = com[2]
		
		if (sp2 == "Spolac"){
			ctf = "neuroid"
		} else {
			ctf = "neuron"
		}

		# load cross-species marker data (from `s01`)
		message(sprintf("csps icc | load %s %s-%s csps data...",focus, sp1, sp2))
		
		# load cell type data
		ctt_sp1_fn = sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp1)
		ctt_sp1 = read.table(ctt_sp1_fn, header = TRUE, comment.char = "", sep = "\t")
		if (sp2 %in% c("Spis","Nvec","Hvul","Dmel","Mmus","Spolac","Mlei","Aque")) {
			ctt_sp2_fn = sprintf("../results_scatlas/data/scdb_outgroups/annot.%s.tsv", sp2)
		} else {
			ctt_sp2_fn = sprintf(ctt_sprtinf_string, ann_fn, spi)
		}
		ctt_sp2 = read.table(ctt_sp2_fn, header = TRUE, comment.char = "", sep = "\t")
		
		# get color annotations for each species
		ann_sp1 = unique(ctt_sp1[,c(focus,"color")])
		ann_sp1 = ann_sp1 [ !duplicated(ann_sp1[,focus]), ]
		ann_sp2 = unique(ctt_sp2[,c(focus,"color")])
		ann_sp2 = ann_sp2 [ !duplicated(ann_sp2[,focus]), ]

		# load csps object
		csps = readRDS(sprintf("%s/csps_icc.%s.all.%s-%s.%s.fc%.2f.rds", inp_fn, focid, sp1, sp2, cor_method, fp_thr))
		
		# matrix
		csps_ovf = csps$cor_matrix
		# pairwise similarities
		if (sprintf("%s|%s", sp2, ctf) %in% colnames(csps_ovf)) {
			csps_ovf = t(csps$cor_matrix) #/ length(unique(unlist(csps$overlap_genes$sp2[[sprintf("%s|%s", sp2, ctf)]])))
			csps_ovf [ csps_ovf < 0 ] = 0
			dv = csps_ovf[sprintf("%s|%s", sp2, ctf),]
			dd = data.frame(to = names(dv), similarity = dv)
			dd$from = sprintf("%s|%s", sp2, ctf)
			dd$from_sps = sp2
			dd$to_sps = sp1
			ovs_d = rbind(ovs_d, dd)
			
		}
		
	}
	
	# plot network
	ovs_d_f = ovs_d [ ovs_d$similarity >= 0.2 & ovs_d$from_sps %in% c(list_bil, list_cni, list_oth),  ]
	ovs_d_f$weight = ovs_d_f$similarity
	ovs_d_f = ovs_d_f [ order(ovs_d_f$similarity), ]
	# ovs_d_f$from = gsub("^[^\\|]+\\|","",ovs_d_f$from)
	jac_g = igraph::graph_from_data_frame(ovs_d_f[,c("to","from","weight") ], directed = FALSE)
	# alluvial::alluvial(ovs_d_f[ , c("from","to") ], freq = ovs_d_f$similarity, col = "blue", alpha = 0.5, blocks = TRUE)
	
	# colors for nodes
	jac_g_v_color = rep("snow4", length(V(jac_g)))
	jac_g_v_color [ grepl("neuron",names(V(jac_g))) ] = "turquoise3"
	jac_g_v_color [ grepl("peptidergic",names(V(jac_g))) ] = "blue3"

	# colors for labels
	jac_g_v_collb = jac_g_v_color
	jac_g_v_collb = scales::alpha(colorspace::darken(jac_g_v_collb, 0.1), 1)

	set.seed(1)
	jac_g_layers = rep(3, length(names(V(jac_g))))
	jac_g_layers [ grepl("neuron", names(V(jac_g))) ] = 2
	jac_g_layers [ grepl("peptidergic", names(V(jac_g))) ] = 1
	jac_g_layers [ grepl("Mlei", names(V(jac_g))) ] = 4
	jac_g_layers [ grepl("Spolac", names(V(jac_g))) ] = 4
	jac_g_layers [ jac_g_layers == 0 ] = 3
	jac_g_lay = igraph::layout_with_sugiyama(jac_g, layers = jac_g_layers)
	# jac_g_lay = igraph::layout_with_fr(jac_g, niter = 5000)

	pdf(sprintf("%s/pairwise_ct_similarity.graph.%s.pdf", out_fn, "neuron"), height = 3, width = 3)
	igraph::plot.igraph(
		jac_g,
		vertex.size  = 8,
		vertex.color = jac_g_v_collb,
		vertex.frame.color = jac_g_v_collb,
		vertex.label.family = "sans",
		vertex.label.cex =  0.2,
		edge.arrow.size = 1,
		layout = jac_g_lay,
		edge.color = scales::alpha("snow3", 0.8),
		edge.width = (((E(jac_g)$weight - 0.2) ^ (1/1.1))) * 8
	)
	legend_weights = seq(0, 1, length.out = 11)
	legend("bottomleft", legend = legend_weights, lwd = ((legend_weights - 0.2) ^ (1/1.1)) * 8, bty = "n", cex = 0.5, title = "cor")
	dev.off()

	
	# save
	write.table(ovs_d, sprintf("%s/pairwise_ct_similarity.graph.%s.csv",out_fn, ctf), sep = "\t", quote = FALSE, row.names = FALSE)
	
}

message("All comparisons done!")



