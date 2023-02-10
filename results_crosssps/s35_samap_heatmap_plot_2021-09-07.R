#### Input ####

# libraries
require("scales")
source("../scripts/helper.R")

# list of comparisons
com_list = list(
	c("Tadh","TrH2"),
	c("Tadh","Hhon"),
	c("Tadh","HoiH23"),
	c("TrH2","Hhon"),
	c("TrH2","HoiH23"),
	c("Hhon","HoiH23")
)

# folders
out_fn = "results_samap_it4/"
ann_fn = "../results_scatlas/results_metacell_it4/"

# colors
heatmap_colors = c("gray95","orange","orangered2","#520c52")

# #### Pairwise heatmaps ####

cli_list = c("mcs","cts")
cln_list = c("metacell","cell_type")

for (c in 1:length(cli_list)) {
	
	cli = cli_list[c]
	cln = cln_list[c]
	
	for (n in 1:length(com_list)) {
		
		# species data
		com = com_list[[n]]
		sp1 = com[1]
		sp2 = com[2]
		
		# load SAMap transfer scores table
		message(sprintf("csps SAMap plot | %s-%s samap heatmap, %s", sp1, sp2, cli))
		sam_m = read.table(sprintf("%s/csps_samap.%s-%s.%s.ntop0.score_matrix.csv", out_fn, sp1, sp2, cli), sep = "\t", header = TRUE, row.names = 1)
		sam_m = t(sam_m)
		
		# drop species names from cell types and clean names
		sam_m = sam_m [ grepl(sp1, rownames(sam_m)),  ]
		sam_m = sam_m [ ,grepl(sp2, colnames(sam_m))  ]
		
		if (cli != "lei") {
			
			# load annotation
			met1 = read.table(sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp1), sep = "\t", header = TRUE)
			met2 = read.table(sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp2), sep = "\t", header = TRUE)
			
			# get clean metadata table
			met1_f = unique(met1 [ , c(cln, "color") ])
			met1_f[,cln] = paste(sp1, met1_f[,cln], sep = "_")
			met2_f = unique(met2 [ , c(cln, "color") ])
			met2_f[,cln] = paste(sp2, met2_f[,cln], sep = "_")
			
			# reorder
			sam_m = sam_m [ met1_f[,1], met2_f[,1] ]
			
		} else {
			
			sam_m_cor_row = hclust(as.dist(1 - cor(t(sam_m)) / 2))$order
			sam_m_cor_col = hclust(as.dist(1 - cor(sam_m) / 2))$order
			# reorder
			sam_m = sam_m [ sam_m_cor_row, sam_m_cor_col ]
			# no annot available
			met1_f = data.frame(cln = rownames(sam_m), color = "gray")
			met2_f = data.frame(cln = colnames(sam_m), color = "gray")
			
		}
		
		# open pdf
		pdf(sprintf("%s/csps_samap.%s-%s.%s.ntop0.score_matrix.pdf", out_fn, sp1, sp2, cli),
				width = 5 + round(ncol(sam_m) / 20), 
				height = 5 + round(nrow(sam_m) / 20))
		
		# heatmap
		csps_hm_sam = csps_plot_annotated_matrix(
			mat = sam_m,
			name = sprintf("Score\n%s",cli),
			heatmap_colors = heatmap_colors,
			use_raster = FALSE,
			min_val = 0.0,
			max_val = 1.0,
			row_title = sp1,
			col_title = sp2,
			fontsize = 6,
			row_annot = list(met1_f),
			col_annot = list(met2_f),
			row_cluster = FALSE,
			col_cluster = FALSE,
			do_dotplot = FALSE, 
			cex_dotplot = 0.02)
		tryCatch(print(csps_hm_sam), error = function(e) { 
			message("Caught an error!") 
		})
		
		# heatmap, only peptidergic
		if (cli == "cts") {
			
			boo_rows = grepl("peptidergic", rownames(sam_m)) & ! grepl("trans", rownames(sam_m))
			boo_cols = grepl("peptidergic", colnames(sam_m)) & ! grepl("trans", colnames(sam_m))
			
			csps_hm_sam = csps_plot_annotated_matrix(
				mat = sam_m [ boo_rows, boo_cols ],
				name = sprintf("Score\n%s",cli),
				heatmap_colors = heatmap_colors,
				use_raster = FALSE,
				min_val = 0.0,
				max_val = 1.0,
				row_title = sp1,
				col_title = sp2,
				fontsize = 6,
				row_annot = list(met1_f [ boo_rows, ]),
				col_annot = list(met2_f [ boo_cols, ]),
				row_cluster = FALSE,
				col_cluster = FALSE,
				do_dotplot = FALSE, 
				cex_dotplot = 0.02)
			tryCatch(print(csps_hm_sam), error = function(e) { 
				message("Caught an error!") 
			})
		}
		
		dev.off()
		
	}
	
	
}


#### All species heatmaps ####

sp1="Tadh"
sp2="TrH2"
sp3="Hhon"
sp4="HoiH23"
cli_list = c("cts")
cln_list = c("cell_type")

for (c in 1:length(cli_list)) {
	
	# clustering info
	cli = cli_list[c]
	cln = cln_list[c]
	
	# load SAMap transfer scores table
	message(sprintf("csps SAMap plot | 4sps samap heatmap, %s", cli))
	sam_m = read.table(sprintf("%s/csps_samap.4sps.%s.ntop0.score_matrix.csv", out_fn, cli), sep = "\t", header = TRUE, row.names = 1)
	sam_m = t(sam_m)
	sam_m = sam_m [ !grepl("trans", rownames(sam_m)), !grepl("trans", colnames(sam_m)) ]
	
	if (cli != "lei") {
		
		# load annotation
		met1 = read.table(sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp1), sep = "\t", header = TRUE)
		met2 = read.table(sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp2), sep = "\t", header = TRUE)
		met3 = read.table(sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp3), sep = "\t", header = TRUE)
		met4 = read.table(sprintf("%s/annotation_mc.%s.it4.reordered.tsv", ann_fn, sp4), sep = "\t", header = TRUE)
		
		# get clean metadata table
		met1_f = unique(met1 [ , c(cln, "color") ])
		met2_f = unique(met2 [ , c(cln, "color") ])
		met3_f = unique(met3 [ , c(cln, "color") ])
		met4_f = unique(met4 [ , c(cln, "color") ])
		met1_f[,cln] = paste(sp1, met1_f[,cln], sep = "_")
		met2_f[,cln] = paste(sp2, met2_f[,cln], sep = "_")
		met3_f[,cln] = paste(sp3, met3_f[,cln], sep = "_")
		met4_f[,cln] = paste(sp4, met4_f[,cln], sep = "_")
		met1_f = met1_f [ !grepl("trans", met1_f[,1]), ]
		met2_f = met2_f [ !grepl("trans", met2_f[,1]), ]
		met3_f = met3_f [ !grepl("trans", met3_f[,1]), ]
		met4_f = met4_f [ !grepl("trans", met4_f[,1]), ]
		
		# reorder
		sam_m = sam_m [ c(met1_f[,1], met2_f[,1], met3_f[,1], met4_f[,1]), c(met1_f[,1], met2_f[,1], met3_f[,1], met4_f[,1]) ]

		# metadata
		metg_f = data.frame(
			cln = c(met1_f[,1], met2_f[,1], met3_f[,1], met4_f[,1]), 
			color = c(met1_f[,2], met2_f[,2], met3_f[,2], met4_f[,2])
		)
		
	} else {
		
		sam_m_cor_row = hclust(as.dist(1 - cor(t(sam_m)) / 2))$order
		sam_m_cor_col = hclust(as.dist(1 - cor(sam_m) / 2))$order
		# reorder
		sam_m = sam_m [ sam_m_cor_row, sam_m_cor_col ]
		# no annot available
		vector_sps = gsub("_.*", "", rownames(sam_m))
		vector_col = rainbow(length(unique(vector_sps)), alpha = 1) [ as.factor(vector_sps) ]
		metg_f = data.frame("species" = vector_sps, color = vector_col)
		
		
	}

	# ones at diagonal?
	diag(sam_m) = 1
	
	# open pdf
	pdf(sprintf("%s/csps_samap.4sps.%s.ntop0.score_matrix.pdf", out_fn, cli),
			width = 5 + round(ncol(sam_m) / 20), 
			height = 5 + round(nrow(sam_m) / 20))
	
	# heatmap
	sam_m_na = sam_m
	sam_m_na [ upper.tri(sam_m_na) ] = NA
	csps_hm_sam = csps_plot_annotated_matrix(
		mat = sam_m_na,
		name = sprintf("Score\n%s",cli),
		heatmap_colors = heatmap_colors,
		use_raster = FALSE,
		min_val = 0.0,
		max_val = 1.0,
		row_title = cln,
		col_title = cln,
		fontsize = 6,
		row_annot = list(metg_f),
		col_annot = list(metg_f),
		row_cluster = FALSE,
		col_cluster = FALSE,
		do_dotplot = FALSE, 
		cex_dotplot = 0.02)
	tryCatch(print(csps_hm_sam), error = function(e) { 
		message("Caught an error!") 
	})
	
	
	
	# heatmap, only peptidergic
	if (cli == "cts") {
		
		sam_m_na = sam_m
		sam_m_na [ upper.tri(sam_m_na) ] = NA
		boo_rows = grepl("peptidergic", rownames(sam_m_na)) & ! grepl("trans", rownames(sam_m_na))
		boo_cols = grepl("peptidergic", colnames(sam_m_na)) & ! grepl("trans", colnames(sam_m_na))
	
		
		csps_hm_sam = csps_plot_annotated_matrix(
			mat = sam_m_na [ boo_rows, boo_cols ],
			name = sprintf("Score\n%s",cli),
			heatmap_colors = heatmap_colors,
			use_raster = FALSE,
			min_val = 0.0,
			max_val = 1.0,
			row_title = sp1,
			col_title = sp2,
			fontsize = 6,
			row_annot = list(metg_f [ boo_rows,]),
			col_annot = list(metg_f [ boo_cols,]),
			row_cluster = FALSE,
			col_cluster = FALSE,
			do_dotplot = FALSE, 
			cex_dotplot = 0.02)
		tryCatch(print(csps_hm_sam), error = function(e) { 
			message("Caught an error!") 
		})
	}
		
	
	dev.off()
	
	
}


message("all done!")


## alluvial plots

library("reshape2")

# Specify id.vars: the variables to keep but not split apart on
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

for (spi in sps_list) {

	sam_m2 = sam_m [ , grepl(spi, colnames(sam_m)) ]
	sam_d2 = reshape2::melt(sam_m2, id.vars = colnames(sam_m2), value.name = "score")
	sam_d2_f = sam_d2 [ sam_d2$score > 0.15, ]
	sam_d2_f = merge(sam_d2_f, metg_f, by.x = "Var2", by.y = "cln", all.x = TRUE, all.y = FALSE)
	# sam_d2_f = sam_d2_f [ !(grepl("_trans_", sam_d2_f[,1]) | grepl("_trans_", sam_d2_f[,2])),  ]
	sam_d2_f = sam_d2_f [ order(sam_d2_f[,1], -sam_d2_f[,3]) , ]

	pdf(sprintf("results_samap_it4/alluvial.%s.pdf", spi), height = 4, width = 7)
	layout(matrix(1:3, ncol = 3))
	for (spj in sps_list [ sps_list != spi ]) {
		sam_d2_ff = sam_d2_f
		sam_d2_ff = sam_d2_ff [ grepl(spj, sam_d2_ff[,"Var1"]), ]
		colnames(sam_d2_ff) = c(spi, spj, "score", "color")
		alluvial::alluvial(sam_d2_ff[,c(1,2)], freq = sam_d2_ff$score, col = sam_d2_ff$color, cex = 0.5, gap.width = 0.4, blocks = TRUE, cw = 0.2)
	}
	dev.off()

}