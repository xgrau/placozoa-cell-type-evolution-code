# sources
library("pheatmap")
library("zoo")

select_top_markers = function(matrix, matrix_thr = 0, n_top_markers = 20, n_markers_rollmean = 2) {
	
	markers = unique(as.vector(unlist( apply(matrix, 2, function(c) { names( head(sort(-c[ c >= matrix_thr ]), n = n_top_markers ) ) } ))))
	markers_order = order(apply(matrix[markers,], 1, function(r) which.max(rollmean(r, n_markers_rollmean) )))
	markers_ordered = markers [ markers_order ]
	return(markers_ordered)
	
}


# read original algs
to = read.table("results_macrosynteny_all/original_algs_simakov22.long.csv", sep = "\t", fill = TRUE, header = FALSE)
colnames(to) = c("og","alg","gene")
to = to [ grepl("^Em", to$gene), ]
to$gene = paste("Emue", to$gene, sep = "_")
# to = to [ , c(2,3) ]

# read hgs
hg = read.table("results_macrosynteny_all/mcl.all.out.txt", sep = "\t", col.names = c("homology_group","gene"))
hg = hg [ grepl("^Emue", hg$gene), ]
hg$gene = gsub("\\.t\\d+$","", hg$gene)

# colorscale
col_blue = colorRampPalette(interpolate="l",c("gray95", "deepskyblue","dodgerblue3","midnightblue"))
numcol=20



# read algs
tt_fl = list.files(path = "results_macrosynteny_all", pattern = "alg.emu.*.csv", full.names = TRUE)

for (tt_fn in tt_fl) {
	
	# read algs and merge
	tt = read.table(tt_fn, sep = "\t", header = TRUE)
	ti = gsub(".csv$","",basename(tt_fn))
	ta = merge(hg, tt, by = "homology_group", all.x = FALSE, all.y = TRUE)
	ta = merge(ta, to, by = "gene", all.x = TRUE, all.y = TRUE)

	# order
	ma = table(ta$alg, ta$ancestral_linkage_group)
	ma_top = select_top_markers(t(ma), 0, 20, 1)

	pdf(sprintf("results_macrosynteny_all/original_algs_simakov22.%s.pdf", ti), width = 8, height = 8)
	pheatmap(
		ma[,ma_top],
		color = col_blue(numcol),
		breaks = seq(0,20,length.out = numcol + 1), 
		cluster_rows = FALSE,
		cluster_cols = FALSE,
		cellwidth = 10, cellheight = 10, na_col = "grey", 
		border_color = "white",
		display_numbers = TRUE,
		number_format = "%i",
		main = "alg overlaps"
	)
	dev.off()
	
}

# to_alg = to[,2]
# to_gen = to[,3:ncol(to)]
# to_gen_l = apply(to_gen,1, function(v) { 
# 	v = v [ v != "" ]
# 	v = sort(unique(v))
# 	v = v [ grep("^Em", v) ]
# 	v = paste("Emue", v, sep = "_")
# })
# names(to_gen_l) = to_alg
