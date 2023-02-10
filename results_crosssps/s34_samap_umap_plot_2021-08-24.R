#### Input ####

# libraries
require("scales")

# list of comparisons
com_list = list(
	c("Tadh","TrH2"),
	c("Tadh","Hhon"),
	c("Tadh","HoiH23"),
	c("TrH2","Hhon"),
	c("TrH2","HoiH23"),
	c("Hhon","HoiH23")
)

out_fn = "results_samap_it4/"

cell_type_colors = c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen","deepskyblue","cadetblue1","mediumblue","darkviolet","violet")
color_palette = colorRampPalette(cell_type_colors)

for (n in 1:length(com_list)) {
	
	# species data
	com = com_list[[n]]
	sps1 = com[1]
	sps2 = com[2]

	# load umap coordinates
	ucs = read.table(sprintf("%s/csps_samap.%s-%s.cts.transfer_umap.csv", out_fn, sps1, sps2), header = TRUE, sep = "\t")
	ix_sps1 = which(ucs$species == sps1)
	ix_sps2 = which(ucs$species == sps2)
	
	# plot
	pdf(sprintf("%s/csps_samap.%s-%s.cts.transfer_umap.pdf",out_fn, sps1, sps2), width = 8, height = 12)
	par(mfrow=c(3,2))
	plot(  ucs$X_umap1[ix_sps2], ucs$X_umap2[ix_sps2], col = scales::alpha("gray", 0.5), cex = 0.5, xlab = "UMAP 1", ylab = "UMAP 2")
	points(ucs$X_umap1[ix_sps1], ucs$X_umap2[ix_sps1], col = scales::alpha("blue", 0.5), cex = 0.5)
	title(main = sprintf("UMAP %s-%s, by species", sps1, sps2))
	plot(0,0, xlab = "", ylab = "", col = NULL, yaxt = NULL, xaxt = NULL, axes = FALSE)
	legend("topleft", legend = c(sps1, sps2), col = c("blue","gray"), pch = 1, cex = 0.7, bty = "n")
	
	# color according to sps1 stage
	factr_vec = as.factor(ucs[,sprintf("%s_cell_type", sps1)][ix_sps1])
	color_vec = color_palette(nlevels(factr_vec)) [ as.numeric(factr_vec) ]
	plot(  ucs$X_umap1[ix_sps2], ucs$X_umap2[ix_sps2], col = scales::alpha("gray", 0.5), cex = 0.5, xlab = "UMAP 1", ylab = "UMAP 2")
	points(ucs$X_umap1[ix_sps1], ucs$X_umap2[ix_sps1], col = color_vec, cex = 0.5)
	title(main = sprintf("UMAP %s-%s, by %s stage", sps1, sps2, sps1))
	plot(0,0, xlab = "", ylab = "", col = NULL, yaxt = NULL, xaxt = NULL, axes = FALSE)
	legend("topleft", legend = levels(factr_vec), col = color_palette(nlevels(factr_vec)), pch = 1, cex = 0.7, bty = "n", ncol = max(ceiling(nlevels(factr_vec) / 20),1))
	
	# color according to sps2 stage
	factr_vec = as.factor(ucs[,sprintf("%s_cell_type", sps2)][ix_sps2])
	color_vec = color_palette(nlevels(factr_vec)) [ as.numeric(factr_vec) ]
	plot(  ucs$X_umap1[ix_sps1], ucs$X_umap2[ix_sps1], col = scales::alpha("gray", 0.5), cex = 0.5, xlab = "UMAP 1", ylab = "UMAP 2")
	points(ucs$X_umap1[ix_sps2], ucs$X_umap2[ix_sps2], col = color_vec, cex = 0.5)
	title(main = sprintf("UMAP %s-%s, by %s stage", sps1, sps2, sps2))
	plot(0,0, xlab = "", ylab = "", col = NULL, yaxt = NULL, xaxt = NULL, axes = FALSE)
	legend("topleft", legend = levels(factr_vec), col = color_palette(nlevels(factr_vec)), pch = 1, cex = 0.7, bty = "n", ncol = max(ceiling(nlevels(factr_vec) / 20),1))
	
	# close
	dev.off()
}
