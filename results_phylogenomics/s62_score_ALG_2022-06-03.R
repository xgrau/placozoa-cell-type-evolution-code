suppressMessages(library("rtracklayer"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/Cross_species_functions.R"))
graphics.off()

# out
out_fn = "results_macrosynteny_all"
dir.create(sprintf("%s/tests_chisq", out_fn), showWarnings = FALSE)


# define species lists
sps_list = read.table("species_list_synteny_blocks.long.txt")[,1]
sps_list = c("Tado", sps_list)

# function
table_to_matrix = function(table) {
  mat = matrix(table, nrow = nrow(table))
  rownames(mat) = rownames(table)
  colnames(mat) = colnames(table)
  return(mat)
}

# bin length and step
# applies to genes belonging to homology groups (others are ignored)
bin_w = 200
bin_s = 20
bin_r = bin_w / bin_s




list_comparisons = list(
	emuastrho = c("Emue","Astrub","Rhoesc"),
	emublarho = c("Emue",  "Bralan","Rhoesc"),
	emucgirho = c("Emue",  "Cgig","Rhoesc"),
	
	emuastami = c("Emue","Astrub","Amil"),
	emublaami = c("Emue",  "Bralan","Amil"),
	emucgiami = c("Emue",  "Cgig","Amil"),
	
	emuasthvu = c("Emue","Astrub","Hvul_v3"),
	emublahvu = c("Emue",  "Bralan","Hvul_v3"),
	emucgihvu = c("Emue",  "Cgig","Hvul_v3"),
	
	emuastnve = c("Emue","Astrub","Nvec_vc1.1"),
	emublanve = c("Emue",  "Bralan","Nvec_vc1.1"),
	emucginve = c("Emue",  "Cgig","Nvec_vc1.1")
)


# homology groups
hom_d = read.table(sprintf("%s/mcl.all.out.txt", out_fn), header = FALSE, col.names = c("homology_group", "transcript"))
hom_d$species = gsub("_.*","", hom_d$transcript)
hom_d$species [ hom_d$species == "Hvul" ] = "Hvul_v3"
hom_d$species [ hom_d$species == "Mlei" ] = "Mlei_v3"
hom_d$species [ hom_d$species == "Nvec" ] = "Nvec_vc1.1"
# hom_d = hom_d [ hom_d$species %in% sps_list , ]
hom_d = hom_d [ !duplicated(hom_d$transcript), ]
rownames(hom_d) = hom_d$transcript

# color palettes
categorical_colors      = colorspace::darken(c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen"))
categorical_colors_dark = colorspace::darken(c("deepskyblue","paleturquoise1","mediumblue","darkviolet","orchid1","thistle1"))
color_palette           = colorRampPalette(categorical_colors)
color_palette_dark      = colorRampPalette(categorical_colors_dark)


for (n in 1:length(list_comparisons)) {
	
	cmpnm = names(list_comparisons)[n]

	# load ancestral linkage groups
	alg_d = read.table(sprintf("%s/alg.%s.csv", out_fn, cmpnm), sep = "\t", header = TRUE)
	list_algs = unique(alg_d$ancestral_linkage_group)

	# color ALGs
	spr_colors = color_palette(n=length(list_algs))
	names(spr_colors) = unique(list_algs)

	# loop over queries	
	list_ovs = list()
	list_ovb = list()
	list_jac = list()
	for (i in 1:length(sps_list)) {
		
		spi = sps_list[i]

		message(sprintf("ALG %s | %s load gff", cmpnm, spi))
		if (spi == "Tadh") {
			gen_i = rtracklayer::readGFFAsGRanges(sprintf("../data/reference/%s_chr_long.annot.gtf", spi))
			gen_i$transcript_id = gsub("^Tadh_v2_","", gen_i$transcript_id)
		} else if (spi == "Tado") {
			gen_i = rtracklayer::readGFFAsGRanges(sprintf("../data/reference/Tadh_long.annot.gtf"))
		} else {
			gen_i = rtracklayer::readGFFAsGRanges(sprintf("../data/reference/%s_long.annot.gtf", spi))
		}
		gen_i = gen_i [ gen_i$type == "transcript" ]
		gen_d = as.data.frame(gen_i)
		gen_d = gen_d [ , c("seqnames","transcript_id") ]
		gen_d$homology_group = hom_d [ gen_d$transcript_id, "homology_group" ]
		gen_d = gen_d [ !is.na(gen_d$homology_group),  ]
		gen_d$seqnames = paste(spi, gen_d$seqnames, sep = ":")
		
		# filter out chromosomes with few genes
		gen_d = gen_d [ gen_d$seqnames %in% names(which(table(gen_d$seqnames) > 50)),  ]
		
		# chrom-level overlap
		# cross-tabulate
		mat_t = table_to_matrix(table(gen_d$seqnames, gen_d$homology_group))
		mat_t [ mat_t > 1 ] = 1

		# overlaps in content between chromosomes and ALGs
		ovs_i = matrix(nrow = length(list_algs), ncol = nrow(mat_t))
		for (ali in 1:length(list_algs)) {
			alg = list_algs[ali]
			hgs_in_alg = alg_d [ alg_d$ancestral_linkage_group == alg, "homology_group" ]
			hgs_in_alg = hgs_in_alg [ hgs_in_alg %in% hom_d [ hom_d$species == spi, "homology_group" ] ]
			mat_t_a = mat_t [ , hgs_in_alg [ hgs_in_alg %in% colnames(mat_t) ]]
			ovs_i [ ali , ] = rowSums(mat_t_a) / ncol(mat_t_a)
		}
		rownames(ovs_i) = list_algs
		colnames(ovs_i) = rownames(mat_t)
		
		# ovs_i_f = ovs_i [ 1:6, 1:5 ]
		
		# reorder chromosomes based on ALG alignment
		ovs_i_f_chrom_order = unlist(unique(apply(ovs_i, 1, function(r) names(which.max(r)))))
		ovs_i_f_chrom_order = c(ovs_i_f_chrom_order, setdiff(colnames(ovs_i), ovs_i_f_chrom_order) )
		ovs_i = ovs_i [ , ovs_i_f_chrom_order ]

		
		# binned jaccard heatmap
		message(sprintf("ALG %s | %s load bins", cmpnm, spi))
		mat_i = Matrix::readMM(file = sprintf("%s/bins/bin.all-%s.mat.csv", out_fn, spi))
		met_i = read.table(sprintf("%s/bins/bin.all-%s.dat.csv", out_fn, spi))
		hgs_i = read.table(sprintf("%s/bins/bin.all-%s.hgs.txt", out_fn, spi))[,1]
		colnames(mat_i) = hgs_i
		rownames(mat_i) = rownames(met_i)
		
		# subset to comparable hgs
		mat_i_f = mat_i
		mat_i_f [ mat_i_f > 1 ] = 1
		
		# jaccard
		message(sprintf("ALG %s | %s measure overlaps", cmpnm, spi))
		jac_i = matrix(nrow = length(list_algs), ncol = nrow(mat_i_f))
		ovb_i = matrix(nrow = length(list_algs), ncol = nrow(mat_i_f))
		for (ali in 1:length(list_algs)) {
			alg = list_algs[ali]
			hgs_in_alg = alg_d [ alg_d$ancestral_linkage_group == alg, "homology_group" ]
			hgs_in_alg = hgs_in_alg [ hgs_in_alg %in% hom_d [ hom_d$species == ifelse(spi == "Tado", "Tadh", spi), "homology_group" ] ]
			jac_along_bin = sapply(1:nrow(mat_i_f), function(b) {
				v = mat_i_f[b,]
				h = names(which(v > 0))
				h = h [ h %in% alg_d$homology_group ]
				jaccard_index(hgs_in_alg, h)
			})
			jac_i [ ali, ] = jac_along_bin
			ovs_along_bin = sapply(1:nrow(mat_i_f), function(b) {
				v = mat_i_f[b,]
				h = names(which(v > 0))
				h = h [ h %in% alg_d$homology_group ]
				length(intersect(hgs_in_alg, h))
			})
			ovb_i [ ali, ] = ovs_along_bin
		}
		rownames(jac_i) = list_algs
		colnames(jac_i) = rownames(mat_i_f)
		rownames(ovb_i) = list_algs
		colnames(ovb_i) = rownames(mat_i_f)
		
		
		# keep only chromosomes present in the overlaps matrix
		jac_i_chrs = stringr::str_split(colnames(jac_i), ":", simplify = TRUE)[,2] %in% stringr::str_split(colnames(ovs_i), ":", simplify = TRUE)[,2]
		jac_i = jac_i [ , jac_i_chrs ]
		ovb_i = ovb_i [ , jac_i_chrs ]
		
		# store matrices
		list_ovs[[spi]] = ovs_i
		list_jac[[spi]] = jac_i
		list_ovb[[spi]] = ovb_i
		
	}


	# loop heatmaps
	list_hms_chr = list()
	list_hms_jac = list()
	list_hms_ovb = list()
	for (i in 1:length(sps_list)) {
		
		# get data from previous loop
		spi = sps_list[i]
		ovs_i = list_ovs[[spi]]
		jac_i = list_jac[[spi]]
		ovb_i = list_ovb[[spi]]
		message(sprintf("ALG %s | loop heatmaps, %s", cmpnm, spi))

		# chromosome-level ALG heatmap
		# color strings
		chrs_row = rownames(ovs_i)
		chrs_col = stringr::str_split(colnames(ovs_i), ":", simplify = TRUE)[,2]
		chrs_row = factor(chrs_row)
		chrs_col = factor(chrs_col, levels = unique(stringr::str_split(colnames(ovs_i), ":", simplify = TRUE)[,2]))
		spi_colors = color_palette_dark(n=nlevels(chrs_col))
		colors_row = spr_colors [ chrs_row ]
		colors_col = spi_colors [ chrs_col ]
		names(colors_row) = chrs_row
		names(colors_col) = chrs_col
		list_hms_chr[[spi]] = plot_complex_heatmap(
			ovs_i, 
			name = "ovs",
			color_mat = c("gray95","lightskyblue2","deepskyblue","dodgerblue3","midnightblue"),
			color_min = 0, 
			color_max = 1,
			categories_row = chrs_row, 
			categories_col = chrs_col,
			colors_row = colors_row,
			colors_col = colors_col,
			cluster_row = FALSE,
			cluster_col = FALSE,
			name_row_show = TRUE,
			name_col_show = TRUE,
			show_legend_row = FALSE,
			show_legend_col = TRUE,
			separate_row = TRUE,
			separate_col = TRUE,
			both_sides_row = FALSE,
			use_raster = FALSE,
			fontsize = 6,
			title_col = spi, title_row = "ALG"
		)

		# heatmap of bin-wise jaccard values
		# color strings
		chrs_row = rownames(jac_i)
		chrs_col = stringr::str_split(colnames(jac_i), ":", simplify = TRUE)[,2]
		chrs_row = factor(chrs_row, levels = unique(chrs_row))
		chrs_col = factor(chrs_col, levels = stringr::str_split(colnames(ovs_i), ":", simplify = TRUE)[,2])
		spi_colors = color_palette_dark(n=nlevels(chrs_col))
		colors_row = spr_colors [ chrs_row ]
		colors_col = spi_colors [ chrs_col ]
		names(colors_row) = chrs_row
		names(colors_col) = chrs_col
		list_hms_jac[[spi]] = plot_complex_heatmap(
			jac_i, 
			name = "jac",
			color_mat = c("gray95","#d6e72e","#6fb600","#003f4d"),
			color_min = 0, 
			color_max = max(0.3, quantile(jac_i, 0.99)),
			categories_row = chrs_row, 
			categories_col = chrs_col,
			colors_row = colors_row,
			colors_col = colors_col,
			cluster_row = FALSE,
			cluster_col = FALSE,
			name_row_show = TRUE,
			name_col_show = TRUE,
			separate_row = TRUE,
			separate_col = TRUE,
			show_legend_row = FALSE,
			show_legend_col = TRUE,
			use_raster = FALSE,
			fontsize = 6,
			title_col = sprintf("%s gene blocks, jaccard", spi), title_row = sprintf("ALG")
		)
		
		# heatmap of bin-wise overlap counts
		# color strings
		chrs_row = rownames(ovb_i)
		chrs_col = stringr::str_split(colnames(ovb_i), ":", simplify = TRUE)[,2]
		chrs_row = factor(chrs_row, levels = unique(chrs_row))
		chrs_col = factor(chrs_col, levels = stringr::str_split(colnames(ovs_i), ":", simplify = TRUE)[,2])
		spi_colors = color_palette_dark(n=nlevels(chrs_col))
		colors_row = spr_colors [ chrs_row ]
		colors_col = spi_colors [ chrs_col ]
		names(colors_row) = chrs_row
		names(colors_col) = chrs_col
		list_hms_ovb[[spi]] = plot_complex_heatmap(
			ovb_i, 
			name = "num",
			color_mat = c("gray95","#d6e72e","#6fb600","#003f4d"),
			color_min = 0, 
			color_max = 20,
			categories_row = chrs_row, 
			categories_col = chrs_col,
			colors_row = colors_row,
			colors_col = colors_col,
			cluster_row = FALSE,
			cluster_col = FALSE,
			name_row_show = TRUE,
			name_col_show = TRUE,
			separate_row = TRUE,
			separate_col = TRUE,
			show_legend_row = FALSE,
			show_legend_col = TRUE,
			use_raster = FALSE,
			fontsize = 6,
			title_col = sprintf("%s gene blocks, counts", spi), title_row = sprintf("ALG")
		)
		
		# load bin coordinates
		met_i = read.table(sprintf("%s/bins/bin.all-%s.dat.csv", out_fn, spi))
		met_i_midpoint = met_i$start #+ (met_i$end - met_i$start) / 2
		names(met_i_midpoint) = rownames(met_i)

		
		# ALG contribution to each chromosome
		message(sprintf("ALG %s | loop ALG contribution, %s", cmpnm, spi))
		pdf(sprintf("%s/alg.%s.chrom-%s.pdf", out_fn, cmpnm, spi), height = 8, width = 18)
		layout(matrix(1:12, nrow = 3, byrow = FALSE))
		for (chr in colnames(ovs_i)) {
			
			# select bins in this chrom
			jac_i_chrs = sprintf("%s:%s" ,spi, stringr::str_split(colnames(jac_i), ":", simplify = TRUE)[,2])
			jac_i_i = as.matrix(jac_i [ , jac_i_chrs == chr ])
			ovb_i_i = as.matrix(ovb_i [ , jac_i_chrs == chr ])
			colnames(jac_i_i) = colnames(jac_i) [jac_i_chrs == chr]
			colnames(ovb_i_i) = colnames(jac_i) [jac_i_chrs == chr]
			
			# which algs to plot
			keep_alg_ixs = which(sapply(1:nrow(jac_i_i), function(r) max(jac_i_i[r,]) >= 0.05))
			keep_alg_nms = rownames(jac_i_i) [ keep_alg_ixs ]
			
			# create matrix
			jac_i_i_m = matrix(jac_i_i [ keep_alg_ixs , ], ncol = ncol(jac_i_i))
			ovb_i_i_m = matrix(ovb_i_i [ keep_alg_ixs , ], ncol = ncol(jac_i_i))
			colnames(jac_i_i_m) = colnames(jac_i_i)
			colnames(ovb_i_i_m) = colnames(jac_i_i)
			rownames(jac_i_i_m) = keep_alg_nms
			rownames(ovb_i_i_m) = keep_alg_nms
			do_plot = dim(jac_i_i_m)[1] > 0 & dim(jac_i_i_m)[2] > 0

			# plot
			if (do_plot) {
				
				# barplot aggregated
				jac_i_i_m_f = t(t(jac_i_i_m) / apply(t(jac_i_i_m), 1, sum))
				jac_i_i_m_f [ is.na(jac_i_i_m_f) ] = 0
				par(lwd = 0.5)
				barplot(jac_i_i_m_f, col = colors_row [ rownames(jac_i_i_m) ], las = 2, cex.names = 0.5, cex.main = 0.7, cex.axis = 0.7, ylab = "Sum Jaccard", cex.lab = 0.7, main = sprintf("ALG contribution to %s", chr), xlim = c(1,80), border = "white")
				par(lwd = 1)
				legend("topright", fill = colors_row [ rownames(jac_i_i_m) ], legend = rownames(jac_i_i_m), bty = "n", cex = 0.5)
				
				# counts in non-overlapping bins
				ovb_i_i_m_nvp = ovb_i_i_m[ , seq(1,ncol(ovb_i_i_m), by = bin_r) ]

				# chisq, observed and expected: non-overlapping bins
				if (!is.null(dim(ovb_i_i_m_nvp))) {
					ovb_chi = chisq.test(ovb_i_i_m_nvp)
					write.table(t(ovb_chi$observed), sprintf("%s/tests_chisq/test.%s.%s.%s.csv", out_fn, cmpnm, spi, chr), quote = FALSE)
					n_bins = ncol(ovb_i_i_m_nvp)
				} else {
					ovb_chi = list()
					ovb_chi$method = "Chisq fails"
					ovb_chi$p.value = 1
					n_bins = 1
				}
				# chisq, observed and expected: overlapping bings
				if (dim(ovb_i_i_m)[1]>1 & dim(ovb_i_i_m)[2]>1) {
					ovb_chi_ovs = chisq.test(ovb_i_i_m)
				} else {
					ovb_chi_ovs = NULL
				}
				
				# line plots
				if (dim(jac_i_i_m)[2] > 1) {
					jac_i_i_m_f = t(t(jac_i_i_m) / apply(t(jac_i_i_m), 1, sum))
					jac_i_i_m_f [ is.na(jac_i_i_m_f) ] = 0
					plot(x=NA, y=NA, xlab = "Mb", cex.names = 0.7, cex.main = 0.7, cex.axis = 0.7, ylab = "Jaccard", cex.lab = 0.7, main = sprintf("ALG contribution to %s", chr), xlim = c(0, max(met_i_midpoint)/1e6), ylim = c(0,1), las = 2)
					for (r in 1:nrow(jac_i_i_m_f)) {
						lines(x=met_i_midpoint[names(jac_i_i_m_f[r,])]/1e6, y=jac_i_i_m_f[r,], col = colors_row [ rownames(jac_i_i_m)[r] ], las = 2)
					}
					title(sub = sprintf("%s p = %.2E (%i non-overlapping bins)\np = %.2E (%i overlapping bins)", ovb_chi$method, ovb_chi$p.value, n_bins, ovb_chi_ovs$p.value, ncol(ovb_i_i_m)), cex.sub = 0.7)
					abline(v=met_i_midpoint[colnames(ovb_i_i_m_nvp)]/1e6, lty = 2, col = "gray80")
					legend("topright", fill = colors_row [ rownames(jac_i_i_m) ], legend = rownames(jac_i_i_m), bty = "n", cex = 0.5)
				} else {
					plot(0,0)
				}
				
				# plot
				if (!is.null(ovb_chi_ovs) & !is.null(dim(ovb_chi_ovs$expected))) {
					
					plot(x=NA, y=NA, xlab = "Mb", cex.names = 0.7, cex.main = 0.7, cex.axis = 0.7, ylab = "log2(O/E)", cex.lab = 0.7, main = sprintf("ALG contribution to %s", chr), xlim = c(0, max(met_i_midpoint)/1e6), ylim = c(0,4), las = 2)
					for (r in 1:nrow(jac_i_i_m_f)) {
						oe = log2(t(ovb_chi_ovs$observed)[,r] / t(ovb_chi_ovs$expected)[,r])
						oe [ oe < -0 ] =  0
						oe [ oe >  4 ] =  4
						oe [ is.na(oe) ] = 4
						lines(x=met_i_midpoint[names(jac_i_i_m_f[r,])]/1e6, y=oe, col = colors_row [ rownames(jac_i_i_m)[r] ], las = 2, pch = 19)
					}
					abline(v=met_i_midpoint[colnames(ovb_i_i_m_nvp)]/1e6, lty = 2, col = "gray80")

					
					
				} else {
					plot(0,0)
				}

			}
			
		}
		dev.off()

	}

	# all heatmaps side by side	
	pdf(sprintf("%s/alg.%s.general.pdf", out_fn, cmpnm), height = 6, width = 21)
	for (spi in sps_list) { 
		message(sprintf("ALG %s | plot heatmaps %s", cmpnm, spi))
		list_hms_chr[[spi]]@column_names_param$anno@height = unit(1.5, "cm")
		list_hms_jac[[spi]]@column_names_param$anno@height = unit(1.5, "cm")
		list_hms_ovb[[spi]]@column_names_param$anno@height = unit(1.5, "cm")
		list_hms_chr[[spi]]@left_annotation_param$width    = unit(0.6, "cm")
		list_hms_chr[[spi]]@heatmap_param$width  = unit(0.4 * ncol(list_hms_chr[[spi]]@matrix) + 0.6,  "cm")
		list_hms_jac[[spi]]@heatmap_param$width  = unit(18, "cm")
		list_hms_ovb[[spi]]@heatmap_param$width  = unit(18, "cm")
		list_hms_chr[[spi]]@heatmap_param$height = unit(12, "cm")
		list_hms_jac[[spi]]@heatmap_param$height = unit(12, "cm")
		list_hms_ovb[[spi]]@heatmap_param$height = unit(12, "cm")
		print(list_hms_chr[[spi]] + list_hms_jac[[spi]] + list_hms_ovb[[spi]])
	}
	dev.off()

}


message("all done!")
