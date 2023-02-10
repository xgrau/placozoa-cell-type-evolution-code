# Load pfam annotations as list, for functional enrichments
#' 
#' @param pfam_architecture_file tab-separated file where first column is gene name and second column is a string of domain architectures for that gene (by default, separated by whitespace)
#' @param do_remove_empty remove genes without annotation (default: TRUE)
#' @param do_unique remove redundant annotations for each gene (default: TRUE)
#' @param architecture_sep character separating the annotations in each row (by default, whitespace) (default: whitespace)
#' 
gsa_enrichment_load_pfam_list = function(pfam_architecture_file, do_unique = TRUE, do_remove_empty = TRUE, architecture_sep = " ") {
	
	# load pfam species i
	pfm_i = read.table(pfam_architecture_file, col.names = c("gene","architecture"), sep = "\t")
	if (do_unique) {
		pfm_i_l = lapply(pfm_i$architecture, function(v) unique(unlist(strsplit(v, split = architecture_sep))) )
		names(pfm_i_l) = pfm_i$gene
	} else {
		pfm_i_l = lapply(pfm_i$architecture, function(v) unlist(strsplit(v, split = architecture_sep)) )
		names(pfm_i_l) = pfm_i$gene
	}
	
	if (do_remove_empty) {
		pfm_i_l = pfm_i_l [ lengths(pfm_i_l) > 0 ]
	}
	
	# output
	return(pfm_i_l)
	
}

# Functional enrichment of gene lists with a hypergeometric test
#' 
#' @param annotation named list, where each entry is a gene, and each entry contains a list of annotations mapped to this gene (e.g. a unique list of pfam domains in the gene)
#' @param genes_fg vector of genes in the foreground (list of genes of interest)
#' @param genes_bg vector of genes in the background (this defaults to all genes in the annotation object, but you can specify a narrower background if needed)
#' @param do_print_output,name_fg,top_markers control output plot
#' @param p_adj_method pvalue adjustment method for p.adj function (default is BH)
#' 
#' @return data.frame with annotation enrichments
#'
gsa_enrichment_hypergeometric = function(
		annotation,
		genes_fg,
		genes_bg = NULL,
		do_print_output = TRUE,
		output_prefix = "enrichment_output",
		name_fg = "fg",
		top_markers = 30,
		p_adj_method = "BH") {
	
	if (is.null(genes_bg)) {
		genes_bg = names(annotation)
	}
	
	# subset list of annotations to genes in foreground and background
	annotation = annotation [ names(annotation) %in% c(genes_fg, genes_bg) ]
	
	# count distributions in background and foreground list
	dist_fg = as.data.frame(table(unlist(annotation [ genes_fg ])))
	dist_bg = as.data.frame(table(unlist(annotation)))
	
	if (nrow(dist_fg) > 0) {
		
		# prepare table
		tab_in_all  = merge(dist_fg, dist_bg, by.x = "Var1", by.y = "Var1", all.x = TRUE, all.y = FALSE)
		tab_in_all  = tab_in_all[apply(tab_in_all!=0, 1, all),]
		colnames(tab_in_all) = c("annot","freq_in_fg","freq_in_bg")
		tab_in_all$total_annot_in_fg = length(unlist(annotation [ genes_fg ]))
		tab_in_all$total_annot_in_bg = length(unlist(annotation))
		
		# hypergeometric test
		tab_in_all$pval = phyper(
			tab_in_all$freq_in_fg - 1,
			tab_in_all$freq_in_bg,
			tab_in_all$total_annot_in_bg - tab_in_all$freq_in_bg,
			tab_in_all$total_annot_in_fg,
			lower.tail = FALSE)
		
		# pvalue adjustment
		tab_in_all$pval_adj = p.adjust(tab_in_all$pval, method = p_adj_method)
		
		# sort result by pvalue
		tab_in_all = tab_in_all [ order(tab_in_all$pval_adj) , ]
		
		# Output 
		if(do_print_output) {
			
			pdf(paste(output_prefix,".",name_fg,".hypergeo.pdf",sep=""), height = 4 + length(top_markers / 10), width = 8)
			# par(mar=c(5,10,5,2))
			layout(matrix(1:3, ncol = 3))
			
			# padding plot			
			plot(NA,NA,xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n")
			
			# pval plot
			ploti = barplot(
				height = -rev(head(log10(tab_in_all$pval_adj), top_markers)),
				names.arg = rev(head(tab_in_all$annot, top_markers)),
				xlim = c(0,5), horiz = TRUE,las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
				main = name_fg,
				sub = sprintf("n=%i/%i input genes with annotations", nrow(dist_fg), length(genes_fg)),
				xlab = "-log10(p)"
			)
			abline(v=log10(0.01),lty=2,lwd=0.5,col="indianred2")
			abline(v=log10(0.05),lty=2,lwd=0.5,col="indianred2")
			text(
				x= 0, ploti,
				labels = sprintf("p=%.1E | n=%i", rev(head(tab_in_all$pval_adj, top_markers)) , rev(head(tab_in_all$freq_in_fg, top_markers))),
				col = "gray30", pos = 4)
			
			# fg/bg fraction
			ploti = barplot(
				height = rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)),
				xlim = c(0,1), horiz = TRUE, las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
				main = "fraction genes in fg and bg",
				xlab = "fraction"
			)
			ploti = barplot(
				height = rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers)),
				xlim = c(0,1), horiz = TRUE, las = 1,col = "indianred1", border = NA, ylim = c(0, top_markers + 10),
				main = "fraction genes in fg and bg",
				xlab = "fraction",
				add = TRUE, xaxt = "n"
			)
			text(
				x= 0, ploti,
				labels = sprintf("fg=%.2f | bg=%.2f", rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)) , rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers))),
				col = "gray30", pos = 4)
			# ploti = barplot(
			# 	height = rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)) / rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers)),
			# 	names.arg = rev(head(tab_in_all$annot, top_markers)),
			# 	horiz = TRUE, las = 1, col = "azure3", border = NA, ylim = c(0, top_markers + 10), log= "x",
			# 	xlim = c(1, 100),
			# 	main = "FC",
			# 	xlab = "FC",
			# )
			# text(
			# 	x=1, ploti,
			# 	labels = sprintf("fc=%.2f", rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)) / rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers))),
			# 	col = "gray30", pos = 4)
			
			
			dev.off()
			
			# table
			write.table(
				tab_in_all,
				file = sprintf("%s.%s.hypergeo.csv", output_prefix, name_fg),
				row.names = FALSE, sep="\t", quote = FALSE)
		}
		
	} else {
		print("skip, no annotations in interest list!")
		tab_in_all = NULL
	}
	return(tab_in_all)
	
}

### Topgo ###

gsa_topgo_load_emapper = function(
	emapper_fn,
	sep_col = "\t",
	sep_gos = ",",
	index_col_GOs = 10,
	index_col_gen = 1
) {
	
	gos_i = read.table(emapper_fn, sep = sep_col)
	gos_i = data.frame("gene" = gos_i[,index_col_gen], "GO" = gos_i[,index_col_GOs])
	gos_i_l = lapply(gos_i$GO, function(v) unlist(strsplit(v, split = ",")) )
	names(gos_i_l) = gos_i$gene
	return(gos_i_l)
	
}


gsa_topgo_enrichment  = function(
		annotation,
		genes_fg,
		genes_bg = NULL,
		do_print_output = TRUE,
		output_prefix = "enrichment_output",
		name_fg = "fg",
		ontologyset = c("BP","MF","CC"),
		tg_test = "fisher",
		tg_algorithm = "elim",
		top_markers = 30,
		nodesize = 10,
		printfile = TRUE,
		p_adj_method="BH") {
	
	library("topGO")
	
	# Input 
	genes_fg = unique(genes_fg)
	
	
	if (!is.null(genes_bg)) {
		genes_bg = genes_fg [ genes_fg %in% names(annotation) ]
	} else  {
		genes_bg = names(annotation)
	}
	genes_fg_ix = factor(as.integer(genes_bg %in% genes_fg))
	names(genes_fg_ix) = genes_bg
	
	# shortened go mappings without empty transcripts
	gomap_nonempty = annotation[lapply(annotation,length)>0]
	
	if (do_print_output) {
		
		pdf(paste(output_prefix,".", name_fg,".topgo.pdf", sep=""), height = (4 + length(top_markers / 10)) * 3, width = 8)
		# par(mar=c(5,10,5,2))
		layout(matrix(1:9, ncol = 3, byrow = TRUE))
		
	}
	
	# test each go domain
	topgo_tau_tot = data.frame()
	
	if (length(genes_fg[genes_fg %in% names(gomap_nonempty)])>1) {
		
		for (ontologyseti in ontologyset) {
			
			# topgo setup
			GOdata = suppressMessages(new("topGOdata", ontology = ontologyseti, allGenes = genes_fg_ix, annot = annFUN.gene2GO, gene2GO = annotation))
			num_interest_feasible = sum(GOdata@feasible & genes_bg %in% genes_fg)
			
			# topGO test 
			topgo_res = suppressMessages(topGO::runTest(GOdata, algorithm = tg_algorithm, statistic = tg_test))
			topgo_tau = suppressMessages(topGO::GenTable(GOdata, pval_test = topgo_res, orderBy = "pval_test",  topNodes = length(usedGO(object = GOdata))))
			topgo_tau$pval_test = as.numeric(topgo_tau$pval_test)
			topgo_tau$pval_adj  = p.adjust(topgo_tau$pval_test, method=p_adj_method)
			topgo_tau$ontology = ontologyseti
			topgo_tau_tot = rbind(topgo_tau_tot,topgo_tau)
			
			if (do_print_output) {
				
				# padding plot
				plot(NA,NA,xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n")
				
				# pval plot
				ploti = barplot(
					height = -rev(head(log(topgo_tau$pval_test,10), top_markers)),
					names.arg = rev(head(paste(topgo_tau$Term, topgo_tau$GO.ID), top_markers)),
					xlim = c(0,5), horiz = TRUE,las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
					main = sprintf("GO:%s\n%s", ontologyseti, name_fg),
					sub = sprintf("n=%i/%i input genes with annotations", num_interest_feasible, length(genes_fg)),
					xlab="-log(p)")
				abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
				abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
				text(
					x= 0, ploti,
					labels = sprintf("p=%.1E | n=%i", rev(head(topgo_tau$pval_test, top_markers)), rev(head(topgo_tau$Significant, top_markers))),
					col = "gray30", pos = 4)
				
				# fraction significant
				ploti = barplot(
					height = rev(head(topgo_tau$Significant / num_interest_feasible, top_markers)),
					xlim = c(0,1), horiz = TRUE, las = 1, col = "azure3", border = NA, ylim = c(0, top_markers + 10),
					main = "fraction genes in fg and expected value",
					xlab = "fraction"
				)
				ploti = barplot(
					height = rev(head(topgo_tau$Expected / num_interest_feasible, top_markers)),
					xlim = c(0,1), horiz = TRUE, las = 1, col = "indianred1", border = NA, ylim = c(0, top_markers + 10),
					add = TRUE, xaxt = "n"
				)
				text(
					x= 0, ploti,
					labels = sprintf("fg=%.2f | bg=%.2f", rev(head(topgo_tau$Significant / num_interest_feasible, top_markers)), rev(head(topgo_tau$Expected / num_interest_feasible, top_markers)) ),
					col = "gray30", pos = 4)
				
				
			}
			
		}
		
	} else {
		print("skip, no annotations in interest list!")
	}
	
	if (do_print_output) {
		# table
		write.table(
			topgo_tau_tot [ topgo_tau_tot$Significant > 0, ],
			file = sprintf("%s.%s.topgo.csv", output_prefix, name_fg),
			row.names = FALSE, sep="\t", quote = FALSE)
		dev.off()
	}
	
	# return
	return(topgo_tau_tot)
	
}



#### VENN DIAGRAMS ####
# plot venn diagrams from lists

venn.two = function(
		list1,
		list2,
		catname1,
		catname2,
		main="Venn",
		col1="green3",
		col2="magenta3",
		eulerbool=TRUE,
		print.mode=c("raw","percent")) {
	
	library(VennDiagram)
	library(grid)
	
	# unique
	list1 = unique(list1)
	list2 = unique(list2)

	# compute set overlaps, intersections, etc
	list_uniq      = base::unique(c(list1, list2))
	list_intersect = base::intersect(list1, list2)
	list_diff_i    = base::setdiff(list1, list2)
	list_diff_j    = base::setdiff(list2, list1)
	union = base::union(list1, list2)
	
	# draw venn
	venn = draw.pairwise.venn(
		area1 = length(list1),
		area2 = length(list2),
		cross.area = length(list_intersect),
		category=c(catname1,catname2),
		fontfamily="Helvetica",
		cat.fontfamily = "Helvetica",
		col=c(col1,col2),ext.text = F,
		cat.col=c(col1,col2),
		cat.cex  = 0.6, cex = 0.5,
		ind = F, print.mode = print.mode,
		cat.dist = 0.01,
		euler.d = eulerbool,
		scaled = eulerbool)
	grid::grid.newpage()
	grid::grid.draw(venn)
	grid::grid.text(sprintf("%s\nn=%i", main, length(union)),x=0.5,y=0.95, gp = gpar(fontsize = 5, fontface = "bold"))
	
	# output
	output = list(
		"catname_i"      = catname1,
		"catname_j"      = catname2,
		"title"          = main,
		"list_uniq"      = list_uniq,
		"list_intersect" = list_intersect,
		"list_diff_i"    = list_diff_i,
		"list_diff_j"    = list_diff_j,
		"venn" = venn,
		"union" = union
	)
	
	return(output)
	
}



venn.three = function(
		list1,
		list2,
		list3,
		catname1 = "a",
		catname2 = "b",
		catname3 = "c",
		main="Venn",
		col1="green3",
		col2="magenta3",
		col3="orange",
		eulerbool=TRUE,
		print.mode=c("raw","percent")) {
	
	library(VennDiagram)
	library(grid)
	
	# unique
	list1 = unique(list1)
	list2 = unique(list2)
	list3 = unique(list3)
	
	# compute set overlaps, intersections, etc
	intersect_n12  = base::intersect(list1, list2)
	intersect_n13  = base::intersect(list1, list3)
	intersect_n23  = base::intersect(list2, list3)
	intersect_n123 = base::intersect(intersect_n12, list3)
	union = base::union(union(list1, list2), list3)
	
	# draw venn
	venn = draw.triple.venn(
		area1 = length(list1),
		area2 = length(list2),
		area3 = length(list3),
		n12 = length(intersect_n12),
		n13 = length(intersect_n13),
		n23 = length(intersect_n23),
		n123 = length(intersect_n123),
		category=c(catname1,catname2,catname3),
		fontfamily="Helvetica",
		cat.fontfamily = "Helvetica",
		col=c(col1,col2,col3),ext.text = F,
		cat.col=c(col1,col2,col3),
		cat.cex  = 0.6, cex = 0.5,
		cat.dist = 0.01,
		ind = F, print.mode = print.mode,
		euler.d = eulerbool, scaled = eulerbool)
	grid::grid.newpage()
	grid::grid.draw(venn)
	grid::grid.text(sprintf("%s\nn=%i", main, length(union)),x=0.5,y=0.95, gp = gpar(fontsize = 5, fontface = "bold"))
	
	# output
	output = list(
		"venn" = venn,
		"catname_1"      = catname1,
		"catname_2"      = catname2,
		"catname_3"      = catname3,
		"intersect_12" = intersect_n12,
		"intersect_13" = intersect_n13,
		"intersect_23" = intersect_n23,
		"intersect_123" = intersect_n123,
		"union" = union
	)
	
	return(output)
	
}



venn.four = function(
		list1,
		list2,
		list3,
		list4,
		catname1 = "a",
		catname2 = "b",
		catname3 = "c",
		catname4 = "d",
		main="Venn",
		col1="green3",
		col2="magenta3",
		col3="orange",
		col4="blue3",
		eulerbool=TRUE,
		print.mode=c("raw","percent")) {
	
	library(VennDiagram)
	library(grid)

	# unique
	list1 = unique(list1)
	list2 = unique(list2)
	list3 = unique(list3)
	list4 = unique(list4)
	
	# compute set overlaps, intersections, etc
	intersect_n12  = base::intersect(list1, list2)
	intersect_n13  = base::intersect(list1, list3)
	intersect_n14  = base::intersect(list1, list4)
	intersect_n23  = base::intersect(list2, list3)
	intersect_n24  = base::intersect(list2, list4)
	intersect_n34  = base::intersect(list3, list4)
	intersect_n123 = base::intersect(intersect_n12, list3)
	intersect_n124 = base::intersect(intersect_n12, list4)
	intersect_n134 = base::intersect(intersect_n13, list4)
	intersect_n234 = base::intersect(intersect_n23, list4)
	intersect_nall = base::intersect(intersect_n123, list4)
	union = base::union(union(list1, list2), union(list3, list4))
	
	# draw venn
	venn = draw.quad.venn(
		area1 = length(list1),
		area2 = length(list2),
		area3 = length(list3),
		area4 = length(list4),
		n12 = length(intersect_n12),
		n13 = length(intersect_n13),
		n14 = length(intersect_n14),
		n23 = length(intersect_n23),
		n24 = length(intersect_n24),
		n34 = length(intersect_n34),
		n123 = length(intersect_n123),
		n124 = length(intersect_n124),
		n134 = length(intersect_n134),
		n234 = length(intersect_n234),
		n1234 = length(intersect_nall),
		category=c(catname1,catname2,catname3,catname4),
		fontfamily="Helvetica",
		cat.fontfamily = "Helvetica",
		col=c(col1,col2,col3,col4),ext.text = F,
		cat.col=c(col1,col2,col3,col4),
		cat.cex  = 0.6, cex = 0.5,
		cat.dist = 0.01,
		ind = F, print.mode = print.mode,
		euler.d = eulerbool, scaled = eulerbool)
	grid::grid.newpage()
	grid::grid.draw(venn)
	grid::grid.text(sprintf("%s\nn=%i", main, length(union)),x=0.5,y=0.95, gp = gpar(fontsize = 5, fontface = "bold"))
	
	# output
	output = list(
		"venn" = venn,
		"catname_1"      = catname1,
		"catname_2"      = catname2,
		"catname_3"      = catname3,
		"catname_4"      = catname4,
		"intersect_12" = intersect_n12,
		"intersect_13" = intersect_n13,
		"intersect_14" = intersect_n14,
		"intersect_23" = intersect_n23,
		"intersect_24" = intersect_n24,
		"intersect_34" = intersect_n34,
		"intersect_123" = intersect_n123,
		"intersect_124" = intersect_n124,
		"intersect_134" = intersect_n134,
		"intersect_234" = intersect_n234,
		"intersect_1234" = intersect_nall,
		"union" = union
		
	)
	
	return(output)
	
}


gsa_enrichment_scatter_plot = function(
	enrichment_table,
	annotation_col = "annot",
	main = "", 
	x_dim = "freq_in_fg", x_lab = "Num genes",
	y_dim = "pval_adj",   y_lab = "p-value",
	filter_y = NULL,
	filter_x = NULL,
	log = "y",
	color = "blue",
	y_lim = NULL,
	x_lim = NULL,
	categories_col = NULL) {
  
	enr = enrichment_table [ , c(annotation_col,x_dim, y_dim, categories_col) ] 
	
	if (!is.null(filter_y)) {
		enr[,annotation_col] [ enr[,y_dim] > filter_y ] = ""
	}
	if (!is.null(filter_x)) {
		enr[,annotation_col] [ enr[,x_dim] < filter_x ] = ""
	}
	
	if (is.null(y_lim)) {
		y_lim = c(0.1, min(enr[,y_dim]))
	}
	if (is.null(x_lim)) {
		x_lim = c(0.1, max(enr[,x_dim]))
	}
	
	if (is.null(categories_col)) {
		plot(enr[,x_dim], enr[,y_dim], pch = 19, col = color, log = "xy", xlim = x_lim, ylim = y_lim, las = 1, xlab = x_lab, ylab = y_lab)
		title(main = main)
		abline(v=filter_x, h=filter_y, lty = 2, col = "gray10")
		text(enr[,x_dim], enr[,y_dim], enr[,annotation_col], col = scales::alpha("thistle4",0.7), cex = 0.8)
	} else {
		for (ontology in levels(enr[,categories_col])) {
			eni = enr [ enr[,categories_col] == ontology,]
			if (nrow(eni)>0) {
				plot(eni[,x_dim], eni[,y_dim], pch = 19, col = color, log = log, xlim = x_lim, ylim = y_lim, las = 1, xlab = x_lab, ylab = y_lab)
				title(main = main, sub = ontology)
				abline(v=filter_x, h=filter_y, lty = 2, col = "gray10")
				text(eni[,x_dim], eni[,y_dim], eni[,annotation_col], col = scales::alpha("thistle4",0.7), cex = 0.8)
			}
		}
	}

}
