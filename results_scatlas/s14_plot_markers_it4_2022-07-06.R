# libraries
library("metacell")
library("scales")
source("../scripts/helper.R")
par(family  = "Arial")

# where to store output?
out_fn = "results_metacell_it4/"

# list of species
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

for (spi in sps_list) {
	
	# init database
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	
	# first load old mc solution
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it4",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it2",run_name))
	
	# load cell type annotations for it4
	ctt_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.tsv", spi)
	ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")
	
	# get mc counts, umifrac
	mc_counts = sca_mc_gene_counts(mc,mat,0)

	# load gene annotations (based on transcripts, need to be changed to genes)		
	gene_annot = read.table(sprintf("../data/reference/%s_long.pep.annotations.csv", spi), sep = "\t", row.names = 1)
	rownames(gene_annot) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = rownames(gene_annot))
	# load TF-specific annotations
	dic_genetf = read.table(sprintf("../data/gene_annotations/tfs.%s_genes.curated.csv", spi), sep = "\t", row.names = 1, col.names = c("gene","OG"))
	dic_genetf_gene_name = gsub("like:","like_",dic_genetf$OG)
	dic_genetf_gene_name = stringr::str_split(dic_genetf_gene_name, pattern = ":", simplify = TRUE)[,2]
	rownames(dic_genetf) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = rownames(dic_genetf))
	names(dic_genetf_gene_name) = rownames(dic_genetf)
	gene_annot_tfs = merge(gene_annot, dic_genetf_gene_name, by.x = 0, by.y = 0, all.x = TRUE, all.y = FALSE)[,"y"]
	gene_annot [ !is.na(gene_annot_tfs) , "V2" ] = gene_annot_tfs [ !is.na(gene_annot_tfs) ]

	# plot tfs and other markers
	dir.create(sprintf("%s/markers/", out_fn), showWarnings = FALSE)
	for (gene_subset in c("ecm","tfs","sig","sigother","sigpkin","siggpcr","flag","mus","chr","neu","rbp","myo","ion","sterol","ppsyn","stereocilin","carbanhy","neuropept","neugaba","hypNPs")) {
		
		if ( file.exists(sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", gene_subset, spi)) ) {
			gene_subset_fn = sprintf("../data/gene_annotations/%s.%s_genes.curated.csv", gene_subset, spi)
			gene_subset_v = read.table(gene_subset_fn, sep = "\t", row.names = 1)
			rownames(gene_subset_v) = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = rownames(gene_subset_v))
		} else {
			gene_subset_fn = sprintf("../data/gene_annotations/%s.%s_genes.txt", gene_subset, spi)
			gene_subset_l = unique(read.table(gene_subset_fn, sep = "\t")[,1])
			gene_subset_l = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi), vector_to_fix = gene_subset_l)
			gene_subset_v = data.frame(row.names = gene_subset_l, annot = paste(gene_annot[gene_subset_l,1], stringr::str_trunc(gene_annot[gene_subset_l,2], width = 60), sep = " | "))
			# if statement to avoid having 1-row lists: will add a random gene
			if (nrow(gene_subset_v) == 1)  {
				gene_subset_v = rbind(gene_subset_v, data.frame(row.names = rownames(mc_counts) [ ! rownames(mc_counts) %in% rownames(gene_subset_v) ] [1], annot = "PADDING"))
			}
		}
		
		# custom list
		# custom_list = c("Tadh_TriadG64164","Tadh_TriadG64165","Tadh_TriadG63123","Tadh_TriadG59086","Tadh_TriadG63996","Tadh_TriadG55983","Tadh_TriadG63128","Tadh_TriadG37971","Tadh_TriadG63088","Tadh_TriadG64380","Tadh_TriadG63140","Tadh_TriadG35864")
		# custom_anns = custom_list
		# gene_subset = "custom"
		# gene_subset_v = data.frame(row.names = custom_list, annot = custom_anns)
		
		# plot
		gene_subset_m = scp_barplot_heatmap_markers(
			mc_object = mc,
			mat_object = mat,
			mc_counts = mc_counts,
			markers_file = gene_subset_v,
			heatmap_colors =  c("white","orange","orangered2","#520c52"),
			output_file_heatmap = sprintf("%s/markers/%s.markers_%s.heatmap.pdf", out_fn, run_name, gene_subset),
			output_file_barplot = sprintf("%s/markers/%s.markers_%s.barplot.pdf", out_fn, run_name, gene_subset),
			T_totumi = 10,
			width = 16,
			height = NULL,
			use_raster = FALSE,
			min_gene_fc = 1.5,
			min_expression_fc = 1,
			max_expression_fc = 3,
			mc_color = mc@colors,
			print_barplots = TRUE
		)
	}
	
	
}

message("All done!")

