#### Input ####

# libraries
library("metacell")
library("scales")
source("../scripts/helper.R")
graphics.off()

## Placozoa ##

# spslist
sps_list = c("Tadh","TrH2","Hhon","HoiH23")
set_list = list(
	global = list(inp_fn = "results_gmod_it4/", out_fn = "results_gmod_it4c/", mc_sprintf_string = "%s_it4", mat_sprintf_string = "%s_it2", ctt_sprtinf_string = "../results_scatlas/results_metacell_it4/annotation_mc.%s.it4.reordered.tsv")
)

# loop
for (set in set_list) {

	out_fn = set$out_fn
	inp_fn = set$inp_fn
	mc_sprintf_string = set$mc_sprintf_string
	mat_sprintf_string = set$mat_sprintf_string
	ctt_sprtinf_string = set$ctt_sprtinf_string

	dir.create(out_fn, showWarnings = FALSE)

	# loop
	for (spi in sps_list) {
		
		# info
		message(sprintf("annotate modules %s", spi))
		run_name = sprintf("scdr_%s", spi)
		
		# load data
		metacell::scdb_init("../results_scatlas/data/scdb/",force_reinit=TRUE)
		mc  = metacell::scdb_mc(sprintf(mc_sprintf_string,run_name))
		# mat = metacell::scdb_mat(sprintf(mat_sprintf_string,run_name))
		
		# load cell type annotations for it4
		ctt_fn = sprintf(ctt_sprtinf_string, spi)
		ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")
		
		# define variable genes
		var_genes = read.table(sprintf("%s/gmod_%s.wgcna_object.var_genes.txt", inp_fn, spi))[,1]
		
		# read in wgcna modules
		mc_wgcna = readRDS(sprintf("%s/gmod_%s.wgcna_object.wgcna.rds", inp_fn, spi))
		mc_wgcna_gmods = readRDS(sprintf("%s/gmod_%s.wgcna_object.gmods.rds", inp_fn, spi))
		mc_wgcna_me = readRDS(sprintf("%s/gmod_%s.wgcna_object.ME.rds", inp_fn, spi))

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
		
		# list of tfs to highlight
		list_tfs = rownames(dic_genetf)
		list_tfs_annot = gene_annot[list_tfs, 1 ]
		colnames(gene_annot) = c("gene_name","pfam")
		
		# plot mc fp
		gmod_plotMCheatmap_annotate_modules(
			expr_matrix = mc@mc_fp,
			gmods = mc_wgcna_gmods,
			me = mc_wgcna_me,
			expr_matrix_colors = mc@colors,
			ex_output_file = sprintf("%s/gmod_%s.gmod_expression.pdf", out_fn, spi),
			an_output_file = sprintf("%s/gmod_%s.gmod_annotation.csv", out_fn, spi),
			me_output_file = sprintf("%s/gmod_%s.gmod_eigengenes.pdf", out_fn, spi),
			ex_width = 24, ex_height = 12,
			me_width = 8, me_height = 6,
			# resolution_rate = 1, 
			eigen_min = 0,
			cor_cutoff_max = NULL,
			annotation = gene_annot,
			highlight_genes = list_tfs, 
			heatmap_colors = c("gray99","#accbcc","#508490","#004066","#000738"),
			highlight_genes_annot = list_tfs_annot)
			
			
		# map modules to mcs and cell types
		ct_vector = ctt$cell_type
		names(ct_vector) = ctt$metacell
		gmod_annotate_modules_to_mc_and_ct(me = mc_wgcna_me, ct_vector = ct_vector, output_fn = sprintf("%s/gmod_%s.gmod_to_cts.csv", out_fn, spi))
		
		# map modules to mcs and broad cell types
		ct_vector = ctt$broad_cell_type
		if (!is.null(ct_vector)) {
			names(ct_vector) = ctt$metacell
			gmod_annotate_modules_to_mc_and_ct(me = mc_wgcna_me, ct_vector = ct_vector, output_fn = sprintf("%s/gmod_%s.gmod_to_bct.csv", out_fn, spi))
		}
			
	}

}
