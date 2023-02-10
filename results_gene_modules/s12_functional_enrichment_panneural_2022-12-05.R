#### Input ####

# libraries
source("../scripts/helper.R")
source("../scripts/geneSetAnalysis.R")

# load expression evol data
top_sum = read.table("results_panneural_markers_broad/matrix.CB.tree_sum.csv", sep = "\t", header = TRUE) 
top_pre = read.table("results_panneural_markers_broad/anc.all.long.ogexpression.CB.dollo_pres.csv", header = TRUE, row.names = 1) 
top_gai = read.table("results_panneural_markers_broad/anc.all.long.ogexpression.CB.dollo_gain.csv", header = TRUE, row.names = 1) 
top_los = read.table("results_panneural_markers_broad/anc.all.long.ogexpression.CB.dollo_loss.csv", header = TRUE, row.names = 1) 

# using mus as referncefor Bilaterian-involving clades
# paths to input data
ref_gen = read.table("results_panneural_markers_broad/markers.neuron.Mmus.txt", sep = "\t", header = TRUE)
ref_gos = gsa_topgo_load_emapper("../data/reference/Mmus_ensembl.GO.csv", index_col_GOs = 2)
ref_pfa = gsa_enrichment_load_pfam_list("../data/reference/Mmus_long.pep.pfamscan_archs.csv")

# suppressMessages(metacell::scdb_init("../results_scatlas/data/scdb_outgroups/",force_reinit=TRUE))
# ref_nul = rownames(metacell::scdb_mc( sprintf("Mmus"))@mc_fp)
# ref_nul = as.character(dictionary_t2g(sprintf("../data/reference/Mmus_long.annot.gtf"), ref_nul, t2g = FALSE))
# ref_nul = ref_nul [ !is.na(ref_nul) ]

graphics.off()


dir.create("results_panneural_markers_broad/enrichments/", showWarnings = FALSE)

for (noi in c("CniBilPla","Bilateria","CniBil","Metazoa")) {
	
	# present
	message(sprintf("enrichments | %s, presence", noi))
	top_ogs_i = top_pre [ , noi ] > 0
	names(top_ogs_i)  = rownames(top_pre)
	top_ogs_i = names(which(top_ogs_i))
	top_gen_i = ref_gen [ ref_gen$orthogroup %in% top_ogs_i , "gene" ]
	enr_gos = gsa_topgo_enrichment(ref_gos, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out"), name_fg = sprintf("%s_pres", noi))
	enr_pfa = gsa_enrichment_hypergeometric(ref_pfa, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out"), name_fg = sprintf("%s_pres", noi))
	
	# joint plot
	if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
		enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
		colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
		enr_pfa_b$ontology = "Pfam"
		enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
		enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
		enr = enr [ enr$pval_test <= 0.1 , ]
		enr = enr [ enr$Significant > 0 , ]
		enr$Term = gsub(" ","_",enr$Term)
		enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
		pdf(sprintf("results_panneural_markers_broad/enrichments/scatter.%s.pdf", sprintf("%s_pres", noi)), width = 16, height = 16)
		layout(matrix(1:4, nrow = 2, byrow = TRUE))
		gsa_enrichment_scatter_plot(enr, main = sprintf("%s_pres", noi), annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
		dev.off()
	}
	
	# gained
	message(sprintf("enrichments | %s, gains", noi))
	top_ogs_i = top_gai [ , noi ] > 0
	names(top_ogs_i)  = rownames(top_gai)
	top_ogs_i = names(which(top_ogs_i))
	top_gen_i = ref_gen [ ref_gen$orthogroup %in% top_ogs_i , "gene" ]
	enr_gos = gsa_topgo_enrichment(ref_gos, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out"), name_fg = sprintf("%s_gain", noi))
	enr_pfa = gsa_enrichment_hypergeometric(ref_pfa, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out"), name_fg = sprintf("%s_gain", noi))
	
	# joint plot
	if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
		enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
		colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
		enr_pfa_b$ontology = "Pfam"
		enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
		enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
		enr = enr [ enr$pval_test <= 0.1 , ]
		enr = enr [ enr$Significant > 0 , ]
		enr$Term = gsub(" ","_",enr$Term)
		enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
		pdf(sprintf("results_panneural_markers_broad/enrichments/scatter.%s.pdf", sprintf("%s_gain", noi)), width = 16, height = 16)
		layout(matrix(1:4, nrow = 2, byrow = TRUE))
		gsa_enrichment_scatter_plot(enr, main = sprintf("%s_gain", noi), annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
		dev.off()
	}

	# lost
	message(sprintf("enrichments | %s, losses", noi))
	top_ogs_i = top_los [ , noi ] > 0
	names(top_ogs_i)  = rownames(top_los)
	top_ogs_i = names(which(top_ogs_i))
	top_gen_i = ref_gen [ ref_gen$orthogroup %in% top_ogs_i , "gene" ]
	enr_gos = gsa_topgo_enrichment(ref_gos, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out"), name_fg = sprintf("%s_loss", noi))
	enr_pfa = gsa_enrichment_hypergeometric(ref_pfa, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out"), name_fg = sprintf("%s_loss", noi))
	
	# joint plot
	if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
		enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
		colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
		enr_pfa_b$ontology = "Pfam"
		enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
		enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
		enr = enr [ enr$pval_test <= 0.1 , ]
		enr = enr [ enr$Significant > 0 , ]
		enr$Term = gsub(" ","_",enr$Term)
		enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
		pdf(sprintf("results_panneural_markers_broad/enrichments/scatter.%s.pdf", sprintf("%s_loss", noi)), width = 16, height = 16)
		layout(matrix(1:4, nrow = 2, byrow = TRUE))
		gsa_enrichment_scatter_plot(enr, main = sprintf("%s_loss", noi), annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
		dev.off()
	}

}



# using tadh as reference for placozoa, pfam
# paths to input data
ref_mus = read.table("results_panneural_markers_broad/markers.neuron.Mmus.txt", sep = "\t", header = TRUE)
top_gen = read.table("results_panneural_markers_broad/markers.peptidergic.Tadh.txt", sep = "\t", header = TRUE)
ref_gos = gsa_topgo_load_emapper("../data/reference/Tadh_ensembl.GO.csv", index_col_GOs = 2)
ref_pfa = gsa_enrichment_load_pfam_list("../data/reference/Tadh_long.pep.pfamscan_archs.csv")

dir.create("results_panneural_markers_broad/enrichments/")

for (noi in c("Placozoa")) {
	
	# present
	message(sprintf("enrichments | %s, presence", noi))
	top_ogs_i = top_pre [ , noi ] > 0
	names(top_ogs_i)  = rownames(top_pre)
	top_ogs_i = names(which(top_ogs_i))
	top_gen_i = top_gen [ top_gen$orthogroup %in% top_ogs_i , "gene" ]
	enr_pfa = gsa_enrichment_hypergeometric(ref_pfa, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out", noi), name_fg = sprintf("%s_pres", noi))
	enr_gos = gsa_topgo_enrichment(ref_gos, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out", noi), name_fg = sprintf("%s_pres", noi))
	
	# joint plot
	if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
		enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
		colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
		enr_pfa_b$ontology = "Pfam"
		enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
		enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
		enr = enr [ enr$pval_test <= 0.1 , ]
		enr = enr [ enr$Significant > 0 , ]
		enr$Term = gsub(" ","_",enr$Term)
		enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
		pdf(sprintf("results_panneural_markers_broad/enrichments/scatter.%s.pdf", sprintf("%s_pres", noi)), width = 16, height = 16)
		layout(matrix(1:4, nrow = 2, byrow = TRUE))
		gsa_enrichment_scatter_plot(enr, main = sprintf("%s_pres", noi), annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
		dev.off()
	}
	
	
	# gained
	message(sprintf("enrichments | %s, gains", noi))
	top_ogs_i = top_gai [ , noi ] > 0
	names(top_ogs_i)  = rownames(top_gai)
	top_ogs_i = names(which(top_ogs_i))
	top_gen_i = top_gen [ top_gen$orthogroup %in% top_ogs_i , "gene" ]
	enr_pfa = gsa_enrichment_hypergeometric(ref_pfa, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out", noi), name_fg = sprintf("%s_gain", noi))
	enr_gos = gsa_topgo_enrichment(ref_gos, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out", noi), name_fg = sprintf("%s_gain", noi))
	
	# joint plot
	if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
		enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
		colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
		enr_pfa_b$ontology = "Pfam"
		enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
		enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
		enr = enr [ enr$pval_test <= 0.1 , ]
		enr = enr [ enr$Significant > 0 , ]
		enr$Term = gsub(" ","_",enr$Term)
		enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
		pdf(sprintf("results_panneural_markers_broad/enrichments/scatter.%s.pdf", sprintf("%s_gain", noi)), width = 16, height = 16)
		layout(matrix(1:4, nrow = 2, byrow = TRUE))
		gsa_enrichment_scatter_plot(enr, main = sprintf("%s_gain", noi), annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
		dev.off()
	}

	
	# lost
	message(sprintf("enrichments | %s, losses", noi))
	top_ogs_i = top_los [ , noi ] > 0
	names(top_ogs_i)  = rownames(top_los)
	top_ogs_i = names(which(top_ogs_i))
	top_gen_i = top_gen [ top_gen$orthogroup %in% top_ogs_i , "gene" ]
	enr_pfa = gsa_enrichment_hypergeometric(ref_pfa, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out", noi), name_fg = sprintf("%s_loss", noi))
	enr_gos = gsa_topgo_enrichment(ref_gos, top_gen_i, output_prefix = sprintf("results_panneural_markers_broad/enrichments/out", noi), name_fg = sprintf("%s_loss", noi))

	# joint plot
	if (nrow(enr_gos) > 0 & !is.null(enr_pfa)) {
		enr_pfa_b = enr_pfa[,c("annot","freq_in_fg","pval")]
		colnames(enr_pfa_b) = c("Term","Significant","pval_test") 
		enr_pfa_b$ontology = "Pfam"
		enr = rbind(enr_gos[,c("Term","Significant","pval_test","ontology")], enr_pfa_b)
		enr$pval_test [ is.na(enr$pval_test) ] = min(enr$pval_test[!is.na(enr$pval_test)])
		enr = enr [ enr$pval_test <= 0.1 , ]
		enr = enr [ enr$Significant > 0 , ]
		enr$Term = gsub(" ","_",enr$Term)
		enr$ontology = factor(enr$ontology, levels = c("BP","MF","CC","Pfam"))
		pdf(sprintf("results_panneural_markers_broad/enrichments/scatter.%s.pdf", sprintf("%s_loss", noi)), width = 16, height = 16)
		layout(matrix(1:4, nrow = 2, byrow = TRUE))
		gsa_enrichment_scatter_plot(enr, main = sprintf("%s_loss", noi), annotation_col = "Term", x_dim = "Significant", y_dim = "pval_test", filter_y = 0.01, filter_x = 2, x_lim = c(1,max(enr$Significant)), categories_col = "ontology", log = "xy")
		dev.off()
	}
	
}

