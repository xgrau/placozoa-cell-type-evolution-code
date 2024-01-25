# libraries
library("scales")
library("metacell")
source("../scripts/helper.R")

# files
out_fn = "results_metacell_it4"

# list of species
sps_list = c("Tadh","TrH2","Hhon","HoiH23")

# loop
for (spi in sps_list) {
	
	# init database
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	
	# first load reordered recluster mc solution
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it4",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it2",run_name))
	
    # metacell all
    ctt_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
    ctt = read.table(ctt_fn, header = TRUE, sep = "\t", comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$broad_cell_type = factor(ctt$broad_cell_type, levels = unique(ctt$broad_cell_type))

	# footprints
	message(sprintf("%s | fp cts", spi))
	cts_fp = sca_cell_type_fp(ctt, mc, mat)
	message(sprintf("%s | fp bct", spi))
	bct_fp = sca_cell_type_fp(ctt[,c("metacell","broad_cell_type","color")], mc, mat)

	# umifrac
	message(sprintf("%s | umifracs", spi))
	cts_mc_v = ctt$cell_type
	bct_mc_v = ctt$broad_cell_type
	names(cts_mc_v) = ctt$metacell
	names(bct_mc_v) = ctt$metacell
	cts_sc = cts_mc_v [ mc@mc ]
	bct_sc = bct_mc_v [ mc@mc ]
	names(cts_sc) = names(mc@mc)
	names(bct_sc) = names(mc@mc)
	cts_counts = sca_mc_gene_counts_noobj(as.matrix(mat@mat), grouping_vector = cts_sc)
	bct_counts = sca_mc_gene_counts_noobj(as.matrix(mat@mat), grouping_vector = bct_sc)
	cts_umifra = sca_mc_gene_umifrac_noobj(cts_counts)
	bct_umifra = sca_mc_gene_umifrac_noobj(bct_counts)

	# second load peptidergic solution
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it4_pep_ord",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it4_pep",run_name))
	
    # metacell pep
    ctt_fn = sprintf("results_metacell_it4_peptidergic/annotation_mc.%s.it4.peptidergic.tsv", spi)
    ctt = read.table(ctt_fn, header = TRUE, sep = "\t", comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	
	# footprints
	message(sprintf("%s | fp cts pep", spi))
	pep_cts_fp = sca_cell_type_fp(ctt, mc, mat)

	# umifrac
	message(sprintf("%s | umifracs pep", spi))
	cts_mc_v = ctt$cell_type
	names(cts_mc_v) = ctt$metacell
	cts_sc = cts_mc_v [ mc@mc ]
	names(cts_sc) = names(mc@mc)
	cts_pep_counts = sca_mc_gene_counts_noobj(as.matrix(mat@mat), grouping_vector = cts_sc)
	cts_pep_umifra = sca_mc_gene_umifrac_noobj(cts_pep_counts)

	# save
	message(sprintf("%s | fp save", spi))
	metacell::scdb_add_mc(sprintf("%s_it4_cts",run_name), cts_fp)
	metacell::scdb_add_mc(sprintf("%s_it4_bct",run_name), bct_fp)
	metacell::scdb_add_mc(sprintf("%s_it4_pep_cts",run_name), pep_cts_fp)
	write.table(cts_counts, sprintf("results_metacell_it4/scdr_%s.matrix.cts_counts.csv", spi), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
	write.table(bct_counts, sprintf("results_metacell_it4/scdr_%s.matrix.bct_counts.csv", spi), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
	write.table(cts_pep_counts, sprintf("results_metacell_it4_peptidergic/scdr_%s.matrix.cts_counts.csv", spi), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
	write.table(cts_umifra, sprintf("results_metacell_it4/scdr_%s.matrix.cts_umifrac.csv", spi), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
	write.table(bct_umifra, sprintf("results_metacell_it4/scdr_%s.matrix.bct_umifrac.csv", spi), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
	write.table(cts_pep_umifra, sprintf("results_metacell_it4_peptidergic/scdr_%s.matrix.cts_umifrac.csv", spi), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
	

}




message("All done!")




