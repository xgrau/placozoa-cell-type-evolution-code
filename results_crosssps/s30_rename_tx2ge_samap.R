source("../scripts/helper.R")

list_blast_tables = list.files("data_samap/", pattern = ".csv.gz", full.names = TRUE)
list_blast_tables = list_blast_tables [ !grepl(".genes.csv.gz", list_blast_tables) ]

for (tab_fn in list_blast_tables) {

	# load table	
	tab = read.table(tab_fn)

	# find out sps to compare	
	sp1 = stringr::str_split(tab[1,1], "_", simplify = TRUE)[,1]
	sp2 = stringr::str_split(tab[1,2], "_", simplify = TRUE)[,1]

	# rename
	message(sprintf("rename %s: %s and %s species", tab_fn, sp1, sp2))
	tao = tab
	tao[,1] = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", sp1), vector_to_fix = tab[,1])
	tao[,2] = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", sp2), vector_to_fix = tab[,2])
	
	# write output
	tao_fn = gsub(".csv.gz$", ".genes.csv", tab_fn)
	message(sprintf("write %s: %s and %s species", tao_fn, sp1, sp2))
	write.table(tao, tao_fn, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
	R.utils::gzip(tao_fn, overwrite = TRUE)
	
	gc()
	
}