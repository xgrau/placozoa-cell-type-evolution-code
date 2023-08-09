# out
out_fn = "results_macrosynteny_all"

# define species lists
# sps_list = read.table("species_list_synteny_blocks3.txt")[,1]
sps_list = c("Choren","Astrub","Rhoesc")

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

# only report ancestral linkage groups with at least that many HGs
min_hg_per_alg = 10

# function
table_to_matrix = function(table) {
  mat = matrix(table, nrow = nrow(table))
  rownames(mat) = rownames(table)
  colnames(mat) = colnames(table)
  return(mat)
}


# homology groups
message("read homology groups")
hom_d = read.table(sprintf("%s/mcl.all.out.txt", out_fn), header = FALSE, col.names = c("homology_group", "transcript"))
hom_d$species = gsub("_.*","", hom_d$transcript)
hom_d$species [ hom_d$species == "Hvul" ] = "Hvul_v3"
hom_d$species [ hom_d$species == "Nvec" ] = "Nvec_vc1.1"
hom_d = hom_d [ !duplicated(hom_d$transcript), ]
rownames(hom_d) = hom_d$transcript

# homology groups as factors
hom_d$homology_group = factor(hom_d$homology_group)

for (n in 1:length(list_comparisons)) {
	
	sps_list = list_comparisons[[n]]
	cmpnm = names(list_comparisons)[n]
	hom_i = hom_d [ hom_d$species %in% sps_list , ]


	# add chromosome of origin of each gene
	for (spi in sps_list) {
		
		message(sprintf("process %s | %s", cmpnm, spi))
		gen_i = rtracklayer::readGFFAsGRanges(sprintf("../data/reference/%s_long.annot.gtf", spi))
		gen_i = gen_i [ gen_i$type == "transcript" ]
		gen_d = as.data.frame(gen_i)
		gen_d = gen_d [ , c("seqnames","transcript_id") ]
		gen_d$homology_group = hom_i [ gen_d$transcript_id, "homology_group" ]
		
		# filter genes without homology group
		gen_d = gen_d [ !is.na(gen_d$homology_group), ]
		
		# filter rare chromosomes
		gen_d = gen_d [ gen_d$seqnames %in% names(which(table(gen_d$seqnames) > 50)), ]
		
		chh_i = gen_d [ , c("homology_group","seqnames") ]
		chh_i = chh_i [ !duplicated(chh_i),  ]
		colnames(chh_i)[2] = "chrom"
		chh_i$chrom = paste(spi, chh_i$chrom, sep = ":")
		
		if (spi == sps_list[1]) {
			chh_t = chh_i
		} else {
			chh_t = merge(chh_t, chh_i, by = "homology_group", all.x = FALSE, all.y = FALSE)
		}
		
	}

	# find common chromosome combinations
	message(sprintf("process %s | find common chromosome combinations", cmpnm))
	chh_t$chrom_combination = sprintf("%s %s %s", chh_t[,2], chh_t[,3], chh_t[,4])

	xta = table(chh_t$homology_group, chh_t$chrom_combination)
	mta = table_to_matrix(xta)

	# keep only homology group - chrom combination pairs that fulfill criteria:
	# - this HG is only associated to one chromosome triplet
	# - chromosome triplet has at least 10 HGs associated to it
	mta_f = mta
	mta_f = mta_f [ rowSums(mta_f) == 1 , ]
	mta_f = mta_f [ , apply(mta_f, 2, function(c) sum(c) >= min_hg_per_alg) ]
	mta_f = mta_f [ rowSums(mta_f) == 1 , ]

	# table of ancestral linkage groups
	alg_d = data.frame(
		homology_group =          rownames(mta_f),
		ancestral_linkage_group = sprintf("ALG%03d", apply(mta_f, 1, which.max)),
		chromosomes_in_ALG      = sapply(1:nrow(mta_f), function(r) { names(which.max(mta_f[r,])) })
	)
	alg_d = alg_d [ order(alg_d$ancestral_linkage_group) , ]

	# save
	message(sprintf("process %s | save ALGs", cmpnm))
	write.table(alg_d, sprintf("%s/alg.%s.csv", out_fn, cmpnm), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

	# plot counts
	pdf(sprintf("%s/alg.%s.sizes.pdf", out_fn, cmpnm), height = 4, width = 6)
	barplot(sort(table(alg_d$ancestral_linkage_group)), las = 2, col = "lightblue", cex = 0.7, cex.names = 0.7, cex.main = 0.7, cex.axis = 0.7, main = sprintf("homology groups per ancestral linkage group\nn=%i", length(table(alg_d$ancestral_linkage_group))), ylab = "num HGs")
	abline(h = min_hg_per_alg, lty = 2)
	dev.off()

}