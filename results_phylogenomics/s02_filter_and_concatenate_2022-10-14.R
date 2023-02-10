# libraries
library("ape")
library("pheatmap")

# input
list_datasets = list(

	cho09 = list(
		ali_fn = "results_broccoli_metc/alignments_trimmed/",
		p4f_fn = "results_broccoli_metc/alignments/",
		ogs_fn = "results_broccoli_metc/statistics_per_OG.filtered.it2.csv",
		sps_fn = "data/sps_list_metc.txt",
		tax_fn = "data/sps_list_metc.taxonomy.csv",
		out_fn = "results_broccoli_metc/mat.cho09",
		clades_interest = c("Bilateria","Cnidaria","Placozoa","Porifera","Ctenophora","Choanoflagellata"),
		mono_clades =  list(c("Bilateria","Cnidaria","Placozoa")),
		global_completeness_fraction = 0.7,
		chom_p = 0.01
	),
	
	
	meto10 = list(
		ali_fn = "results_broccoli_meto/alignments_trimmed/",
		p4f_fn = "results_broccoli_meto/alignments/",
		ogs_fn = "results_broccoli_meto/statistics_per_OG.filtered.it2.csv",
		sps_fn = "data/sps_list_meto.txt",
		tax_fn = "data/sps_list_meto.taxonomy.csv",
		out_fn = "results_broccoli_meto/meto10",
		clades_interest = c("Bilateria","Cnidaria","Placozoa","Porifera","Ctenophora"),
		mono_clades =  list(c("Bilateria","Cnidaria","Placozoa")),
		global_completeness_fraction = 0.7,
		chom_p = 0.01
	)
	
	
)


# loop over datasets
for (dataset in list_datasets) {
	
	# get dataset info
	ali_fn = dataset[["ali_fn"]]
	p4f_fn = dataset[["p4f_fn"]]
	ogs_fn = dataset[["ogs_fn"]]
	sps_fn = dataset[["sps_fn"]]
	tax_fn = dataset[["tax_fn"]]
	out_fn = dataset[["out_fn"]]
	chom_p = dataset[["chom_p"]]
	clades_interest = dataset[["clades_interest"]]
	mono_clades     = dataset[["mono_clades"]]
	global_completeness_fraction = dataset[["global_completeness_fraction"]]
	
	# start
	ali_list = list.files(path = ali_fn, pattern = "*ltt.fasta$", full.names = TRUE)
	ogs = read.table(ogs_fn, header = TRUE)
	rownames(ogs) = ogs$orthogroup
	ogs_list = ogs$orthogroup
	sps_list = as.character(read.table(sps_fn)[,1])
	
	# load taxonomy
	tax = read.table(tax_fn)
	colnames(tax) = c("species","clade")
	tax$clade = factor(tax$clade, levels = unique(tax$clade))
	tax_v = tax$clade
	names(tax_v) = tax$species
	
	# keep only aligned orthogroups
	ogs_al = gsub("\\..*","",basename(ali_list))
	ogs_list = ogs_list [ ogs_list %in% ogs_al ]
	
	# keep only species from clades of interest
	sps_list_f = sps_list [ tax_v %in% clades_interest ]
	
	# prepare output data
	out = matrix(data = NA, nrow = length(sps_list_f), ncol = length(ogs_list))
	rownames(out) = sps_list_f
	colnames(out) = ogs_list
	
	
	#### Whole dataset ####
	
	for (ogi in ogs_list) {
		
		# file to load  
		fai_fn = sprintf("%s/%s.ltt.fasta", ali_fn, ogi)
		tre_fn = sprintf("%s/%s.iqtree.treefile", p4f_fn, ogi)
		
		# if file exists, load sequences
		if (file.exists(fai_fn)) {
			
			# read data
			# message(sprintf("%s\t| Load alignment file: %s", ogi, fai_fn))
			fai = ape::read.FASTA(file = fai_fn, type = "AA")
			names(fai) = stringr::str_split(names(fai), pattern = "_", simplify = TRUE)[,1]
			fai_f = fai [ names(fai) %in% sps_list_f ]

			# check if global completeness is sufficient
			check_global_completeness = length(fai_f) >= length(sps_list_f) * global_completeness_fraction
			if (!check_global_completeness) {
				message(sprintf("%s\t| Drop marker: too many missing species (%i out of %i)", ogi, length(fai_f), length(sps_list_f) ))
			}
			ogs[ogi,"occupancy_in_interest_clades"] = length(fai_f) / length(sps_list_f)
			ogs[ogi,"n_species_in_interest_clades"] = length(fai_f)
			ogs[ogi,"check_global_completeness"] = check_global_completeness

			# check if alignment is long enough			
			seq_length = length(as.character(fai[[1]]))
			check_seq_length = seq_length >= 50
			if (!check_seq_length) { 
				message(sprintf("%s\t| Drop marker: too short (%i aa)", ogi, seq_length))
			}
			ogs[ogi,"check_seq_length"] = check_seq_length

			# check if there's enough coverage in the clades of interest
			tax_v_i = table(tax_v [ names(fai_f) ])
			check_taxon_wise_completeness = all(tax_v_i [ clades_interest ] > 3)
			if (!check_taxon_wise_completeness) {
				message(sprintf("%s\t| Drop marker: taxon-wise incompleteness", ogi))
			}
			ogs[ogi,"check_taxon_wise_completeness"] = check_taxon_wise_completeness
			
			# check if key clades are monophyletic
			if (!is.null(mono_clades)) {
				tre = ape::read.tree(tre_fn)
				tre$tip.label = gsub("_.*","", tre$tip.label)
				# tre = phangorn::midpoint(tre)
				key_taxon_monophyly_v = c()
				for (keyc in mono_clades) {
					spout = tax$species [ !tax$clade %in% keyc & tax$species %in%  tre$tip.label ][1]
					if (length(spout) > 1) {
						tre = ape::root(tre, spout)
					} else {
						tre = phangorn::midpoint(tre)	
					}
					key_taxa = tax [ tax$clade %in% keyc, "species" ]
					key_taxa = key_taxa [ key_taxa %in% tre$tip.label ]
					key_taxa_is_mono = ape::is.monophyletic(tre , key_taxa)
					key_taxon_monophyly_v = c(key_taxon_monophyly_v, key_taxa_is_mono)
				}
				if (sum(key_taxon_monophyly_v) >= length(mono_clades)) {
					check_key_taxon_monophyly = TRUE
				} else {
					check_key_taxon_monophyly = FALSE
				}
				if (!check_key_taxon_monophyly) {
					message(sprintf("%s\t| Drop marker: taxon-wise non-monophyly", ogi))
				}
			} else {
				check_key_taxon_monophyly = TRUE
			}

			# keep?
			keep = all(check_global_completeness, check_seq_length, check_taxon_wise_completeness, check_key_taxon_monophyly)

			# if keep...
			if (keep) {
				
				# get sequences from list object
				message(sprintf("%s\t| Keep marker", ogi))
				fai_all_sps = lapply(sps_list_f, function(s) {
					if (s %in% names(fai_f)) {
						paste(as.character(fai_f[names(fai_f) == s])[[1]], collapse = "")
					} else {
						paste(rep("-", seq_length), collapse = "")
					}
				})
				names(fai_all_sps) = sps_list_f
				fai_all_sps_mat = unlist(as.matrix(fai_all_sps)[,1])
				
				# concatenate
				out[,ogi] = fai_all_sps_mat
				
			} else {
				
				# drop column
				out = out[, - which(colnames(out) == ogi)]
				ogs_list = ogs_list [ ogs_list != ogi ]
				
			}
			
		} else {
			message(sprintf("%s | No alignment file: %s doesn't exist", ogi, fai_fn))
		}
		
	}
	
	
	# save dataset: as table
	write.table(out, file = sprintf("%s.all.tsv", out_fn), sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
	# save dataset: as fasta
	out_sequences = apply(out, 1, function(x) paste(x[!is.na(x)], collapse = ""))
	output_file = file(sprintf("%s.all.fasta", out_fn), open="w")
	for (i in 1:length(sps_list_f)) {
		
		writeLines(sprintf(">%s", sps_list_f[i]), con = output_file)
		writeLines(out_sequences[i], con = output_file)
		
	}
	close(output_file)
	
	# save as charset file for MARE
	output_file = file(sprintf("%s.all.charset.txt", out_fn), open="w")
	p_init = 1
	for (i in 1:ncol(out)) {
		
		writeLines(sprintf("charset %s = %i - %i ;", colnames(out)[i], p_init, p_init + nchar(out[1,i]) - 1), con = output_file)
		p_init = p_init + nchar(out[1,i])
		
	}
	close(output_file)
	
	
	### Topology consistency test ####
	
	
	### Compositional test ####
	
	# load results from p4 compositional heterogeneity tests
	# p<0.05 means not compositionally heterogeneous, ie homogeneous (PASS dataset)
	chom_d = rep(NA, length = length(ogs$orthogroup))
	names(chom_d) = ogs$orthogroup
	for (ogi in ogs$orthogroup) {
		if (file.exists(sprintf("%s/%s.p4.txt", p4f_fn, ogi))) {
			chom_d[ogi] = read.table(sprintf("%s/%s.p4.txt", p4f_fn, ogi))[,2]
		} else {
			chom_d[ogi] = 0
		}
	}
	chom_d [ is.na(chom_d) ] = 0
	
	# is this marker compositionally homogeneous?
	chom_d_pass = chom_d >= chom_p
	chom_d_fail = chom_d <  chom_p
	
	# save compositionally homogeneous pass
	# save dataset: as table
	out_f = out [ , chom_d_pass [ ogs_list ] ]
	write.table(out_f, file = sprintf("%s.comhom.tsv", out_fn), sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
	# save dataset: as fasta
	out_sequences = apply(out_f, 1, function(x) paste(x[!is.na(x)], collapse = ""))
	output_file = file(sprintf("%s.comhom.fasta", out_fn), open="w")
	for (i in 1:length(sps_list_f)) {
		
		writeLines(sprintf(">%s", sps_list_f[i]), con = output_file)
		writeLines(out_sequences[i], con = output_file)
		
	}
	close(output_file)
	
	# save as charset file for MARE
	output_file = file(sprintf("%s.comhom.charset.txt", out_fn), open="w")
	p_init = 1
	for (i in 1:ncol(out_f)) {
		
		writeLines(sprintf("charset %s = %i - %i ;", colnames(out_f)[i], p_init, p_init + nchar(out_f[1,i]) - 1), con = output_file)
		p_init = p_init + nchar(out_f[1,i])
		
	}
	close(output_file)
	
	
	# save compositionally homogeneous fail
	# save dataset: as table
	out_f = out[ , chom_d_fail [ ogs_list ] ]
	write.table(out_f, file = sprintf("%s.comhet.tsv", out_fn), sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
	# save dataset: as fasta
	out_sequences = apply(out_f, 1, function(x) paste(x[!is.na(x)], collapse = ""))
	output_file = file(sprintf("%s.comhet.fasta", out_fn), open="w")
	for (i in 1:length(sps_list_f)) {
		
		writeLines(sprintf(">%s", sps_list_f[i]), con = output_file)
		writeLines(out_sequences[i], con = output_file)
		
	}
	close(output_file)
	
	# save as charset file for MARE
	output_file = file(sprintf("%s.comhet.charset.txt", out_fn), open="w")
	p_init = 1
	for (i in 1:ncol(out_f)) {
		
		writeLines(sprintf("charset %s = %i - %i ;", colnames(out_f)[i], p_init, p_init + nchar(out_f[1,i]) - 1), con = output_file)
		p_init = p_init + nchar(out_f[1,i])
		
	}
	close(output_file)
	
	
	#### Filters summary ####
	
	# statistics
	ogs$compositional_homogeneity_test      = chom_d [ ogs$orthogroup ]
	ogs$compositional_homogeneity_test_pass = chom_d [ ogs$orthogroup ] >= chom_p
	
	# write table
	write.table(ogs, sprintf("%s.stats.tsv", out_fn), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
}