# libraries
library("ape")
library("stringr")

# input
list_datasets = list(
  
  metc = list(
    sps_list = read.table("data/sps_list_metc.txt")[,1],
    tre_fn = "results_broccoli_metc/possvm_prefiltering/",
    ali_fn = "results_broccoli_metc/alignments_prefiltering/",
    out_fn = "results_broccoli_metc/alignments/",
    ous_fn = "results_broccoli_metc/statistics_per_OG.filtered.csv"
  ),
  meto = list(
    sps_list = read.table("data/sps_list_meto.txt")[,1],
    tre_fn = "results_broccoli_meto/possvm_prefiltering/",
    ali_fn = "results_broccoli_meto/alignments_prefiltering/",
    out_fn = "results_broccoli_meto/alignments/",
    ous_fn = "results_broccoli_meto/statistics_per_OG.filtered.csv"
  )
  
)

# probability threshold for closest homolog distribution
# for each gene, we retrieve its k nearest neighbours in the current phylogeny
# then, we check how often this species has the k species as its closest neighbours in the whole tree collection
# if this happens infrequently (< p_knn), ignore the gene, it's an outlier!
f_knn = 0.01
n_knn = 1


for (dataset in list_datasets) {
  
  # input
  sps_list = dataset[["sps_list"]]
  tre_fn = dataset[["tre_fn"]]
  ali_fn = dataset[["ali_fn"]]
  out_fn = dataset[["out_fn"]]
  ous_fn = dataset[["ous_fn"]]
  
  # list files
  ogs_fl = list.files(path = tre_fn, pattern = "*.groups.csv$", full.names = TRUE)
  tre_fl = list.files(path = tre_fn, pattern = "*.newick$", full.names = TRUE)
  ali_fl = list.files(path = ali_fn, pattern = "OG_[0-9]*.fasta$", full.names = TRUE)
  
  # work only with ogs where all data is available
  tre_done = gsub("\\..*","",basename(tre_fl))
  ali_done = gsub("\\..*","",basename(ali_fl))
  ali_fl = ali_fl [ ali_done %in% tre_done ]
  
  # distribution of cross-species distance expectations
  exm_m = matrix(ncol = length(sps_list), nrow = length(sps_list) * length(ogs_fl))
  for (i in 1:length(ogs_fl)) {
    
    if ((i - 1) %% 100 == 0 | i == length(ogs_fl) | i == 1) {
      message(sprintf("collect distribution of distances to nearest neighbours, %i trees processed...", i - 1))
    }
    
    # read classification and tree
    ogi = gsub("\\..*","",basename(ogs_fl[i]))
    tre = ape::read.tree(tre_fl[i])
    
    # clean tree and get pairwise distances
    tre$tip.label = gsub("\\|OG\\d+\\|$","",tre$tip.label)
    tre_d = ape::cophenetic.phylo(tre)
    sps_tips = stringr::str_split(rownames(tre_d), "_.*", simplify = TRUE)[,1]
    tre_d_s = tre_d [ !duplicated(sps_tips), !duplicated(sps_tips) ]
    rownames(tre_d_s) = sps_tips [ !duplicated(sps_tips) ]
    colnames(tre_d_s) = rownames(tre_d_s)
    
    # store in matrix
    exm_i = matrix(nrow = length(sps_list), ncol = length(sps_list))
    colnames(exm_i) = sps_list
    rownames(exm_i) = sps_list
    exm_i [ rownames(tre_d_s), colnames(tre_d_s) ] = tre_d_s
    diag(exm_i) = NA
    
    # save
    coord_s = 1 + (i - 1) * length(sps_list)
    coord_e = coord_s + length(sps_list) - 1
    exm_m [ coord_s : coord_e , ] = exm_i
    
  }
  
  
  exm_m_d = as.data.frame(exm_m)
  colnames(exm_m_d) = sps_list
  exm_m_d$species = rep(sps_list, length(ogs_fl) )

  # format as list for faster reference
  exm_l  = vector(mode = "list", length = length(sps_list))
  exm_lt = vector(mode = "list", length = length(sps_list))
  names(exm_l)  = sps_list
  names(exm_lt) = sps_list
  for (spi in sps_list) {
    exm_l_spi = exm_m_d [ exm_m_d$species == spi, colnames(exm_m_d) != "species" ]
    exm_l[[spi]]  = exm_l_spi
    exm_lt[[spi]] = sort(table( unlist( apply(exm_l_spi,1, function(v) names(head(sort(v), n_knn))) )), decreasing = TRUE)
  }
  
  
  
  # loop per marker
  stats_filt = data.frame()
  for (i in 1:length(ogs_fl)) {
    
    # read classification and tree
    ogi = gsub("\\..*","",basename(ogs_fl[i]))
    ogt = read.table(ogs_fl[i], sep = "\t", header = TRUE)
    tre = ape::read.tree(tre_fl[i])
    
    # read alignment
    ali = Biostrings::readAAStringSet(ali_fl[i], format = "fasta")
    ald = data.frame(row.names = names(ali), length = lengths(ali))
    
    # most abundant og (sps-wise)
    ogt$species = stringr::str_split(ogt$gene, "_", simplify = TRUE)[,1]
    ogt$species = factor(ogt$species, levels = sps_list)
    largest_og = names(sort(rowSums(table(ogt$orthogroup, ogt$species) > 0), decreasing = TRUE)[1])
    ogt$is_largest_og = ogt$orthogroup == largest_og
    message(sprintf("%s | %i gene(s) are out-paralogs", ogi, sum(!ogt$is_largest_og)))
    
    # clean tree and get pairwise distances
    tre$tip.label = gsub("\\|OG\\d+\\|$","",tre$tip.label)
    tre_d = ape::cophenetic.phylo(tre)
    tre_h = ape::node.depth.edgelength(tre) [ 1:length(tre$tip.label)  ]
    names(tre_h) = tre$tip.label
    
    # is duplicated?
    sps_with_dups = names(which(table(ogt$species [ ogt$is_largest_og ]) > 1))
    ogt$is_duplicated_in_og = ogt$species %in% sps_with_dups & ogt$is_largest_og
    message(sprintf("%s | %i gene(s) are in-paralogs", ogi, sum(ogt$is_duplicated_in_og)))
    
    # decision metrics:
    # distance to root
    ogt$dist_to_root = tre_h [ ogt$gene ]
    
    # length
    ogt$prot_length = ald [ ogt$gene, ]
    
    
    # which are the n nearest neighbour species of each gene?
    closest_neighbours = lapply(ogt$gene, function(g) { 
      tre_di = sort(tre_d [ g , ])
      sps_t = stringr::str_split(names(tre_di), "_", simplify = TRUE)[,1]
      sps_g = stringr::str_split(g, "_", simplify = TRUE)[,1]
      tre_di = tre_di [ names(tre_di) != g & sps_t != sps_g & !duplicated(sps_t) ]
      ge_closest = names(head(tre_di, n_knn))
      sp_closest = stringr::str_split(ge_closest, "_", simplify = TRUE)[,1]
    })
    names(closest_neighbours) = ogt$gene
    ogt$closest_neighbours = sapply(1:length(closest_neighbours), function(m) paste(sort(closest_neighbours[[m]]), collapse = ","))
    
    # how common it is to find these nearest neighbours in the expectation matrix?
    ogt$closest_neighbours_freq = sapply(1:length(ogt$species), function(n) { 
      spi = as.character(ogt[n,"species"])
      exm_i_a = exm_lt[[spi]]
      exm_s = sum(exm_i_a [ closest_neighbours[[n]] ])
      exm_s [ is.na(exm_s) ] = 0
      exm_p = exm_s / sum(exm_i_a)
      return(exm_p)
    })
    ogt$is_outlier_to_closest = ogt$closest_neighbours_freq < f_knn
    message(sprintf("%s | %i outliers gene(s) based on most frequent nearest neighbours", ogi, sum(ogt$is_outlier_to_closest & ogt$is_largest_og)))
    
    # which sequences we keep:
    ogt_o = ogt [ ogt$is_largest_og & ! ogt$is_outlier_to_closest , ]
    ogt_o = ogt_o [ order(ogt_o$species) , ]
    rownames(ogt_o) = ogt_o$gene
    
    # if there's any species with paralogs that are not outliers...
    sps_with_dups = names(which(table(ogt_o$species) > 1)) 
    if (length(sps_with_dups) > 0) {
      for (spi in sps_with_dups) {
        
        cgi = ogt_o [ ogt_o$species == spi , ] 
        cgi = cgi [ order(-cgi$prot_length, cgi$dist_to_root),]
        cgi_drop = rownames(cgi) [ duplicated(cgi$species) ]
        ogt_o = ogt_o [ ! rownames(ogt_o) %in% cgi_drop, ]
        
      }
    }
    
    # save
    message(sprintf("%s | %i out of %i genes kept (%i out of %i species)", ogi, nrow(ogt_o), nrow(ogt), length(unique(ogt_o$species)), length(unique(ogt$species)) ) )
    write.table(ogt, sprintf("%s/%s.keep.txt", out_fn, ogi), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    Biostrings::writeXStringSet(ali[ogt_o$gene], sprintf("%s/%s.fasta", out_fn, ogi))
    
    # keep stats
    stats_filt = rbind(
      stats_filt,
      data.frame(
        orthogroup = ogi,
        genes_input = nrow(ogt),
        genes_output = nrow(ogt_o),
        species_input = length(unique(ogt$species)),
        species_output = length(unique(ogt_o$species)),
        occupancy_output = length(unique(ogt_o$species)) / length(sps_list),
        kept_orthogroup = ogt_o$orthogroup[1]
      ))
    
  }
  
  write.table(stats_filt, ous_fn, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  message(sprintf("%i gene trees processed!", length(ogs_fl)))
  message("done!")
  
}
