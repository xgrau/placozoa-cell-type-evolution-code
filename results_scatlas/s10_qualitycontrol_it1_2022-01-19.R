# libraries
library("metacell")
source("../scripts/helper.R")
par(family  = "Arial")

# where to store output?
out_fn = "results_metacell_it1/"

# load dual classification of cells
dou_fl = list.files(pattern = "dual_samples.*class.tsv", path = "data/", full.names = TRUE)
dou = data.frame()
for (dou_fn in dou_fl) {
    doi = read.table(dou_fn, header = TRUE)
    doi$dataset = gsub(".class.tsv", "", gsub("dual_samples.","", basename(dou_fn)))
    dou = rbind(dou, doi)
}


# list of species
sps_list = c("Tadh","TrH2","Hhon","HoiH23")
# sps_list = c("TrH2")

for (spi in sps_list) {
	
    # metacell and single cell annots
    mct_fn = sprintf("results_metacell_it1/annotation_mc.%s.it1.tsv", spi)
    sca_fn = sprintf("results_metacell_it1/scdr_%s.matrix.sc_annot.csv", spi)
    # load
    mct = read.table(mct_fn, header = TRUE, sep = "\t", comment.char = "")
    sca = read.table(sca_fn, header = TRUE, sep = "\t", comment.char = "")

    # merge with doublet data
    sca_m = merge(sca, dou, by.x = "cell", by.y = "long_cell_name", all.x = TRUE, all.y = FALSE)
    sca_m$metacell = factor(sca_m$metacell, levels = mct$metacell)
    
    # crosstabulate metacells and clicktag datasets
    ctdat_list = unique(sca_m$dataset [ !is.na(sca_m$dataset) ])
    # sca_m$dataset [ is.na(sca_m$dataset) ] = "single"
    # sca_m$dataset = factor(sca_m$dataset, levels = c("single",ctdat_list))
    ctdat_per_mc = table(sca_m$metacell, sca_m$dataset)
    ctdat_per_mc_max = apply(ctdat_per_mc, 1, max)
    ctdat_per_mc_sum = apply(ctdat_per_mc, 1, sum)
    ctdat_per_mc_max_dualmaxfrac = ctdat_per_mc_max / ctdat_per_mc_sum

    # plot
    pdf(sprintf("results_metacell_it1/quality_control_doublets.%s.pdf", spi), height = 5, width = 40)
    # fraction of dual sample cells (from ct runs)
    fraction_dual_sample_cells = table(sca_m$metacell, is.na(sca_m$ct_relative_size_ft), useNA = "ifany")
    # fraction_dual_sample_cells = fraction_dual_sample_cells / rowSums(fraction_dual_sample_cells)
    barplot(t(fraction_dual_sample_cells),
            col = c("blue","gray"), xlab = "metacell",
            main = sprintf("%s | fraction of dual-sample cells (from CT runs)", spi), las = 2, cex.names = 0.7)
    legend("topright", legend = c("dual sample","single sample"), fill = c("blue","gray"), bty = "n", cex = 0.8)
    # amongst ct cells, relative size of each sps ct
    boxplot(
        sca_m$ct_relative_size_ft ~ sca_m$metacell,
        col = mct$color, log = "y",
        xlab = "metacell",
        main = sprintf("%s | ratio 1st to 3rd CT", spi),
        ylab = "ratio", las = 2, ylim = c(1,1000), cex.axis = 0.7)
    abline(h=c(5,8,10), lty = 2, col = "darkred")
    # amongst ct cells, relative size of each sps ct
    boxplot(
        sca_m$ct_relative_size_fs ~ sca_m$metacell,
        col = mct$color, 
        xlab = "metacell",
        main = sprintf("%s | ratio 1st to 2nd CT", spi),
        ylab = "ratio", las = 2, outline = FALSE, cex.axis = 0.7)
    # amongst ct cells, distribution of total ct counts
    boxplot(
        sca_m$ct_total_counts ~ sca_m$metacell,
        col = mct$color, log = "y",
        xlab = "metacell",
        main = sprintf("%s | distribution of total ct counts", spi), las = 2, cex.axis = 0.7)
        
    # mean ratio first-to-third
    thrs_ratio = 10
    mean_ratio = data.frame(aggregate(ct_relative_size_ft ~ metacell, data = sca_m, function(i) quantile(i, 0.5)))
    fraction_dual = fraction_dual_sample_cells[,1] / rowSums(fraction_dual_sample_cells)
    # how many cells in that metacell come from clicktagged samples?
    mean_ratio$fraction_dual = 0
    mean_ratio$fraction_dual = fraction_dual [ mean_ratio[,1] ]
    # within cells that come from clicktagged samples, how many come from the same sample? (larest)
    mean_ratio$fraction_same_within_dual = 0
    mean_ratio$fraction_same_within_dual = ctdat_per_mc_max_dualmaxfrac [ mean_ratio[,1] ]
    
    
    # find mcs to drop:
    # * mean ratio first-to-third ratio > 10
    # * >50% cells from a clicktagged sample
    # * >90% cells from clicktagged samples come from the same sample
    ixs_to_drop = which(mean_ratio$fraction_same_within_dual > 0.9 & mean_ratio$fraction_dual > 0.5 & mean_ratio$ct_relative_size_ft > thrs_ratio)
    mcs_to_drop = mean_ratio[,1] [ ixs_to_drop ]
    scs_to_drop = unique(sca [ sca$metacell %in% mcs_to_drop , "cell"])
    message(sprintf("drop mcs %s: %s (%i mcs, %i cells)", spi, paste(mcs_to_drop, collapse = ","), length(mcs_to_drop), length(scs_to_drop)))
    
    # print(head( mean_ratio [ order(mean_ratio[,2], decreasing = TRUE),] , n = 15 ))
    
    dev.off()
    gc()
    
  
}

message("All done!")



