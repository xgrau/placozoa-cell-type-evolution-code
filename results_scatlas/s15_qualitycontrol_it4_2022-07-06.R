# libraries
library("metacell")
library("scales")
source("../scripts/Downstream_functions.R")
source("../scripts/Modified_functions.R")
source("../scripts/helper.R")
par(family  = "Arial")

# where to store output?
out_fn = "results_metacell_it4/"

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

for (spi in sps_list) {
	
    # metacell and single cell annots
    mct_fn = sprintf("results_metacell_it4/annotation_mc.%s.it4.reordered.tsv", spi)
    sca_fn = sprintf("results_metacell_it4/scdr_%s.matrix.sc_annot.csv", spi)
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
    pdf(sprintf("results_metacell_it4/quality_control_doublets.%s.pdf", spi), height = 5, width = 40)
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
    message(sprintf("drop mcs %s: %s", spi, paste(mcs_to_drop, collapse = ",")))
    # print(head( mean_ratio [ order(mean_ratio[,2], decreasing = TRUE),] , n = 15 ))
    
    dev.off()
    gc()
    
    # Identify metacells that are intermediate between annotated cell types

    # cell types
    ctt = mct
    list_cell_types = unique(ctt$cell_type)
    list_cell_types = list_cell_types [ !grepl("trans", list_cell_types) & !grepl("peptidergic", list_cell_types) ]
    
    # init database
	metacell::scdb_init("data/scdb/",force_reinit=TRUE)
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it4",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it2",run_name))
    mc_counts = sca_mc_gene_counts(mc, mat)
    mc_umifrac = sca_mc_gene_umifrac(mc, mc_counts)
    mc_ct = sca_cell_type_fp(ctt[,c(1,4,3)], mc, mat)
    
    list_cell_types = unique(ctt$broad_cell_type)
    list_cell_types = list_cell_types [ list_cell_types != "trans" & list_cell_types != "unknown" ]
    
    pdf(sprintf("results_metacell_it4/quality_control_celltypes.%s.pdf", spi), height = 9, width = 9)
    layout(mat = matrix(1:9, nrow = 3))
    for (ni in 1:length(list_cell_types)) {
        for (nj in 1:length(list_cell_types)) {
            if (ni < nj){
                
                # mcs in each cell type
                cti = list_cell_types[ni]
                ctj = list_cell_types[nj]
                mc_cti = ctt [ ctt$cell_type == cti , "metacell" ]
                mc_ctj = ctt [ ctt$cell_type == ctj , "metacell" ]
                
                # find markers, cell type i
                mk_cti = data.frame(cti = mc_ct@mc_fp[,cti], cto = apply(mc@mc_fp[, ! colnames(mc@mc_fp) %in% mc_cti], 1, mean))
                mk_cti = mk_cti [ mk_cti$cti > 1.1 & mk_cti$cto < 1.1, ]
                mk_cti = mk_cti [ order(mk_cti$cti, decreasing = TRUE) , ]
                mk_cti = rownames(head(mk_cti, 200))
                # find markers, cell type j
                mk_ctj = data.frame(ctj = mc_ct@mc_fp[,ctj], cto = apply(mc@mc_fp[, ! colnames(mc@mc_fp) %in% mc_ctj], 1, mean))
                mk_ctj = mk_ctj [ mk_ctj$ctj > 1.1 & mk_ctj$cto < 1.1, ]
                mk_ctj = mk_ctj [ order(mk_ctj$ctj, decreasing = TRUE) , ]
                mk_ctj = rownames(head(mk_ctj, 200))
                # scores
                sc_cti = scale(colSums(mc_umifrac [ mk_cti, ]))
                sc_ctj = scale(colSums(mc_umifrac [ mk_ctj, ]))
                # plot
                plot(
                    sc_cti, sc_ctj, col = ctt$color, 
                    xlab = cti, ylab = ctj, pch = 19, 
                    # xlim = c(-2,8), ylim = c(-2,8),
                    main = sprintf("%s %s\n%i %s & %i %s markers", cti, spi, length(mk_cti), cti, length(mk_ctj), ctj), 
                    cex.main = 0.7, cex.axes = 0.7
                )
                text(sc_cti, sc_ctj, labels = ctt$metacell, col = alpha("black",0.9), cex = 0.5)
                abline(h=0, v=0, a=0, b= 1, lty = 2)
                
            }
        }

    }
    plot(0,0, col = NA)
    legend("topleft", legend = unique(ctt$cell_type), col = unique(ctt$color), cex = 0.6, bty = "n", pch = 19, ncol = 4)
    dev.off()
    
    
    library("umap")

    # mcs
    mcs_terminal = ctt [ ctt$broad_cell_type != "trans", "metacell" ]
    mcs_puttrans = ctt [ ctt$broad_cell_type == "trans", "metacell" ]
    # cols
    col_terminal = ctt [ ctt$broad_cell_type != "trans", "color" ]
    col_puttrans = ctt [ ctt$broad_cell_type == "trans", "color" ]

    # reference umap
    mc_umifrac_t = mc_umifrac [ , mcs_terminal ]
    um_mcumi_t = umap::umap(d=t(mc_umifrac_t))

    # predicted umap
    mc_umifrac_p = mc_umifrac [ , mcs_puttrans ]
    um_mcumi_p = predict(um_mcumi_t, t(mc_umifrac_p))

    # plot
    pdf(sprintf("results_metacell_it4/quality_control_umap.%s.pdf", spi), height = 9, width = 9)
    plot(um_mcumi_t$layout, col = col_terminal, pch = 19)
    points(um_mcumi_p, col = col_puttrans, pch = 19)
    text(um_mcumi_t$layout, colnames(mc_umifrac_t), cex = 0.5, col = alpha("black", 0.9))
    text(um_mcumi_p, paste0("t",colnames(mc_umifrac_p)), cex = 0.5, col = alpha("darkred", 0.9))
    dev.off()

}

message("All done!")



