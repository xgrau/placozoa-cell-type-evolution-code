#### Input ####

# libraries
library("optparse")

option_list = list(
	# input and output
	optparse::make_option(
		c("-g", "--gtf"),
		help="Path to input GTF file (mandatory).",
		type="character", 
		default=NULL, 
		metavar="character"),	
	optparse::make_option(
		c("-p", "--peaks"),
		help="Path to MACS2 broad peaks file (including both + and - peaks).",
		type="character", 
		default=NULL, 
		metavar="character"),
	optparse::make_option(
		c("-o", "--out"),
		help="Path to output GTF file (contains only genes). Default is `new_genes.gtf`",
		type="character", 
		default="new_genes.gtf",
		metavar="character"),
	# logical options
	optparse::make_option(
		c("-t", "--trim_5p"),
		help="If set, trim 5' ends of genes that are overlapping with 3' ends of other genes. Default is not to (unset).",
		type="logical",
		default = FALSE,
		action = "store_true",
		metavar="character"),
	optparse::make_option(
		c("-a", "--add_orphan"),
		help="If set, add orphan peaks to the output GTF file. Default is not to (unset).",
		type="logical",
		default = FALSE,
		action = "store_true",
		metavar="character"),
	# other options
	optparse::make_option(
		c("-m", "--max_dist"),
		help="Max distance between a gene and a downstream peak that'll be assigned to this gene. Default is 5000 bp. Peaks beyond this distance are ignored or designated as orphan.",
		type="numeric",
		default=5e3,
		metavar="numeric"),
	optparse::make_option(
		c("-q", "--qval_thr"),
		help="q-value threshold to filter peaks in broadPeaks output (from MACS2).",
		type="numeric",
		default=1e-3,
		metavar="numeric"),
	optparse::make_option(
		c("-r", "--prefix_orphan"),
		help="Prefix for orphan peaks that will be added as genes in the output GTF. Default is NULL (none).",
		type="character",
		default=NULL,
		metavar="character"),
	optparse::make_option(
		c("-x", "--orphan_merge_distance"),
		help="Merge orphan peaks within this distance to each other. Default is NULL (median of gene length).",
		type="numeric",
		default=NULL,
		metavar="character")
)

# parse options
opp = OptionParser(option_list=option_list)
opt = parse_args(opp)

# check
if (is.null(opt$gtf) | is.null(opt$peaks)){
	print_help(opp)
	stop("Provide a GTF file with genes (-g/--gtf) and a BED file with peaks (-p/--peaks)", call.=FALSE)
}

# input files
gtf_fn = opt$gtf
p3p_fn = opt$peaks
out_fn = opt$out
# input logicals
trim_5p = opt$trim_5p
add_orphan = opt$add_orphan
# other options
orphan_prefix = opt$prefix_orphan
max_dist = opt$max_dist
qval_thr = opt$qval_thr
orphan_merge_distance = opt$orphan_merge_distance

# log start
message(sprintf("extend 3' | Input annotation: %s", gtf_fn))

# libraries
suppressPackageStartupMessages(library("GenomicRanges", quietly = TRUE))
suppressPackageStartupMessages(library("rtracklayer", quietly = TRUE))

#### Functions ####

# function to trim 5' ends of overlapping genes from a GRanges object
# containing gene ranges. Specifically, it removes 5' regions of genes 
# that overlap with 3' regions of other genes.
trim_5p_ends = function(gti) {
	
	# find overlapping genes and remove 5' regions of the overlaps 
	ovs_self = GenomicRanges::findOverlaps(gti, ignore.strand = TRUE, type = "any", drop.self = TRUE, drop.redundant = FALSE)
	
	if (length(ovs_self) > 0)	{ 
		
		# set aside genes that don't need to be modified (they don't overlap with any other gene)
		gti_f = gti [ !1:length(gti) %in% ovs_self@from ]
		
		# find genes that do overlap with other genes
		gti_self_ovs = gti[ c(ovs_self@from) ]
		gti_self_ovs_disj = GenomicRanges::disjoin(gti_self_ovs, with.revmap=TRUE, ignore.strand = TRUE)
		gti_self_ovs_disj = gti_self_ovs_disj [ lengths(gti_self_ovs_disj$revmap) == 2 ]
		
		# find ends of genes that are internally overlapping
		gti_self_ovs_ends = find_gene_ends(gti_self_ovs)
		
		# select overlapping regions that are nearest to the end of one gene
		keep_ixs = GenomicRanges::nearest(gti_self_ovs_disj, gti_self_ovs_ends)
		
		# keep closest genes as is
		gti_self_keep = gti_self_ovs [ keep_ixs ]
		
		# for non-closest genes, subtract the overlapping region
		gti_self_mods = sort(gti_self_ovs [ ! 1:length(gti_self_ovs) %in% keep_ixs ])
		gti_self_ovs_disj = sort(gti_self_ovs_disj)
		gti_self_subt = GenomicRanges::setdiff(gti_self_mods, gti_self_ovs_disj, ignore.strand = TRUE)
		start(gti_self_mods) = start(gti_self_subt)
		end(gti_self_mods) = end(gti_self_subt)
		
		# add modified gene ranges to the set-aside genes
		message(sprintf("extend 3' | Trim overlapping 5' ends | n = %i / %i modified genes", length(gti_self_mods), length(gti)))
		gto = c(gti_f, gti_self_keep, gti_self_mods)
		
	} else {
		
		message(sprintf("extend 3' | Trim overlapping 5' ends | n = %i / %i modified genes", 0, length(gti)))
		gto = gti
		
	}
	
	# return
	return(gto)
	
}


# create a GRanges object of gene 3' ends (with width = 1)
find_gene_ends = function(gti) {
	
	# init: get ends of genes (width 1)
	gti_ends = GenomicRanges::GRanges(
		seqnames = gti@seqnames,
		gene_id = gti@elementMetadata$gene_id,
		ranges = IRanges::IRanges(start = end(gti), end = end(gti)),
		strand = strand(gti)
	)
	
	# for genes in the minus strand, get start coordinate instead
	ixs_minus = as.character(strand(gti)) == "-"
	start(gti_ends) [ ixs_minus ] = start(gti) [ ixs_minus ]
	width(gti_ends) = 1
	
	# return
	return(gti_ends)
	
}

# extend the 3' regions of genes to include the span of downstream RNA-seq peaks, 
# including peaks that overlap with said gene and orphan peaks that are located
# downstream of the original gene body but within a certain `max_dist`.
extend_genes_to_p3p = function(gti, p3p, max_dist = 5000) {
	
	# find ends of genes
	gti_ends = find_gene_ends(gti)
	
	# find peaks that overlap with genes or ends of genes
	ovs_p3p_gti_e = GenomicRanges::findOverlaps(p3p, gti_ends, ignore.strand = FALSE, type = "any")
	ovs_p3p_gti_a = GenomicRanges::findOverlaps(p3p, gti, ignore.strand = TRUE, type = "any")
	
	## Extend genes with orphan peaks ##
	# find orphan peaks
	p3p_orphan = p3p [ ! 1:length(p3p) %in% c(ovs_p3p_gti_e@from, ovs_p3p_gti_a@from) ]
	# find distance between peaks and nearest gene end
	hit_d = GenomicRanges::distanceToNearest(p3p_orphan, gti_ends, ignore.strand = FALSE, select = "all")
	hit_t = data.frame(
		gene = gti_ends [ hit_d@to ]$gene_id,
		gene_end_pos = end(gti_ends [ hit_d@to ]),
		gene_strand = strand(gti_ends [ hit_d@to ]),
		peak = p3p_orphan [ hit_d@from ]$name,
		peak_pos = mid(p3p_orphan [ hit_d@from ]),
		peak_strand = strand(p3p_orphan [ hit_d@from ]),
		dist = hit_d@elementMetadata$distance
	)
	
	# check if this peak is downstream from gene end
	hit_t$peak_is_downstream = NA
	is_plus = hit_t$peak_strand == "+"
	hit_t[ is_plus,"peak_is_downstream"] = hit_t[ is_plus,"gene_end_pos"] < hit_t[ is_plus,"peak_pos"]
	hit_t[!is_plus,"peak_is_downstream"] = hit_t[!is_plus,"gene_end_pos"] > hit_t[!is_plus,"peak_pos"]
	
	# keep only peaks taht are downstream to the gene of interest
	hit_t = hit_t [ hit_t$peak_is_downstream, ]
	
	# apply max distance threshold
	hit_t = hit_t [ hit_t$dist <= max_dist,]
	
	# keep most distant peak
	hit_t = hit_t [ order(hit_t$gene, -hit_t$dist),  ]
	hit_t = hit_t [ !duplicated(hit_t$gene) ,]
	
	# extend genes that are near orphan downstream peaks
	gti_with_orph = gti [ hit_t$gene ]
	p3p_with_orph = p3p [ hit_t$peak ]
	gto_with_orph = GenomicRanges::punion(gti_with_orph, p3p_with_orph, fill.gap = TRUE)
	
	# keep gene and transcript ides
	gto_with_orph$gene_id = gti_with_orph$gene_id
	gto_with_orph$transcript_id = gti_with_orph$transcript_id
	
	
	## Extend genes with overlapping peaks ##
	# find new ends of genes
	gti_ends = find_gene_ends(gti)
	
	# find peaks overlapping with gene ends
	gti_with_end_ovs = gti [ ovs_p3p_gti_e@to ]
	p3p_with_end_ovs = p3p [ ovs_p3p_gti_e@from ]
	
	# extend genes that overlap with said peaks
	gto_with_ovs = GenomicRanges::punion(gti_with_end_ovs, p3p_with_end_ovs)
	
	# keep gene and transcript ides
	gto_with_ovs$gene_id = gti_with_end_ovs$gene_id
	gto_with_ovs$transcript_id = gti_with_end_ovs$transcript_id
	
	# remove genes that had been expanded by orphan peaks
	gto_with_ovs = gto_with_ovs [ !names(gto_with_ovs) %in% names(gto_with_orph) ]
	
	## Create new ranges ##	
	# names(gto_with_ovs) = NULL
	# names(gto_with_orph) = NULL
	# # concatenate	
	# gto = c(gto_with_ovs, gto_with_orph)
	# # select longest gene
	# # reduce
	# gto = GenomicRanges::reduce(gto, )
	
	# concatenate modified genes
	gto               = c(gto_with_orph,               gto_with_ovs)
	gto$gene_id       = c(gto_with_orph$gene_id,       gto_with_ovs$gene_id)
	gto$transcript_id = c(gto_with_orph$transcript_id, gto_with_ovs$transcript_id)
	gto$source        = "extend3p"
	
	# add them to unmmodified genes
	gtt = gti
	gtt [ names(gtt) %in% names(gto) ] = NULL # remove modified entries
	gtt$source = "original"
	gtt = c(gtt, gto) # reconcatenate
	gtt = sort(gtt)   # sort
	gtt$type = "gene"
	names(gtt) = NULL
	
	# return
	message(sprintf("extend 3' | Extend 3' ends | n = %i / %i modified genes", length(gto),length(gtt)))
	return(list(new_genes = gtt, genes_expanded_by_orphan = gto_with_orph, genes_expanded_by_overlaps = gto_with_ovs))
	
}

# define unassigned peaks as additional orphan genes
unassigned_p3p_as_orphan = function(gti, p3p, orphan_prefix = NULL, orphan_merge_distance = NULL) {
	
	# get prefix for orphan peaks
	if (is.null(orphan_prefix)) {
		orphan_string = "orphan"
	} else {
		orphan_string = paste(orphan_prefix, "orphan", sep = "_")
	}
	
	# merge orphan peaks within this distance (median gene length by default)
	if (is.null(orphan_merge_distance)) {
		orphan_merge_distance = median(width(gti))
	}
	
	# find orphan peaks
	ovs_any = GenomicRanges::findOverlaps(p3p, gti, ignore.strand = TRUE, type = "any")
	p3p_orphan = p3p [ ! 1:length(p3p) %in% ovs_any@from ]
	
	# remerge orphan peaks within this distance
	message(sprintf("extend 3' | Add orphan peaks | found %i peaks, merge at %s bp", length(p3p_orphan), orphan_merge_distance))
	p3p_orphan_e = p3p_orphan + orphan_merge_distance / 2
	p3p_orphan_e = GenomicRanges::reduce(p3p_orphan_e)
	p3p_orphan   = p3p_orphan_e - orphan_merge_distance / 2
	
	# reformat as gtf
	p3p_orphan$gene_id = paste(orphan_string, 1:length(p3p_orphan), sep = "_")
	p3p_orphan$transcript_id = paste(orphan_string, 1:length(p3p_orphan), sep = "_")
	names(p3p_orphan) = NULL
	p3p_orphan$type = "gene"
	p3p_orphan$source = "orphan"
	
	# concatenate
	message(sprintf("extend 3' | Add orphan peaks | n = %i new features", length(p3p_orphan)))
	gto = sort(c(gti, p3p_orphan))
	
	# output 
	return(gto)
	
}

#### Main ####

# read input GTF with genes
gtf = rtracklayer::import(gtf_fn, format = "gtf")
# get original genes (or transcripts)
if ("gene" %in% gtf$type) {
	gti = gtf [ gtf$type == "gene", ]
} else if ("mRNA" %in% gtf$type) {
	gti = gtf [ gtf$type == "mRNA", ]
} else if ("transcript" %in% gtf$type) {
	gti = gtf [ gtf$type == "transcript", ]
}
message(sprintf("extend 3' | Load n = %i genes", length(gti)))

# read peaks
p3t = read.table(p3p_fn, sep = "\t", header = FALSE, col.names = c("chr","start","end","name","score","strand","foldchange","minlog10pval","minlog10qval"))
message(sprintf("extend 3' | Load n = %i peaks", nrow(p3t)))
p3t$qval = 10 ^ (- p3t$minlog10qval)
if (!is.null(qval_thr)) {
	p3t = p3t [ p3t$qval < qval_thr, ]
	message(sprintf("extend 3' | Filter peaks at q<%s, n = %i remain", as.character(qval_thr), nrow(p3t)))
}
p3p = GenomicRanges::GRanges(
	seqnames = p3t$chr,
	name = p3t$name,
	ranges = IRanges::IRanges(start = p3t$start + 1, end = p3t$end + 1),
	strand = p3t$strand
)

# add names
names(gti) = gti@elementMetadata$gene_id
names(p3p) = p3p$name

# extend genes to 3p peaks
gti_extended = suppressWarnings(extend_genes_to_p3p(gti, p3p, max_dist = max_dist))
gti = gti_extended$new_genes

# add orphan peaks as extra genes if available
if (add_orphan) {
	gti = unassigned_p3p_as_orphan(gti, p3p, orphan_prefix = orphan_prefix, orphan_merge_distance = orphan_merge_distance)
}

# trim 5p regions of genes that are overlapping with 3p regions of other genes
if (trim_5p) {
	gti = trim_5p_ends(gti)
}

# output
rtracklayer::export(gti, con = out_fn, format = "gtf")

# log data sources
pdf(gsub("\\.gtf$",".pdf",out_fn), height = 4, width = 4)
pie(table(gti$source), labels = sprintf("%s\nn=%i", names(table(gti$source)), table(gti$source)), cex = 0.8, col = c("springgreen3","gray95","violet"))
title(main = sprintf("n=%i features", length(gti)))
sink = dev.off()

message(sprintf("extend 3' | Done! %i reannotated features in %s", length(gti), out_fn))

# for logging, export which genes have been extended as overlaps
# rtracklayer::export(gti_extended$genes_expanded_by_overlaps, con = sprintf("%s/%s_3pextended.ovs.gtf", out_fo, spi), format = "gtf")
# rtracklayer::export(gti_extended$genes_expanded_by_orphan, con = sprintf("%s/%s_3pextended.orphan.gtf", out_fo, spi), format = "gtf")


