# libraries
library("metacell")
source("../scripts/helper.R")
par(family  = "Arial")

# initial metacell iteration stored here
out_fn = "results_metacell_it0/"

# species list
sps_list = c("Tadh","TrH2","Hhon","HoiH23")
mcs_to_drop_list = list(
	"Tadh" = c(68,69,125,129,130,131,132,133,134),
	"TrH2" = c(),
	"Hhon" = c(254),
	"HoiH23" = c()
)

for (spi in sps_list) {

	# first mc solution
	metacell::scdb_init("data/scdb/", force_reinit=TRUE)
	run_name = sprintf("scdr_%s", spi)
	mc = metacell::scdb_mc(sprintf("%s_it1",run_name))
	mat = metacell::scdb_mat(sprintf("%s_it1",run_name))

	# load cell type annotations
	ctt_fn = sprintf("data/annotation_mc.%s.it1.tsv", spi)
	ctt = read.table(ctt_fn, header = TRUE, comment.char = "", sep = "\t")

	# get metacells and cells to drop
	mcs_to_drop = as.character(mcs_to_drop_list[[spi]])
	
	cells_to_drop = names(mc@mc) [ mc@mc %in% as.numeric(mcs_to_drop) ]
	cells_to_drop_mc = mc@mc [ mc@mc %in% as.numeric(mcs_to_drop) ]
		
	# save new mat object for it1
	message(sprintf("save mat.%s_it2...", run_name))
	mcell_mat_ignore_cells(
		new_mat_id = sprintf("%s_it2",run_name), 
		mat_id = sprintf("%s_it1",run_name),
		ig_cells = union(cells_to_drop, mat@ignore_cells))
	
	# log
	message(sprintf("%s drop: %i cells", run_name, length(cells_to_drop)))
	mat0 = scdb_mat(sprintf("%s_it1", run_name))
	mat1 = scdb_mat(sprintf("%s_it2", run_name))
	message(sprintf("%s old size: %i cells", run_name, ncol(mat0@mat)))
	message(sprintf("%s new size: %i cells", run_name, ncol(mat1@mat)))

}

message("all done!")
