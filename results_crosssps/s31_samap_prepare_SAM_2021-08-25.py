# libraries
import pandas as pd
import numpy as np
import sys
from scipy import io
# path to samap
sys.path.append("/users/asebe/xgraubove/Programes/samap_directory/samap")
import samalg
# path to modified samap functions
sys.path.append("../scripts")

# output
out_fn = "results_samap_it4/"
sps_list = ["Tadh","TrH2","Hhon","HoiH23"]

# create SAM object
def create_sam_from_objects(mat, met):
	if np.any(np.isin(mat.columns[0] , ["Genes","Gene", "gene", "genes", 'Unnamed: 0'])):
		mat.index = mat.iloc[:,0]
		mat = mat.drop(columns=mat.columns[0])
	mat = mat.transpose()
	sam = samalg.SAM(counts=mat)   # initialise SAM dataset
	sam.preprocess_data()          # log-transform and filter data
	# add annotations
	for i in range(met.shape[1]):
		sam.adata.obs[met.columns[i]] = met.iloc[:,i].values
	return sam


#### Prepare single-species SAM objects ####

# reference species
for sps in sps_list :

	# count sparse matrix
	mati_fn = "../results_scatlas/results_metacell_it4/scdr_%s.matrix.sc_umi.mtx.gz" % sps
	celi_fn = "../results_scatlas/results_metacell_it4/scdr_%s.matrix.sc_umi.annot_cells.txt" % sps
	geni_fn = "../results_scatlas/results_metacell_it4/scdr_%s.matrix.sc_umi.annot_genes.txt" % sps
	# cell metadata
	meti_fn = "../results_scatlas/results_metacell_it4/scdr_%s.matrix.sc_annot.csv" % sps

	# read
	print("# loading %s" % (sps))
	matt = io.mmread(mati_fn)
	matt_d = pd.DataFrame.sparse.from_spmatrix(matt)
	# cells and genes
	celt = pd.read_csv(celi_fn, sep="\t", header = None, names = ["cell"])
	gent = pd.read_csv(geni_fn, sep="\t", header = None, names = ["gene"])
	matt_d.columns = celt["cell"].values
	matt_d.index = gent["gene"]
	# metadata
	mett = pd.read_csv(meti_fn, sep="\t")

	# create SAM objects for concatenated dataset
	print("# run SAM %s" % (sps))
	sam = create_sam_from_objects(mat = matt_d, met = mett)
	sam_f = sam
	sam_f.run()

	# save precomputed sam object
	print("# save SAM %s" % (sps))
	sam_f.save_anndata("%s/data.sam_object.%s.h5ad" % (out_fn, sps))

	# convert matt to sparse dataframe
	matt_s = matt_d.astype(pd.SparseDtype("int",0))

	# save concatenated dataset
	print("# Saving tables %s...", sps)
	met_fn = "%s/data.sam_object.%s.sc_annotations.tsv" % (out_fn, sps)
	mat_fn = "%s/data.sam_object.%s.sc_counts.tsv.pkl" % (out_fn, sps)
	mett.to_csv(met_fn, sep = "\t", index=False)
	matt_s.to_pickle(mat_fn)

print("all done!")
