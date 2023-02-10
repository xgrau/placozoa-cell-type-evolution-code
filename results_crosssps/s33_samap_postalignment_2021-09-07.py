# libraries
import sys
import numpy as np
import pandas as pd
# path to samap
sys.path.append("/users/asebe/xgraubove/Programes/samap_directory/samap")
import samap
from samap.utils import save_samap, load_samap
# import modified functions
sys.path.append("../scripts")
import samap_mods as smods

# output
out_fn = "results_samap_it4/"

# list of pairwise species comparisons
com_list = [
	[ "Tadh", "TrH2"   ],
	[ "Tadh", "Hhon"   ],
	[ "Tadh", "HoiH23" ],
	[ "TrH2", "Hhon"   ],
	[ "TrH2", "HoiH23" ],
	[ "Hhon", "HoiH23" ]
]

dict_clu = {
	"lei" : "leiden_clusters",
	"mcs" : "metacell",
	"cts" : "cell_type"
}

### Pairwise ###

for n,com in enumerate(com_list):
	
	# query species data
	sid1 = com[0]
	sid2 = com[1]
	
	### Load data ###
	
	for cli in [ "mcs","cts" ] :

		# get cluster type label (leiden_clusters for `lei`)		
		clin = dict_clu[cli]
		
		# load
		print("%s-%s | SAMAP load object, all homologs with %s as clusters" % (sid1, sid2, clin))
		samm = load_samap("%s/data.samap_object.all.%s.%s-%s" % (out_fn, cli, sid1, sid2))

		### Mapping scores ###
		
		# print("%s-%s | SAMAP alignment scores, %s" % (sid1,sid2, clin))
		map_score = samap.analysis.get_mapping_scores(sm = samm, keys = { sid1: clin, sid2: clin }, n_top = 0)
		map_score[0].to_csv("%s/csps_samap.%s-%s.%s.ntop0.score_ranked.csv" % (out_fn, sid1, sid2, cli), sep = "\t")
		map_score[1].to_csv("%s/csps_samap.%s-%s.%s.ntop0.score_matrix.csv" % (out_fn, sid1, sid2, cli), sep = "\t")
		map_score = samap.analysis.get_mapping_scores(sm = samm, keys = { sid1: clin, sid2: clin }, n_top = 100)
		map_score[0].to_csv("%s/csps_samap.%s-%s.%s.ntop100.score_ranked.csv" % (out_fn, sid1, sid2, cli), sep = "\t")
		map_score[1].to_csv("%s/csps_samap.%s-%s.%s.ntop100.score_matrix.csv" % (out_fn, sid1, sid2, cli), sep = "\t")
		
		### Shared genes ###
		# find shared genes across species
		print("%s-%s | SAMAP find shared gene pairs, %s" % (sid1,sid2,cli))
		gpf = samap.analysis.GenePairFinder(sm = samm, keys = { sid1: clin, sid2: clin })
		gpf_all_pairs = gpf.find_all(thr=0.05)
		gpf_all_pairs.to_csv("%s/csps_samap.%s-%s.%s.markers_shared.csv" % (out_fn, sid1, sid2, cli), index = False, sep = "\t")
		
		# ### Other ###
		# # store gene pairs
		# print("%s-%s | SAMAP store gene pair weights" % (sid1,sid2))
		# gps = samm.samap.adata.uns["gene_pairs"]
		# gps_sps1 = [ n.split(";")[0] for n in gps ]
		# gps_sps2 = [ n.split(";")[1] for n in gps ]
		# gpd = pd.DataFrame({ "species_1" : gps_sps1, "species_2" : gps_sps2, "weight" : samm.samap.adata.uns["edge_weights"] })
		# gpd.to_csv("%s/csps_samap.%s-%s.markers.pairs.csv" % (out_fn, sid1, sid2), index = False, sep = "\t")
		
		# coordinates of merged projection
		print("%s-%s | SAMAP store joint UMAP coordinates" % (sid1,sid2))
		sad = pd.concat([ samm.samap.adata.obs, samm.samap.adata.obsm.to_df()[["X_umap1","X_umap2"]] ], axis = 1)
		sad.to_csv("%s/csps_samap.%s-%s.%s.transfer_umap.csv" % (out_fn, sid1, sid2, cli), index = False, sep = "\t")




### 4sps ###

sps_list = [ "Tadh","TrH2","Hhon","HoiH23"]

for cli in [ "lei","mcs","cts" ] :

	# get cluster type label (leiden_clusters for `lei`)		
	clin = dict_clu[cli]
	
	# get dictionary of keys
	key_dict = dict()
	for spi in sps_list:
		key_dict[spi] = clin
	
	# load
	print("%s | SAMAP load object, all homologs with %s as clusters" % ("4sps", clin))
	samm = load_samap("%s/data.samap_object.all.%s.%s" % (out_fn, cli, "4sps"))

	print("%s | SAMAP alignment scores, %s" % ("4sps", clin))
	map_score = samap.analysis.get_mapping_scores(sm = samm, keys = key_dict, n_top = 0)
	map_score[0].to_csv("%s/csps_samap.%s.%s.ntop0.score_ranked.csv" % (out_fn, "4sps", cli), sep = "\t")
	map_score[1].to_csv("%s/csps_samap.%s.%s.ntop0.score_matrix.csv" % (out_fn, "4sps", cli), sep = "\t")
	map_score = samap.analysis.get_mapping_scores(sm = samm, keys = key_dict, n_top = 100)
	map_score[0].to_csv("%s/csps_samap.%s.%s.ntop100.score_ranked.csv" % (out_fn, "4sps", cli), sep = "\t")
	map_score[1].to_csv("%s/csps_samap.%s.%s.ntop100.score_matrix.csv" % (out_fn, "4sps", cli), sep = "\t")
	
	print("%s | SAMAP store joint UMAP coordinates" % ("4sps"))
	sad = pd.concat([ samm.samap.adata.obs, samm.samap.adata.obsm.to_df()[["X_umap1","X_umap2"]] ], axis = 1)
	sad.to_csv("%s/csps_samap.%s.%s.transfer_umap.csv" % (out_fn, "4sps", cli), index = False, sep = "\t")



print("all done!")
