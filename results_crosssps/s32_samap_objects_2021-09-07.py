# libraries
import sys
import numpy as np
import pandas as pd
import samap
import samalg
from samap.mapping import SAMAP
from samap.utils import save_samap, load_samap

# output
out_fn = "results_samap_it4/"
bla_fn = "data_samap/blastdb/"

# list of pairwise species comparisons
sps_list = ["Tadh","TrH2","Hhon","HoiH23"]
com_list = [
	[ "Tadh", "TrH2"   ],
	[ "Tadh", "Hhon"   ],
	[ "Tadh", "HoiH23" ],
	[ "TrH2", "Hhon"   ],
	[ "TrH2", "HoiH23" ],
	[ "Hhon", "HoiH23" ]
]


### Pairwise analyses ###

for n,com in enumerate(com_list):
	
	# query species data
	sid1 = com[0]
	sid2 = com[1]

	### Load data ###

	# load SAM objects for species 1 and 2
	print("# loading %s data" % sid1)
	sam1 = samalg.SAM()
	sam1.load_data("%s/data.sam_object.%s.h5ad" % (out_fn, sid1))
	print("# loading %s data" % sid2)
	sam2 = samalg.SAM()
	sam2.load_data("%s/data.sam_object.%s.h5ad" % (out_fn, sid2))
	# dictionary of sam files
	dict_sam = { 
		sid1 : sam1, 
		sid2 : sam2 
	}
	
	# run samap, all homlogs with leiden clusters
	print("# load SAMAP, all homologs with leiden as clusters %s-%s" % (sid1,sid2))
	samm_al = SAMAP(dict_sam, f_maps = bla_fn, keys = { sid1 : "leiden_clusters", sid2 : "leiden_clusters" })
	samm_al.run()
	save_samap(samm_al, "%s/data.samap_object.all.lei.%s-%s" % (out_fn, sid1, sid2))

	# run samap, all homlogs with metacells
	print("# load SAMAP, all homologs with metacells as clusters %s-%s" % (sid1,sid2))
	samm_am = SAMAP(dict_sam, f_maps = bla_fn, keys = { sid1 : "metacell", sid2 : "metacell" })
	samm_am.run()
	save_samap(samm_am, "%s/data.samap_object.all.mcs.%s-%s" % (out_fn, sid1, sid2))

	# run samap, all homlogs with cell types
	print("# load SAMAP, all homologs with cell types as clusters %s-%s" % (sid1,sid2))
	samm_ac = SAMAP(dict_sam, f_maps = bla_fn, keys = { sid1 : "cell_type", sid2 : "cell_type" })
	samm_ac.run()
	save_samap(samm_ac, "%s/data.samap_object.all.cts.%s-%s" % (out_fn, sid1, sid2))

print("all pairs done!")



### Fourway analysis ###

# sam dict
dict_sam = dict()
for sps in sps_list:
	dict_sam[sps] = samalg.SAM()
	dict_sam[sps].load_data("%s/data.sam_object.%s.h5ad" % (out_fn, sps))

# run samap, all homlogs with leiden clusters
print("# load SAMAP, all homologs with leiden as clusters 4-way")
dict_clu = dict()
for sps in sps_list:
	dict_clu[sps] = "leiden_clusters"
samm_al = SAMAP(dict_sam, f_maps = bla_fn, keys = dict_clu)
samm_al.run()
save_samap(samm_al, "%s/data.samap_object.all.lei.4sps" % (out_fn))

# run samap, all homlogs with metacells
print("# load SAMAP, all homologs with metacells as clusters 4-way")
dict_clu = dict()
for sps in sps_list:
	dict_clu[sps] = "metacell"
samm_am = SAMAP(dict_sam, f_maps = bla_fn, keys = dict_clu)
samm_am.run()
save_samap(samm_am, "%s/data.samap_object.all.mcs.4sps" % (out_fn))

# run samap, all homlogs with cell types
dict_clu = dict()
for sps in sps_list:
	dict_clu[sps] = "cell_type"
print("# load SAMAP, all homologs with cell types as clusters 4-way")
samm_ac = SAMAP(dict_sam, f_maps = bla_fn, keys = dict_clu)
samm_ac.run()
save_samap(samm_ac, "%s/data.samap_object.all.cts.4sps" % (out_fn))

print("4sps done!")

