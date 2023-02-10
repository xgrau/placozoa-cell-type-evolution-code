import sys
import p4 as p4
import argparse
import logging

# argument parser
arp = argparse.ArgumentParser()

# add the arguments to the parser
arp.add_argument("-a", dest="alig",  required=True,                help="path to input alignment (fasta)", type=str)
arp.add_argument("-t", dest="tree",  required=True,                help="path to input tree (newick)", type=str)
arp.add_argument("-o", dest="out",   required=False, default=None, help="path to output", type=str)
arp.add_argument("-m", dest="model", required=False, default="lg", help="empirical model to test (default is 'lg', can also accept 'wag', 'jtt', and more)", type=str)
arp.add_argument("-g", dest="gamma", required=False, default=4,    help="number of gamma categories to set, default is 4; set to 0 to ignore", type=int)
arp.add_argument("-n", dest="nsims", required=False, default=100,  help="number of simulations (default is 100)", type=int)
# arp.add_argument("-i", dest="invar", required=False, action="store_true", help="use this flag to set the model to use invariant sites.")
arp.add_argument("-b", dest="brlen", required=False, action="store_true", help="use this flag to reoptimise branch lengths given the above model; leave unset to use previously computed branch lengths")
arp = arp.parse_args()

# set logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)-5.5s]\t%(message)s", handlers=[ logging.StreamHandler() ] )
p4.var.doCheckForDuplicateSequences = False

# work
# load alignment
logging.info("Load %s alignment" % arp.alig)
p4.read(arp.alig)
a = p4.var.alignments[0]

# load tree
logging.info("Load %s tree" % arp.tree)
p4.read(arp.tree)
t = p4.var.trees[0]
d = p4.Data()
t.data = d

# define evolutionary model model (LG+G4)
logging.info("Use %s empirical model" % arp.model)
t.newComp(free=1, spec="empirical")
t.newRMatrix(free=0, spec=arp.model)
# should I use gamma categories?
if arp.gamma > 0:
	logging.info("Use %i gamma categories" % arp.gamma)
	t.setNGammaCat(nGammaCat=arp.gamma)
	t.newGdasrv(free=1, val=0.5)
# should I use invariant sites? always:
logging.info("Use invariant sites")
t.setPInvar(free=0, val=0.0)

# optimise tree branches and calculate likelihood
if arp.brlen:
	logging.info("Optimise branch lenghts under %s+G%i..." % (arp.model, arp.gamma))
	t.optLogLike()
else:
	logging.info("Use branch lengths from %s..." % arp.tree)
	

# simulation
logging.info("Compositional test with %i simulations..." % arp.nsims)
sim_chi_p = t.compoTestUsingSimulations(nSims = arp.nsims)

# store output
if arp.out is not None:
	handle = open(arp.out, mode = "w")
	handle.write(arp.out + "\t" + str(sim_chi_p) + "\n")
	handle.close()

logging.info("Test %s done!" % arp.alig)
