#collapse a tree into N seqments with equal branch lengths


import dendropy
import scipy.sparse as sparse
import pickle
import numpy as np
import sys

sys.setrecursionlimit( 10 **8 )



treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
#load tree
tree = dendropy.Tree.get( path=treefile, schema='newick')
treecut = 100000000
thresh = tree.length()/treecut
print(thresh)
#label sectors
def process_node(node , thresh ,total= 0, sector= None,  verbose = False):
		#assign the most parsimonious char from children
		if node.parent_node:
			node.sector = sector
		else:
			#root node
			total = 0
			node.sector = 0
		for child in node.child_nodes():
			#add each child's branch to the total
			total += child.edge.length
			if total > thresh:
				total = 0
				process_node(child, thresh, total, node.sector + 1 )
			else:
				process_node(child, thresh , total,  node.sector  )
process_node( tree.seed_node , thresh )
for i,n in enumerate(tree.nodes()):
	n.matrow = i
#sectors 2 mat, coordinates are matrow , sector 
#data = 1
row = [n.matrow for n in tree.nodes()]
col = [n.sector for n in tree.nodes()]
data = np.ones( (len(tree.nodes())))
sectormat = sparse.csc_matrix( (data,(row,col)) )
print(sectormat.shape)
print(sectormat)
with open( treefile+'sector.pkl'  , 'wb' ) as matout:
	matout.write(pickle.dumps(sectormat))
