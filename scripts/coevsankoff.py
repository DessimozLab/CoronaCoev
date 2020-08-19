import numpy as np
from collections import Counter
import multiprocessing as mp
import pandas as pd
from Bio import SeqIO
import pandas as pd
from colour import Color
import pickle
import time
import h5py
import dendropy
from Bio import AlignIO , SeqIO
from scipy import sparse
import pickle
import multiprocessing as mp
import copy
import sys
import os

sys.setrecursionlimit( 10 **5 )

runName = 'COEV_cleantree_mk4'

NCORE = 50
treefile = '/home/cactuskid13/covid/lucy_mk2/30_07_2020/gisaid_hcov-2020_07_30.QC.NSoutlier.filter.deMaiomask.aln.EPIID.HF.treefile'
alnfile = '/home/cactuskid13/covid/lucy_mk2/30_07_2020/gisaid_hcov-2020_07_30.QC.NSoutlier.filter.deMaiomask.EPIID.HF.noambig.HomoplasyFinder_input.aln'


print('mk5')

tree = dendropy.Tree.get(
    path=treefile,
    schema='newick')


msa = AlignIO.read(alnfile , format = 'fasta')
if os.path.exists(alnfile +'.h5'):
    with h5py.File(alnfile +'.h5', 'r') as hf:
        align_array = hf['MSA2array'][:]
        #implement w np unique could be faster
else:
    print('aln 2 numpy ')
    align_array = np.array([ list(rec.upper())  for rec in msa], np.character)
    print('done')
    print('dumping to hdf5')
    with h5py.File(alnfile +'.h5', 'w') as hf:
        hf.create_dataset("MSA2array",  data=align_array)
    print('done')

sites = {}

for col in range(align_array.shape[1]):
    if col% 1000  == 0:
        print(col)
    (unique, counts) = np.unique(align_array[:,col].ravel() , return_counts=True)
    sites.update({col:dict(zip(list(unique), list(counts)))})
informativesites = [ s for s in sites if len(set( sites[s].keys()) -set([b'-',b'N']) ) > 2  ]
print(len(informativesites))

for l in tree.leaf_nodes()[0:10]:
    print(str(l.taxon))

for i in range(10):
    print(align_array[i,:])

def clipID(ID):
    return ID.replace('|',' ').replace('_',' ').replace('/',' ').strip()

IDs = {i:clipID(rec.id) for i,rec in enumerate(msa)}
#IDs = {i:rec.id for i,rec in enumerate(msa)}
IDindex = dict(zip( IDs.values() , IDs.keys() ) )

print( [(t,IDindex[t]) for t in list(IDindex.keys())[0:10]] )

print('rev3')
import pdb; pdb.set_trace()


#def clipID(ID):
#    return ''.join( [ s +'|' for s in str(ID).split('|')[:-1] ])[:-1].replace('_',' ')

#def clipID(ID):
#    return ID.replace('|',' ').replace('_',' ').replace('/',' ').strip()


 #'msa :::: hCoV-19/Iceland/447/2020|EPI_ISL_424471|2020-03-20'

# 'tree :::: hCoV-19 Scotland EDB399 2020 EPI ISL 425998 2020-03-30'

#prep tree for sankoff
for i,n in enumerate(tree.nodes()):
    n.matrow = i
    n.symbols = None
    n.scores = None
    n.event = None
    n.char = None
matsize = len(tree.nodes())
print(matsize)
print('nodes')
#sites = { col: Counter( msa[:,col] ) for col in range(len(msa[1,:])) }
#place an event on a column with multiple symbols
allowed_symbols = { b'A', b'C', b'G' , b'T' }

def process_node_smallpars_1(node):
    #go from leaves up and generate character sets
    if node.symbols is None:
        for child in node.child_nodes():
            if child.symbols is None:
                process_node_smallpars_1(child)
        symbols = set.intersection( * [ child.symbols for child in node.child_nodes( ) ] )
        if len(symbols) == 0:
            symbols = set.union( * [ child.symbols for child in node.child_nodes( ) ] )

        node.symbols = symbols
        node.scores = { }
        for c in allowed_symbols:
            if c not in node.symbols:
                #add trnasition mat here if needed
                score = min(  [ child.scores[c] for child in node.child_nodes()])+1
            else:
                score = min(  [ child.scores[c] for child in node.child_nodes() ] )
            node.scores[c] = score

def process_node_smallpars_2(node):
    #assign the most parsimonious char from children
    if node.char is None:
        node.char = min(node.scores, key=node.scores.get)
        if node.parent_node:
            if node.parent_node.char == node.char:
                node.event = 0
            else:
                if node.scores[node.parent_node.char] == node.scores[node.char] :
                    node.char = node.parent_node.char
                    node.event = 0
                else:
                    node.event = 1
        else:
            node.event = 0
        for child in node.child_nodes():
            if child.char is None:
                process_node_smallpars_2(child)
def calculate_small_parsimony( t, aln_column , row_index , verbose  = True ):
    missing = 0
    #assign leaf values
    for l in t.leaf_nodes():
        l.event = 0
        l.scores = { c:10**10 for c in allowed_symbols }
        if str(l.taxon).replace("'", '') in row_index:
            char = aln_column[row_index[str(l.taxon).replace("'", '')]]
            if char.upper() in allowed_symbols:
                l.symbols = { char }
                l.scores[char] = 0
            else:
                #ambiguous leaf
                l.symbols =  allowed_symbols
        else:
            missing += 1
            char = None
            l.symbols =  allowed_symbols
            if verbose == True:
                print( 'err ! alncol: ', l.taxon , aln_column  )
        l.char = min(l.scores, key=l.scores.get)
    if verbose == True:
        print('done init')

        
        print(missing)
        print('missing in aln')
        for n in t.leaf_nodes()[0:5]:
            print(n)
            print(n.symbols)
            print(n.scores)
            print(n.char)
            print(n.event)


    #up
    process_node_smallpars_1(t.seed_node)
    #down
    process_node_smallpars_2(t.seed_node)

    eventindex = [ n.matrow for n in t.nodes()   if n.event > 0 ]
    if verbose == True:
        print(eventindex)

    events = np.zeros( (len(t.nodes()),1) )
    events[eventindex] = 1
    return events

def process(q,retq, iolock ,tree , IDindex   ):
    #calculate tree events
    with iolock:
        print('init worker')
    count = 0
    while True:
        stuff = q.get()
        if stuff is None:
            break
        index,aln_column = stuff

        retq.put((index,calculate_small_parsimony( copy.deepcopy(tree) , aln_column , IDindex ) ) )
        if count % 1000 == 0:
            with iolock:
                print(index)
                print('worker')
        count+= 1
    print('done')

def mat_creator(retq,matsize,iolock,runName = ''):
    with iolock:
        print('init matcreator')
    #collect distances and create final mat
    #distmat = np.zeros(  )
    M1 = sparse.lil_matrix((matsize[0],matsize[1] ))
    count = 0
    init = False
    t0 = time.time()
    while True:
        r = retq.get()
        if r is not None:
            init = True

        count+=1
        if r is None and init == True:
            break
        col,events = r
        M1[np.nonzero(events)[0] ,col] = 1
        if count == 0:
            with iolock:
                print(np.nonzero(events))
        if time.time()-t0> 120 :
            t0 = time.time()
            with iolock:
                print('saving')
                print(col)
                print(np.nonzero(events))

                with open( runName + 'coevmatrev3.pkl' , 'wb') as coevout:
                    coevout.write(pickle.dumps(M1))
                print('done')

def main():


    #reset tree
    for i,n in enumerate(tree.nodes()):
        n.matrow = i
        n.symbols = None
        n.scores = None
        n.event = None
        n.char = None
    #our results will be events on tree nodes for each column
    matsize = (len(tree.nodes()),align_array.shape[1])
    q = mp.Queue(maxsize=NCORE*1000)
    retq = mp.Queue(maxsize= NCORE*100 )
    iolock = mp.Lock()
    #start workers
    pool = mp.Pool(NCORE, initializer=process, initargs=(q,retq, iolock ,tree , IDindex ))
    #start saver
    p = mp.Process(target=mat_creator, args=(retq,matsize, iolock, alnfile + runName ))
    p.start()
    for i,k in enumerate(informativesites):
        s1 = align_array[:,k]
        q.put( (k,s1) )
    for _ in range(NCORE):  # tell workers we're done
        q.put(None)
    #retq.put(None)
    pool.close()
    pool.join()
    #retq.put(None)
    p.join()

if __name__ == '__main__':
    main()
