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


tree = dendropy.Tree.get(
    path='./lanford/ft_TBE.tree.txt',
    schema='newick')
print(dir(tree))

for l in tree.leaf_nodes()[0:10]:
    print(str(l.taxon))
print(len(tree.leaf_nodes()))
print('leaves')
for i,n in enumerate(tree.nodes()):
    n.matrow = i
    n.symbols = None
    n.scores = None
    n.event = None
    n.char = None


matsize = len(tree.nodes())
print(matsize)
print('nodes')

with h5py.File('./UKdata/aln.h5', 'r') as hf:
    align = hf['MSA2array'][:]
print(align.shape)

msa = AlignIO.read('./gisaid/msa_0612.lenfilter.fasta' , format = 'fasta')
def clipID(ID):
    return ''.join( [ s +'|' for s in str(ID).split('|')[:-1] ])[:-1].replace('_',' ')
IDs = {i:clipID(rec.id) for i,rec in enumerate(msa)}
IDindex = dict(zip( IDs.values() , IDs.keys() ) )
print( [(t,IDindex[t]) for t in list(IDindex.keys())[0:10]] )

#sites = { col: Counter( msa[:,col] ) for col in range(len(msa[1,:])) }
sites = {}
with h5py.File('./UKdata/aln.h5', 'r') as hf:
    align_array = hf['MSA2array'][:]
    #implement w np unique could be faster
    for col in range(align_array.shape[1]):
        if col% 1000  == 0:
            print(col)
        (unique, counts) = np.unique(align_array[:,col].ravel() , return_counts=True)
        sites.update({col:dict(zip(list(unique), list(counts)))})

informativesites = [ s for s in sites if len(set( sites[s].keys()) -set([b'-',b'N']) ) > 2  ]


#place an event on a column with multiple symbols
allowed_symbols = { b'A', b'C', b'G' , b'T' , b'-', b'N'}

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
                node.event = 1
        else:
            node.event = 0
        for child in node.child_nodes():
            if child.char is None:
                process_node_smallpars_2(child)
def calculate_small_parsimony( t, aln_column , row_index , verbose  = False ):
    missing = 0
    #assign leaf values
    for l in t.leaf_nodes():
        l.event = 0
        try:
            char = aln_column[row_index[str(l.taxon).replace("'" , '' )]]
        except KeyError:
            missing += 1
            char = None
        if char in allowed_symbols:
            l.symbols = { char }
        else:
            #ambiguous leaf with N
            l.symbols =  allowed_symbols
        l.scores = { c:0 if c in l.symbols else 10**10 for c in allowed_symbols }
        l.char = min(l.scores, key=l.scores.get)
    process_node_smallpars_1(t.seed_node)

    #down
    process_node_smallpars_2(t.seed_node)
    if verbose == True:
        print('done init')
        print(missing)
        print('missing in aln')
        for n in t.nodes()[0:5]:
            print(n)
            print(n.symbols)
            print(n.scores)
            print(n.char)
            print(n.event)
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

def mat_creator(retq,matsize,iolock):
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
        #import pdb; pdb.set_trace()
        M1[np.nonzero(events)[0] ,col] = 1
        if count == 0:
            with iolock:
                print(np.nonzero(events))

        if time.time()-t0> 120 :
            t0 = time.time()
            with iolock:
                print('saving')
                print(col)
                with open( 'coevmat.pkl' , 'wb') as coevout:
                    coevout.write(pickle.dumps(M1))
                print('done')

def main():
    NCORE = 60

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
    p = mp.Process(target=mat_creator, args=(retq,matsize, iolock))
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
