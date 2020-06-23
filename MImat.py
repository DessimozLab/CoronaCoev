#calculat hug MI mat from aln

import numpy as np

from matplotlib import pyplot as plt


import itertools
import gzip
import multiprocessing as mp
import lzma


import csb
import pandas as pd
from Bio import SeqIO
import pandas as pd
from colour import Color
import pickle

import h5py

from Bio import AlignIO
import dendropy

recalc_sites = False
NCORE = 55

startk , startl = ( 6100, 8700 )


print('mod3')








tree = dendropy.Tree.get(


    path='./UKdata/cog_global_2020-05-08_tree.newick',
    schema='newick')

print(dir(tree))



for l in tree.leaf_nodes()[0:10]:
    print(l)

msa = AlignIO.read('./UKdata/cog_2020-05-08_alignment.fasta' , format = 'fasta')
print(msa)
align_array = np.array([list(rec) for rec in msa], np.character)
print(align_array)
meta_data = pd.read_csv('UKdata/cog_2020-05-08_metadata.csv')
print(meta_data)


"""
def set_sites( treenode , aln ):
    return set[  msa[seq , col] ) for col in range(len(msa)) ]
"""
if recalc_sites == True:
    from collections import Counter
    #sites = { col: Counter( msa[:,col] ) for col in range(len(msa[1,:])) }
    sites = {}
    for col in range(len(msa[1,:])):
        sites.update({col:Counter(msa[:,col])})
        if col% 1000  == 0:
            print(col)



    #sequences = { ID:Counter( msa[i,:]) for i,ID in enumerate(IDs) }
    sequences = {}
    for i in range(len(msa)):
        sequences.update({i:Counter(msa[i,:])})
        if i% 1000  == 0:
            print(i)

    #remove uniformative sites from counter
    #send the informative ones to MI analysis
    informativesites = [ s for s in sites if len(set( sites[s].keys()) -set(['-','N']) ) > 1  ]

    with open('./UKdata/site_seq_stats.pkl' , 'wb') as pickleout:
        pickleout.write(pickle.dumps([sites,informativesites, IDs , sequences]) )


else:
    pass

with open('./UKdata/site_seq_stats.pkl' , 'rb') as pickleout:
    sites,informativesites, IDs , sequences = pickle.loads(pickleout.read())
#find colum mutual info using compression distance
#find interprot and intraprot interactions

lzma_filters = my_filters = [
    {
      "id": lzma.FILTER_LZMA2,
      "preset": 9 | lzma.PRESET_EXTREME,
      "dict_size":len(msa[:,0]) * 40, # a big enough dictionary, but not more than needed, saves memory
      "lc": 3,
      "lp": 0,
      "pb": 0, # assume ascii
      "mode": lzma.MODE_NORMAL,
      "nice_len": 273,
      "mf": lzma.MF_BT4
    }
]

def clen(s):
    return len( lzma.compress(s, format=lzma.FORMAT_RAW, filters= lzma_filters) )
    #return len(gzip.compress(s))

def compress_dist(pargs):
    i,j,s1,s2 = pargs
    strings = [s1.ravel() , s2.ravel() , np.stack([s1,s2]).ravel() ]
    ls1 , ls2 , ls1_2 = map(clen, strings)
    return  (i,j,(ls1_2 - min(ls1,ls2) )/ max( ls1,ls2 ) )

def process(q,retq, iolock):
    #calculate compression distances
    with iolock:
        print('init worker')
    from time import sleep
    while True:
        stuff = q.get()
        if stuff is None:
            break
        retq.put(compress_dist(stuff))
    print('done')
def mat_creator(retq,matsize,iolock):
    with iolock:
        print('init matcreator')
    #collect distances and create final mat
    calculations = (matsize**2 - matsize) / 2
    distmat = np.zeros((matsize,matsize))
    count = 0

    with h5py.File('./UKdata/alnMI2.h5', 'a') as hf:
        try:
            hf.create_dataset("alnMI",  data=distmat)
        except:
            pass
        while True:
            r = retq.get()
            count+=1
            if r is None:
                break
            row,col,cdist = r
            hf['alnMI'][row,col] = cdist
            if count% 10000 == 0 :
                with iolock:
                    print(count/calculations)
                    print((row,col))
                    hf.flush()
    print('done saver')


if __name__ == '__main__':
    q = mp.Queue(maxsize=NCORE*5000)
    retq = mp.Queue(maxsize=NCORE*5000)
    iolock = mp.Lock()
    pool = mp.Pool(NCORE, initializer=process, initargs=(q,retq, iolock))
    p = mp.Process(target=mat_creator, args=(retq,len(msa[1,:]), iolock))
    p.start()
    for k,i in enumerate(informativesites):
        for l,j in enumerate(informativesites):
            if k < l and k > startk and l > startl :
                s1 = align_array[:,i].ravel()
                s2 = align_array[:,j].ravel()
                q.put( (i,j,s1,s2) )
                if k % 100 == 0 and l %100 == 0:
                    print((k,l))
                    print((i,j))

    for _ in range(NCORE):  # tell workers we're done
        q.put(None)

    #pool.join()

    #retq.put(None)
    #pool.close()
    p.join()
