#!/usr/bin/env python
# coding: utf-8

# In[14]:


import sparse


import subprocess
import shlex
import pandas as pd
import numpy as np

import glob
import os
import wget
import requests
import glob
import time
import dask
import h5py

from Bio import SeqIO
import pickle
import sys
sys.setrecursionlimit( 10 **9 )
from sklearn.cluster import *
from sklearn.metrics import roc_curve , precision_recall_curve , auc

import scipy
import copy
from numpy import linalg as LA
from matplotlib import pyplot as plt
import random
import h5py
import itertools
import dendropy
import os
import psutil
import seaborn as sns


overwrite = False
jk_iterations = 5
os.environ['MKL_ENABLE_INSTRUCTIONS'] = 'AVX2'


# In[15]:


filedir = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/'

files = glob.glob( '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/apr_4_2022/mmsa_2022-04-04/2022-04-04_masked.fa_02022-06-14T07:53:44.200443small_test_coevmats.pkl')
print(files)
#treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
#alnfile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln'
alnfile = filedir + 'datasets/covid_data/msa_0730/msa_0730.fasta'

treefile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/feb_2021/GISAID-hCoV-19-phylogeny-2021-02-21/global.tree'
modeldir = filedir+'datasets/covid_data/structs/'

alnh5 = alnfile+'.h5'
#ts = '2021-08-08T11:16:34.358764'
#ts = '2021-08-08T14:37:59.736512'
#events = alnfile+'*'+ts+'*'
eventmats = files


# In[16]:


nucleotide_mutation = None
AA_mutation = None
for mat in eventmats:
    with open( mat , 'rb') as pklin:
        mats = pickle.loads(pklin.read())
        print(mats)
        if AA_mutation is None:
            nucleotide_mutation = mats[1]
            AA_mutation = mats[0]
        else:
            nucleotide_mutation += mats[1]
            AA_mutation += mats[0]
print(nucleotide_mutation)
print(AA_mutation)


# In[17]:


tree = dendropy.Tree.get( path=treefile, schema='newick' )


# In[18]:


#setup the internal nodes for the fitch algo
for i,l in enumerate(tree.nodes()):
    l.event = {}
    l.scores = {}
    l.symbols = None
    l.char= None
    l.matrow = i
    


# In[19]:



if os.path.exists( alnfile + '_IDs.pkl'):
    with open( alnfile + '_IDs.pkl' , 'rb') as idxin:
        IDindex = pickle.loads(idxin.read())
    IDs = dict(zip( IDindex.values() , IDindex.keys() ) )
else:
    
    msa = SeqIO.parse(alnfile , format = 'fasta')
    def clipID(ID):
        return ''.join( [ s +'|' for s in str(ID).split('|')[:-1] ])[:-1].replace('_',' ') 
    IDs = {i:rec.id for i,rec in enumerate(msa)}
    IDindex = dict(zip( IDs.values() , IDs.keys() ) )
    print( [(t,IDindex[t]) for t in list(IDindex.keys())[0:10]] )
    with open( alnfile + '_IDs.pkl' , 'wb') as idxout:
        idxout.write(pickle.dumps(IDindex))


# In[20]:


AAmat = AA_mutation.sum(axis = 2).to_scipy_sparse()





NT_mat = nucleotide_mutation.sum(axis = 2).to_scipy_sparse()


# In[24]:


for i,n in enumerate(tree.nodes()):
    n.matrow = i
    n.symbols = None
    n.scores = None
    n.event = None
    n.char = None

matsize = len(tree.nodes())
print(matsize)
print('nodes')
#blur w connectivity mat
blurfactor =  .25
connectmat = scipy.sparse.csc_matrix((len(tree.nodes()), len(tree.nodes() ) ) )
index = np.array([ [n.matrow, c.matrow ] for n in tree.nodes() for c in n.child_nodes()])
lengths = np.array([ c.edge_length for n in tree.nodes() for c in n.child_nodes()])
total_len = np.sum(lengths)
connectmat[index[:,0],index[:,1]] = 1
connectmat[index[:,1],index[:,0]] = 1


diag = [ i for i in range(connectmat.shape[0])]
connectmat[diag,diag] = 1
connectmat = scipy.sparse.coo_matrix(connectmat)




struct_annot = alnfile+'struct_blastout.txt'

proteo_annot  = alnfile +'annotation.csv'


annotations = pd.read_csv(  struct_annot , header = None )
annotationp = pd.read_csv(  proteo_annot , header = 0 )

annotationsp = annotationp.drop(columns = ['prots', 'genes', 'Unnamed: 0'])

annotations.columns = [ 'qseqid' , 'sseqid' , 'qlen' ,  'slen' , 'qstart' , 'qend' ,  'qframe' , 'evalue' ] 

annotations['struct'] = annotations.sseqid.map( lambda x : x.split(':')[0].split('/')[-1].replace('.pdb','') )
annotations['chain'] = annotations.sseqid.map( lambda x : x.split(':')[1] )
print(annotations)
print(annotationp)
annotation = pd.concat([annotations, annotationp])
annotation = annotation[ annotation['evalue'] < 10**-3 ]
annotation['coords'] = list(zip( annotation.qstart , annotation.qend ) )
annotation.sort_values('qstart')
print(annotation)
sub = annotation.drop_duplicates(subset=['coords'])
codons = []
for idx, row in annotation.iterrows():
    codons+= [ c for c in range(row.qstart-1,row.qend-1,3)]
codons = list( set(codons) )
codons = sorted(codons)




print('flatten matrix ' ) 
AA_mat = AA_mutation.sum(axis = 2)
NT_mat = nucleotide_mutation.sum(axis = 2)




###rewrite codon compilation for sparse
import sparse
def retcodons( AAmat , NTmat,codons, verbose = True ):
    #aa mutations for each pos
    AAmat_sub = sparse.stack(  [ AAmat[ : , codon:codon+2 ].sum(axis = 1) for codon in codons  ] , axis = 1 )    
    #add the frames for each in a stack 
    if verbose == True:
        print('done AA')
    NTmat_sub = sparse.stack([  sparse.stack( [ NTmat[:, codon + frame ] for frame in [0,1,2] ] , axis = 1 )  for codon in codons ] , axis = 1)
    if verbose == True:
        print('done NT')
    return AAmat_sub , NTmat_sub
codonmat =  retcodons( AA_mat , NT_mat, codons )
print(codonmat)


print(annotation)
structmats = {}
start_stop ={}
print(nucleotide_mutation)
print(AA_mutation)
from datasketch import MinHashLSHForest, WeightedMinHash , MinHash
#Ntaxa known. avoid allvsall
from datasketch import WeightedMinHashGenerator
import dask.array as da 
# WeightedMinHashGenerator requires dimension as the first argument
#flatten aa mat
print('done' )
def blur_cols( mat , blurmat , niter = 10 , distributed = True ):
    for i in range(niter):
        print(i)
        mat += blurmat.dot(mat)
    return mat
#transform dask array
#AAmat_dask = da.from_array(AA_mat)
#blurmat = da.from_array(connectmat)
print('blurring')
AA_blurmmat = blur_cols(codonmat[0].to_scipy_sparse(), connectmat, 20 )
NT_blurmmat = blur_cols(codonmat[1].sum(axis = 2).to_scipy_sparse(), connectmat, 10 )
# In[32]:
print(AA_blurmmat)




import h5py
import multiprocessing as mp

def hash_cols( cols  ):
    for dset in cols:
        col = cols[dset]
        nz = col.nonzero()[0]
        mh = MinHash(num_perm=256 , seed = 0 )
        [mh.update( (str(e)+dset).encode() ) for e in list( nz ) ]
    return mh
print('strating pool')
pool = mp.Pool(processes=mp.cpu_count()-1 )
print('done')

batch = mp.cpu_count()-1
res = []
thiscols = []
forest = MinHashLSHForest(num_perm=256)
print('starting')
mapper = {}
revmapper = {}
count = 0
with h5py.File("hashsigs_both.hdf5", "w") as f:
    dset = f.create_dataset("signatures", maxshape = (None, 256), shape=(len(codons) ,256) , dtype='int64')
    print('created dataset')
    for i,col in enumerate(codons):
        if i % 100 == 0:
            print(i, i / len(codons))
        res.append( pool.apply_async(hash_cols , [  { 'AA':AA_blurmmat[:,i] , 
            'NT':NT_blurmmat[:,i]    } ] ) )
        thiscols.append(col)
        if len(res)> batch :
            res = [ r.get() for r in res ]
            for n,col in enumerate(res):
                forest.add(thiscols[n], res[n])
                if dset.shape[0]< count-1:
                    dset.resize((count+10 ,256 ))  
                    print(dset.shape)
                dset[count,:] = res[n].digest()
                mapper[count] = thiscols[n]
                revmapper[ thiscols[n] ] =count 
                count +=1
                f.flush()
                
            res = []
            thiscols = []
    print(dset,dset[0:10],f)
pool.close()
with open('lshforest_both.pkl' , 'wb') as forestout:
    forestout.write(pickle.dumps(forest))

mapper = {}
for i,col in enumerate(codons):
    mapper[i]= col
revmapper = dict( zip( mapper.values(), mapper.keys()))
with open('lsh_codonmapper.pkl' , 'wb') as forestout:
    forestout.write(pickle.dumps(revmapper))


