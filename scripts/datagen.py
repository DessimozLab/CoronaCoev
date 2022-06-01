#generate data for graphnet
import pandas as pd
import numpy as np
import glob
import os
import wget
import requests
import glob

from Bio import SeqIO
import pickle
import sys
sys.setrecursionlimit( 10 **9 )

from sklearn.metrics import roc_curve , precision_recall_curve , auc

import scipy
import copy

import h5py
import itertools
import dendropy


import torch.nn.functional as F
import torch_geometric.transforms as T
from torch_geometric.nn import ChebConv
from torch_geometric.nn import  to_hetero
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.data import HeteroData ,InMemoryDataset
import copy
import time
import random
import sparse


overwrite = True
restrictAA = True
testgen = True
traingen = True
reload = False




#column to struct mapping
annotation = pd.read_csv( alnfile +'struct_annotation.csv'  )

#treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
#alnfile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln'
alnfile = filedir + 'datasets/covid_data/msa_0730/msa_0730.fasta'
treefile = filedir + 'datasets/covid_data/msa_0730/global.tree'
modeldir = filedir+'datasets/covid_data/structs/'

alnh5 = alnfile+'.h5'
#ts = '2021-08-08T11:16:34.358764'
ts = '2021-08-08T14:37:59.736512'
events = alnfile+'*'+ts+'*'
eventmats = glob.glob(events)


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


tree = dendropy.Tree.get( path=treefile, schema='newick')
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
connectmat = scipy.sparse.csc_matrix((len(tree.nodes()), len(tree.nodes() ) ) )
index = np.array([ [n.matrow, c.matrow ] for n in tree.nodes() for c in n.child_nodes()])
lengths = np.array([ c.edge_length for n in tree.nodes() for c in n.child_nodes()])
total_len = np.sum(lengths)
connectmat[index[:,0],index[:,1]] = 1
connectmat[index[:,1],index[:,0]] = 1
#connectmat = connectmat.todense()
diag = [ i for i in range(connectmat.shape[0])]
connectmat[diag,diag] = 1
#connectmat = connectmat.todense()
#connectmat = scipy.sparse.csc_matrix(connectmat)
#np.fill_diagonal(connectmat , 1)
connectmat = scipy.sparse.coo_matrix(connectmat)

#read in contact maps
with open(modeldir + 'contactmaps' , 'rb') as connectout:
    subthresh_thresh , subthresh_connected = pickle.loads(connectout.read())
#read in phylogeny features
with open('template_features.pkl' , 'rb') as template_in:
    connectmat_up, connectmat_down, connectmat_diag, template_features = pickle.loads(template_in.read())
#read in sector mats
with open('sectormat.pkl' , 'rb') as sector_in:
    sectormat = pickle.loads(sector_in.read())
os.environ['MKL_ENABLE_INSTRUCTIONS'] = 'AVX2'

#create reduced alphabet mapping
#AA transitions
murphy12 = [('L','V','I','M'), ('C'), ('A'), ('G'), ('S','T'), ('P'), ('F','Y'), ('W'), ('E','Q'), ('D','N'), ('K','R'), ('H') ]
murphy12 = { c:i for i,cset in enumerate(murphy12) for c in cset   }
print('murphy12',murphy12)

ProteinAlphabet = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
allowed_AA_transitions = [ c1+c2 for c1 in ProteinAlphabet for c2 in ProteinAlphabet  if c1!= c2]

print('ntransitions' ,len(set(allowed_AA_transitions)))
new_transitions = [ (murphy12[c1],murphy12[c2]) for c1 in ProteinAlphabet for c2 in ProteinAlphabet  if c1!= c2  ]
new_transitions = { tup:i for i,tup in enumerate(set(new_transitions)) }
print('nmurphy transitions' , len(set(new_transitions)))
print(allowed_AA_transitions[0:100] , '...etc...')
transitiondict_AA = {  c : i  for i,c in enumerate( allowed_AA_transitions )  }
rev_transitiondict_AA = dict( zip(transitiondict_AA.values(), transitiondict_AA.keys()))

#reduce alphabet to reduce AA transition dimensionality
def restrictAA_transitions(AAmat, rev_transitiondict_AA, murphy12 , new_transitions , restricted_transitions = 138 , verbose = True ):
    restricted = None
    for z in range(AAmat.shape[2]):
        
        transition = rev_transitiondict_AA[z]
        new_transition = new_transitions[ (murphy12[transition[0]] , murphy12[transition[1]] )]
        data = AAmat[:,:,z].data
        coords =  AAmat[:,:,z].coords
        if coords.shape[1]>0:    
            newcoords = np.vstack( [coords, np.ones((1,coords.shape[1]))*new_transition]).astype('int')
            if restricted is not None:
                restricted  += sparse.COO( coords =  newcoords , data = data
                                            , shape = (AAmat.shape[0] , AAmat.shape[1] , restricted_transitions )  )
            else:
                restricted  =  sparse.COO( coords = newcoords , data = data  
                                            , shape = (AAmat.shape[0] , AAmat.shape[1] , restricted_transitions )  )
    return restricted

#return codons of an amino acid sequence from annotation
def retcodons( AAmat , NTmat, qstart, qend, verbose = True ):
    #aa mutations for each pos
    AAmat_sub = sparse.stack(  [ AAmat[ : , codon:codon+2 , : ].sum(axis = 1) for codon in range(qstart-1, qend-1 , 3 )  ] , axis = 1 )    
    #add the frames for each in a stack 
    if verbose == True:
        print('done AA')
    NTmat_sub = sparse.stack([  sparse.stack( [ NTmat[:, codon + frame , : ] for frame in [0,1,2] ] , axis = 1 )  for codon in range(qstart-1, qend-1 , 3 ) ] , axis = 1)
    if verbose == True:
        print('done NT')
    return AAmat_sub , NTmat_sub

print(annotation)
structmats = {}
start_stop ={}
print(nucleotide_mutation)
print(AA_mutation)
for i,row  in annotation.iterrows():
    print(row.qstart)
    if row.struct not in structmats:
        structmats[row.struct]={}    
    if (row.qstart , row.qend) not in start_stop:
        structmats[row.struct][row.chain ] = retcodons(AA_mutation, nucleotide_mutation  , row.qstart , row.qend)
        if restrictAA == True:
            structmats[row.struct][row.chain ] = ( restrictAA_transitions(structmats[row.struct][row.chain ][0] , rev_transitiondict_AA, murphy12 , new_transitions , restricted_transitions = 138 ) , structmats[row.struct][row.chain ][1] )
        start_stop[(row.qstart , row.qend)] = structmats[row.struct ][ row.chain ]
    else: 
        structmats[row.struct][row.chain ] = start_stop[(row.qstart , row.qend)]
        


def tree2Single_sparse_graph_updown(tree):
    N = len(tree.nodes())
    #mimic the fitch algo
    #propagate up and down in separate graphs
    index_up = np.vstack([ [n.matrow, c.matrow ] for n in tree.nodes() for c in n.child_nodes()])
    index_down = np.vstack([ [c.matrow, n.matrow ] for n in tree.nodes() for c in n.child_nodes()])
    connectmat_up = scipy.sparse.lil_matrix(( N ,  N ) )
    connectmat_down = scipy.sparse.lil_matrix(( N ,  N ) )
    connectmat_up[index_up[:,0],index_up[:,1]] = 1 
    connectmat_down[index_down[:,0],index_down[:,1]] = 1 
    diag = [[n,n] for n in range(N)]
    connectmat_diag=scipy.sparse.lil_matrix(( N ,  N ) )
    connectmat_diag[diag,diag] = 1 
    ntime = np.array([ n.distance_from_root() for n in tree.nodes()])
    mtime = np.amax(ntime)
    ntime/=mtime
    levels = np.array([ n.level() for n in tree.nodes() ] , dtype='double')
    levels /= np.amax(levels)
    Norm_nchild= np.array( [ len(n.child_nodes()) for n in tree.nodes() ] ,dtype='double' )
    mchild =np.amax(Norm_nchild)
    Norm_nchild/=mchild 
    Norm_nsister= np.array( [ len(n.sister_nodes()) for n in tree.nodes() ] ,dtype='double' )
    msis =np.amax(Norm_nsister)
    Norm_nsister/=msis    
    template_features = np.stack([ntime ,  Norm_nchild , Norm_nsister ]).T    
    return connectmat_up, connectmat_down, connectmat_diag, template_features

def sparse2pairs(sparsemat, matrows = None):
    if matrows :
        sparsemat = sparsemat[matrows,:]
        sparsemat = sparsemat[:,matrows]
    sparsemat = scipy.sparse.find(sparsemat)
    return np.vstack([sparsemat[0],sparsemat[1]])



def ret_pytorch_sample(subfeatures, connect_up , connect_down , sub_diag ,  overview , overview_rev , intra , toss ):
    data = HeteroData()
    #add input data
    data['phylonodes_up'].x = torch.tensor( subfeatures )
    data['phylonodes_down'].x =torch.tensor( subfeatures )
    data['sectornode'].x =torch.tensor(  np.zeros((1,1)) )
    #up down fitch net
    data['phylonodes_up', 'phylolink_up', 'phylonodes_up'].edge_index = torch.tensor(connect_up ,  dtype=torch.long )
    data['phylonodes_down', 'phylolink_down', 'phylonodes_down'].edge_index = torch.tensor(connect_down ,  dtype=torch.long )             
    data['phylonodes_up', 'phylolink_up_down', 'phylonodes_down'].edge_index = torch.tensor( sub_diag ,  dtype=torch.long )
    data['phylonodes_down', 'phylolink_down_up', 'phylonodes_up'].edge_index = torch.tensor( sub_diag ,  dtype=torch.long )
    #pooling connections
    data['phylonodes_down', 'informs', 'sectornode'].edge_index = torch.tensor(overview ,  dtype=torch.long )
    data['phylonodes_up', 'informs', 'sectornode'].edge_index = torch.tensor(overview ,  dtype=torch.long )

    #pooling connections
    data['sectornode',  'informs', 'phylonodes_down' ].edge_index = torch.tensor(overview_rev,  dtype=torch.long )
    data['sectornode',  'informs', 'phylonodes_up'].edge_index = torch.tensor(overview_rev ,  dtype=torch.long )
    #categories are intra or interprotein contacts or nothing
    if intra == True:
        cateforical = np.array([0,1])
    else:
        cateforical = np.array([1,0])
    data['phylonodes_down'].y =torch.tensor( np.ones((subfeatures.shape[0],1) ) *toss ,  dtype=torch.long )
    data['phylonodes_up'].y =torch.tensor( np.ones((subfeatures.shape[0],1)) * toss ,  dtype=torch.long )
    #todo change to categorical
    data['sectornode'].y =torch.tensor(  np.ones((1,1))*toss  ,  dtype=torch.long )
    data = T.AddSelfLoops()(data)
    data = T.NormalizeFeatures()(data)
    return data


def create_data_updown_transitions(intra ,pairs, AA_tensor, NT_tensor, sectormat, connectmat_up, connectmat_down, connectmat_diag, 
 template_features , posi_percent = .5 , nsamples = 10 , min_nodes = 100,  q = None , iolock= None,  verbose = True, loop= True , sectors_chunk = True ):
    #upward and downward connected phylo nodes
    Nsectors = 0
    Nnodes = 0
    allcols =list( np.arange(NT_tensor.shape[1]) )
    while True:
        toss = scipy.stats.bernoulli.rvs(posi_percent, loc=0, size=1, random_state=None)
        label = np.ones((1,1))*toss
    
        if verbose == True:
            print('posi/nega',toss)
            print('posi/nega',label)

        if toss == 0:
            col1 = random.choice(allcols)
            col2 = col1
            while col1 == col2 and (col1,col2) not in pairs:
                col2 = random.choice(allcols)
            labels = np.zeros((template_features.shape[0],))           
        else:
            #positive sample
            pairtuple = random.choice(pairs)
            col1 = pairtuple[0]
            col2 = pairtuple[1]
        if verbose == True:
            print(col1, col2)
        nodeAAfeatures = sparse.stack( [AA_tensor[:,col1,:] ,AA_tensor[:,col2,:] ] , axis = 2  ).reshape((AA_tensor.shape[0],-1)).to_scipy_sparse()
        nt_cols = []
        for pos in [0,1,2]:
            nodeNTfeatures = sparse.stack( [NT_tensor[:,col1,pos,:] ,NT_tensor[:, col2,pos,:] ] , axis = 2)
            nt_cols.append(nodeNTfeatures.reshape((nodeNTfeatures.shape[0],-1)).to_scipy_sparse() )
        
        nodefeatures = scipy.sparse.hstack([nodeAAfeatures,scipy.sparse.coo_matrix(template_features)]+nt_cols)
        #slice the features into sectors and yield sectors
        if sectors_chunk == True:
            nodefeatures = scipy.sparse.lil_matrix(nodefeatures)
            count = 0
            sectors = list(np.arange(sectormat.shape[1])) 
            for sample in range(nsamples):
                sector = random.choice(sectors)
                rows = scipy.sparse.find(sectormat[:,sector])[0]
                if rows.shape[0]> min_nodes:
                    count +=1
                    if verbose == True:
                        print('rows',rows.shape)

                #node features
                subfeatures = nodefeatures[rows,:].todense()
                #phylonode connections
                sub_connect_up = connectmat_up[rows,:]
                sub_connect_up = sub_connect_up[:,rows]
                connect_up = sparse2pairs(sub_connect_up)
                sub_connect_down = connectmat_down[rows,:]
                sub_connect_down = sub_connect_down[:,rows]
                connect_down = sparse2pairs(sub_connect_down)
                sub_diag = connectmat_diag[rows,:]
                sub_diag = sub_diag[:,rows]
                sub_diag = sparse2pairs(sub_diag)
                #aggregator node
                overview = scipy.sparse.lil_matrix( (subfeatures.shape[0], 2 ) )
                overview[:,0] = 1
                overview_rev = sparse2pairs(overview.T)
                overview = sparse2pairs(overview)
                data = ret_pytorch_sample(subfeatures, connect_up , connect_down , sub_diag ,  overview , overview_rev , intra , toss )
                if q:
                    q.put(data)
                else:
                    yield data
        else:
            #yield the entire graph
            #not working with sparse yet...
            #node features
            #change to torch sparse
            values = nodefeatures.data
            indices = np.vstack((nodefeatures.row, nodefeatures.col))

            i = torch.LongTensor(indices)
            v = torch.FloatTensor(values)
            shape = nodefeatures.shape
            subfeatures = torch.sparse_coo_tensor(i, v, torch.Size(shape))
            subfeatures = subfeatures.coalesce()
            #phylonode connections
            sub_connect_up = connectmat_up
            connect_up = sparse2pairs(sub_connect_up)
            sub_connect_down = connectmat_down
            connect_down = sparse2pairs(sub_connect_down)
            sub_diag = connectmat_diag
            sub_diag = sparse2pairs(sub_diag)
            #aggregator node
            overview = scipy.sparse.lil_matrix( (subfeatures.shape[0], 2 ) )
            overview[:,0] = 1
            overview_rev = sparse2pairs(overview.T)
            overview = sparse2pairs(overview)
            data = ret_pytorch_sample(subfeatures, connect_up , connect_down , sub_diag ,  overview , overview_rev , intra , toss )
            if q:
                q.put(data)
            else:
                yield data
                



#reduce alphabet to reduce AA transition dimensionality
def restrictAA_transitions(AAmat, rev_transitiondict_AA, murphy12 , new_transitions , restricted_transitions = 138 , verbose = True ):
    restricted = None
    for z in range(AAmat.shape[2]):
        
        transition = rev_transitiondict_AA[z]
        new_transition = new_transitions[ (murphy12[transition[0]] , murphy12[transition[1]] )]
        data = AAmat[:,:,z].data
        coords =  AAmat[:,:,z].coords
        if coords.shape[1]>0:    
            newcoords = np.vstack( [coords, np.ones((1,coords.shape[1]))*new_transition]).astype('int')
            if restricted is not None:
                restricted  += sparse.COO( coords =  newcoords , data = data
                                            , shape = (AAmat.shape[0] , AAmat.shape[1] , restricted_transitions )  )
            else:
                restricted  =  sparse.COO( coords = newcoords , data = data  
                                            , shape = (AAmat.shape[0] , AAmat.shape[1] , restricted_transitions )  )
    return restricted




###rewrite codon compilation for sparse

def retcodons( AAmat , NTmat, qstart, qend, verbose = True ):
    #aa mutations for each pos
    AAmat_sub = sparse.stack(  [ AAmat[ : , codon:codon+2 , : ].sum(axis = 1) for codon in range(qstart-1, qend-1 , 3 )  ] , axis = 1 )    
    #add the frames for each in a stack 
    if verbose == True:
        print('done AA')
    NTmat_sub = sparse.stack([  sparse.stack( [ NTmat[:, codon + frame , : ] for frame in [0,1,2] ] , axis = 1 )  for codon in range(qstart-1, qend-1 , 3 ) ] , axis = 1)
    if verbose == True:
        print('done NT')
    return AAmat_sub , NTmat_sub


print(annotation)
structmats = {}
start_stop ={}
print(nucleotide_mutation)
print(AA_mutation)
for i,row  in annotation.iterrows():
    print(row.qstart)
    if row.struct not in structmats:
        structmats[row.struct]={}    
    if (row.qstart , row.qend) not in start_stop:
        structmats[row.struct][row.chain ] = retcodons(AA_mutation, nucleotide_mutation  , row.qstart , row.qend)
        if restrictAA == True:
            structmats[row.struct][row.chain ] = ( restrictAA_transitions(structmats[row.struct][row.chain ][0] , rev_transitiondict_AA, murphy12 , new_transitions , restricted_transitions = 138 ) , structmats[row.struct][row.chain ][1] )
        start_stop[(row.qstart , row.qend)] = structmats[row.struct ][ row.chain ]
    else: 
        structmats[row.struct][row.chain ] = start_stop[(row.qstart , row.qend)]




def gen_dataset(data_name , structmats , allpairs ,  reload = True  ):
    pairsample = 100
    count = 0
    samples = []
    print('writing dataset to' , data_name)
    for struct in subthresh_connected:
        for i,chains in enumerate(subthresh_connected[struct]):
            AA_tensor , NT_tensor = structmats[struct][chains]
            intra = len(chains)>1
            pairs = list(allpairs[struct][chains])
            contact_gen = create_data_updown_transitions(intra ,pairs, AA_tensor, NT_tensor, sectormat, connectmat_up, connectmat_down, 
                    connectmat_diag, template_features , posi_percent = .5 ,  q = None , iolock= None,  verbose = False, loop= True  )
            for data in contact_gen:
                samples.append(data)
                if len(samples) < 10:
                    print(data)
                if len(samples)%500==0 and len(samples)> 0:
                    print(len(samples))
                    count += 1
                    
                    with open(data_name + str(count) +'.pkl', 'wb' ) as datain:
                        datain.write(pickle.dumps(samples))
                        
                    samples = []



if traingen == True:
    trainfile = alnfile + 'train'
    gen_dataset(trainfile , structmats , allpairs ,  reload  )


if testgen == True:
    testfile = alnfile + 'test.pkl'
    gen_dataset(testfile , structmats , allpairs ,  reload )