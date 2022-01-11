import numpy as np

import pickle
import collections
import itertools
from scipy.stats import bernoulli
import math
import random
import pandas as pd
import scipy
import dendropy
import numpy as np
import datetime
import gzip
import sys
import copy
import os
import glob
import pdb
import random
import scipy.sparse as sparse


def yeildBags( mat ):
  nzrows = list( np.where( mat.sum( axis = 1 ) > 0 )[0])
  random.shuffle(nzrows)
  #shuffle row order

  for row in nzrows:
    index = np.argwhere( mat[row,:] )
    if index.shape[0]>1:
        index = list(index[:,1])
        yield  index , mat[row,:].todense() , row

def yield_nega( sampling , nnega , pow = .75 , simple = False):
    negatives = []
    terms = list(sampling.keys())
    while True:
        neg1 = [ random.choice(terms) for i in range(int(2 * nnega)) ]
        if simple == False:
            neg1 = [ n for n in neg1 if random.uniform(0, 1) < sampling[n]]
        if len(neg1)%2 == 1:
            neg1 = neg1[:-1]
        neg2 = np.array(neg1[0:int(len(neg1)/2)])
        neg1 = np.array(neg1[int(len(neg1)/2):])
        ret = np.vstack( [neg1,neg2] ).T
        yield ret

def yield_posi( sampling , mat , nposi  , iter_row= 10 , itermax = 100,  fancy = False , verbose = False):
    cols = itertools.cycle(yeildBags(mat ) )
    positives =None
    while True:
        index , matrow , rownum = next(cols)
        pairs = np.array([[c[0],c[1]] for i,c in enumerate( itertools.combinations( index , 2 ) )  ] )
        vals= dict( zip( index, [ p for p in matrow[0,index].flat] ) )
        pairvalues = np.vectorize(vals.get)(pairs)
        if fancy == True:
            #balance prob of event w sampling ( more for columns that dont show up often )
            #prob needs something a bit fancier here to reflect the overall amount of events
            matrow = np.multiply(  matrow , [ 1- sampling[i] for i in index ] )
            matrow /= np.amax(matrow)            
        else:
        #proportional to the highest proba
            matrow /= np.amax(matrow)
        #sample in each row a few times
        if verbose == True:
            print(rownum)
        for k in range(iter_row):
            iterations = 0
            while True:
                #sample as a function of the event confidence and sampling
                rand = np.random.uniform(low=0.0, high=1.0, size= pairvalues.shape[0] )
                ar1 = np.array( rand <  pairvalues[:,0]  , dtype = np.bool_ )
                ar2 = np.array( rand <  pairvalues[:,1] , dtype = np.bool_ )
                select = np.bitwise_and(ar1 ,ar2 )
                iterations +=1
                if select.any() or iterations > itermax:
                    break
            if select.any():
                posi = np.array(pairs)
                posi = posi[select,:]
                if verbose == True:
                    print(posi, posi.shape , select)
                if positives is not None:
                    positives = np.vstack( [posi, positives])
                else:
                    positives = np.array(posi)
                if positives.shape[0] > nposi:
                    yield positives
                    positives = None
        
def yield_samples( mat , sampling ,index , pow= .75 , split =.75 ,nsamples = 1000):

    #pick over represented columns more often in negatives
    nnega = int(nsamples*split)
    nposi = int(nsamples*(1-split))

    negagen = yield_nega(sampling , nnega , pow = .75)
    posigen = yield_posi( sampling , mat , nposi  )    
    while True:
        posi = next( posigen )
        nega = next( negagen )
        samples = np.vstack([posi, nega])
        samples = np.vectorize(index.get)(samples)
        labels = np.vstack( [np.ones((posi.shape[0],1)) , np.zeros((nega.shape[0],1))] )
        samples = np.hstack([ samples, labels ])
        x1 = samples[:,0]
        x2 = samples[:,1]
        y = samples[:,2]
        yield [x1,x2],y

def make_neg_sampling( mat , z =.01 , pow = .75 ):
    #proba of keeping a column in neg samples
    sumv = mat.sum(axis = 0) 
    nonzero = np.argwhere( sumv )[:,1]
    sumv = sumv[0,nonzero]
    sumv = np.multiply( np.power( sumv / z , pow ) +1 , z / sumv ) 
    sumv = sumv/np.amax(sumv)
    sampling = dict(zip( list(nonzero), list(sumv.flat)))
    index = dict( zip( sampling.keys(), range(len(sampling.keys()))))
    return  sampling , index


#load and add bootstraps
treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
alnfile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln'

#alnfile = '/scratch/dmoi/datasets/covid_data/msa_0730/msa_0730.fasta'
#treefile = '/scratch/dmoi/datasets/covid_data/msa_0730/global.tree'

alnh5 = alnfile+'.h5'

#modelfile = alnfile + 'embedding_newfile_TF.h5'
#modelfile = alnfile + 'embedding_simpleneg_TF.h5'

modelfile = alnfile + 'embedding_15TF.h5'


ts = '2021-08-08T11:16:34.358764'
#ts = '2021-08-08T14:37:59.736512'

overwrite_mat = True
retrain = False
preprocess = False
blur_iterations = 30
vector_dim = 15

if overwrite_mat == True or not os.path.exists( alnfile + '_blurmat.pkl'):
    #blur w connectivity mat
    if preprocess == False:
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
        AA_mutation/=len(eventmats)

        #load sparse ND array
        from scipy.sparse import coo_matrix
        print(AA_mutation.shape)
        for i in range( AA_mutation.shape[2] ):
            if i == 0:
                AAmat =  AA_mutation[:,:,i]
            else:
                AAmat +=  AA_mutation[:,:,i]
        AAmat = AAmat.to_scipy_sparse()

    else:
        print('loading filtered mat')
        with open(alnfile + 'AAmat_sum.pkl' , 'rb' ) as pklin:
            AAmat = pickle.loads( pklin.read() )

    AAmat = AAmat/AAmat.max()
    sys.setrecursionlimit(10**6)
    tree = dendropy.Tree.get(
        path=treefile,
        schema='newick')
    treelen = tree.length()
    treenodes = len(tree.nodes())

    print('nodes',treenodes)
    print('length',treelen)
    for i,n in enumerate(tree.nodes()):
        n.matrow = i
    matsize = len(tree.nodes())
    print(matsize)
    connectmat = scipy.sparse.csc_matrix((len(tree.nodes()), len(tree.nodes() ) ) )
    index = np.array([ [n.matrow, c.matrow ] for n in tree.nodes() for c in n.child_nodes()])
    #lengths = np.array([ c.edge_length for n in tree.nodes() for c in n.child_nodes()])
    #total_len = np.sum(lengths)
    blurfactor = .25
    connectmat[index[:,0],index[:,1]] = blurfactor
    connectmat[index[:,1],index[:,0]] = blurfactor
    connectmat = scipy.sparse.coo_matrix(connectmat)
    #blur matrix
    blurmat = copy.deepcopy(AAmat)
    for blur in range(blur_iterations):
        print(blur)
        blurmat += connectmat.dot(blurmat)
        print('done')

    #generate sampling
    sampling , index = make_neg_sampling( blurmat )

    with open( alnfile + '_blurmat.pkl' , 'wb' ) as matout:
        matout.write( pickle.dumps([sampling,index , blurmat])) 
else:
    with open( alnfile + '_blurmat.pkl' , 'rb' ) as matout:
        sampling,index ,blurmat = pickle.loads(matout.read() ) 


nterms = len(sampling)


samplegen = yield_samples( blurmat , sampling, index , pow= .75 , split =.5 ,nsamples = 10000)
print(nterms)





from keras.models import *
from keras.optimizers import *
from keras.layers import *
from keras.metrics import *
from keras.regularizers import *
from keras.callbacks import *
from tensorflow.compat.v1 import ConfigProto , Session
from tensorflow.compat.v1.keras import backend as K
config = ConfigProto()
config.gpu_options.allow_growth = True
config.gpu_options.per_process_gpu_memory_fraction= 0.95
K.set_session(Session(config=config) )


if retrain == False:
 
    #word2vec model to be trained
    input_target = Input((1,) , name='target_in')
    input_context = Input((1,) , name='context_in')
    embedding = Embedding( nterms, vector_dim, input_length=1, name='embedding' , embeddings_initializer='glorot_uniform' )
    target = embedding(input_target)
    target = Reshape((vector_dim, 1), name='target')(target)

    context = embedding(input_context)
    context = Reshape((vector_dim, 1) , name='context' )(context)

    dot_product = dot([target, context] , axes=1 , normalize = False)
    dot_product = Reshape((1,))(dot_product)


    # add the sigmoid output layer
    output = Dense(1, activation='sigmoid' , name = 'out')(dot_product)

    # create the primary training model
    #o = Adadelta(lr=1.0, rho=0.95)
    o = RMSprop(lr=0.025, rho=0.9)
    #o = Adagrad(lr=0.000075)
    model = Model(inputs=[input_target,input_context], outputs=[output])
    model.compile(loss='binary_crossentropy', optimizer=o , metrics = [ 'binary_accuracy'])
    embedder = Model( inputs=[input_target], outputs=[target] )

if retrain == True:
    model = print('Load the model..')
    model = load_model(modelfile)
    o = RMSprop(lr=0.025, rho=0.9)
    #o = Adadelta(lr=0.0001)
    #o = Nadam(lr=0.002, beta_1=0.9, beta_2=0.999)
    model.compile(loss='binary_crossentropy', optimizer=o , metrics = [ 'binary_accuracy'])


mc = ModelCheckpoint(modelfile, monitor = 'loss', mode = 'min', verbose = 1, save_best_only = False)
lr = ReduceLROnPlateau(monitor='loss', factor=0.5, patience= 20 , min_lr=0.000001 , verbose = 1)
tb = TensorBoard(log_dir='./logs_mk2',  update_freq='epoch')
#for sample in samplegen:

model.fit( samplegen , verbose=1, callbacks=[  tb, lr , mc ], steps_per_epoch = 100 , epochs = 1000000 )
model.save(modelfile)
