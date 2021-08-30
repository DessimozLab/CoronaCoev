import numpy as np
from keras.models import *
from keras.optimizers import *
from keras.layers import *
from keras.metrics import *
from keras.regularizers import *
from keras.callbacks import *
import tensorflow as tf
import pickle
import collections
import itertools
from scipy.stats import bernoulli
import math
import random
import pandas as pd
import numpy as np

import datetime
import gzip

import sys
import os

def yeildBags(mat, samples, pow):
  #
  nzrows = list( np.where( mat.sum( axis = 1 ) > 0 )[1] )
  for row in nzrows:
    index = np.argwhere( mat[row,:] )[1]
    yield  index , mat[row, index]  
  return cols


def prunesamples(samples , sampling ):

    #remove samples in probabilistic way
    #select overrepresented GO terms more often in negative sampling
    ar1 = np.array([ random.uniform(0, 1) > p  for p in [ sampling[ s[0] ] for s in samples ] ] , dtype = np.bool_ )
    ar2 = np.array([ random.uniform(0, 1) > p  for p in [ sampling[ s[0] ] for s in samples ] ] , dtype = np.bool_ )
    
    select = np.bitwise_and(ar1,ar2)
    samples = np.array(samples)[select,:]
    return samples

def gen_nega(sampling , nnega , pow = .75):
    while len(negatives)< int( 2* nnega ):
        neg1 = [ random.choice(terms) for i in range(nsamples) ]
        neg1 = [ n for n in neg1 if count[n]>0 and random.uniform(0, 1) < sampling[index[n]]**pow ]
        negatives +=neg1

def yield_posi( sampling , mat  ):

    positives = itertools.cycle(yeildBags(mat))
    while len(negatives)< int( 2* nnega ):


def makesamples( mat , sampling , pow= .75 , split =.5 ,nsamples = 1000):
    #generator function to loop through gaf generating samples...
    
    terms = list(index.keys())
    negatives = []
    #pick over represented columns more often in negatives
    

    



    for i,dat in enumerate(positives):
                samples = []
                maxiter = 100
                i = 0
                while len(samples) <1 and i < maxiter:
                    #favor less common words
                    samples =  [  index[s] for s in dat['GO'] if s in index and sampling[index[s]] > random.uniform(0,1)  and  count[s] > 20  ]
                    i += 1
                if i == maxiter:
                    samples =  [  index[s] for s in dat['GO'] if s in index ]

                posi = np.array([  [  c[0] ,c[1]  ]+ [1]  for c in itertools.combinations( samples , 2 )  if  c[0] != c[1] ] )
                samples =  np.vstack([posi,nega])


                if samples.shape[1]>1:
                    x1 = samples[:,0]
                    x2 = samples[:,1]
                    y = samples[:,2]
                    yield [x1,x2],y
            else:
                pass
        except ValueError:
            pass

#load and add bootstraps

treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
alnfile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln'
#alnfile = '/scratch/dmoi/datasets/covid_data/msa_0730/msa_0730.fasta'
#treefile = '/scratch/dmoi/datasets/covid_data/msa_0730/global.tree'
alnh5 = alnfile+'.h5'
ts = '2021-08-08T11:16:34.358764'
#ts = '2021-08-08T14:37:59.736512'


if overwrite_mat or not os.path.exists( alnfile + '_blurmat.pkl'):
    #blur w connectivity mat

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

    #load sparse ND array
    from scipy.sparse import coo_matrix
    print(AA_mutation.shape)
    for i in range( AA_mutation.shape[2] ):
        if i == 0:
            AAmat =  AA_mutation[:,:,i]
        else:
            AAmat +=  AA_mutation[:,:,i]
            

    AAmat = AAmat.to_scipy_sparse()

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
        n.symbols = None
        n.scores = None
        n.event = None
        n.char = None
    matsize = len(tree.nodes())
    print(matsize)
    connectmat = scipy.sparse.csc_matrix((len(tree.nodes()), len(tree.nodes() ) ) )
    index = np.array([ [n.matrow, c.matrow ] for n in tree.nodes() for c in n.child_nodes()])
    lengths = np.array([ c.edge_length for n in tree.nodes() for c in n.child_nodes()])
    total_len = np.sum(lengths)
    connectmat[index[:,0],index[:,1]] = 1
    connectmat[index[:,1],index[:,0]] = 1
    connectmat = scipy.sparse.coo_matrix(connectmat)
    #blur matrix
    for blur in iterations:
        blurmat += blurmat.dot(connectmat)
    #generate sampling
    with open( alnfile + '_blurmat.pkl' , 'wb' ) as matout:
        matout.write( pickle.dumps(sampling,mat)) 
else:
    with open( alnfile + '_blurmat.pkl' , 'rb' ) as matout:
        sampling,mat = pickle.loads(matout.read() ) 





config = tf.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction= 0.95
K.set_session(tf.Session(config=config))

#filter the index
nterms = len(index)

if retrain == False:
    #dimensionality of GO space
    vector_dim = 5
    #word2vec model to be trained
    input_target = Input((1,) , name='target_in')
    input_context = Input((1,) , name='context_in')

    embedding = Embedding(nterms, vector_dim, input_length=1, name='embedding' , embeddings_initializer='uniform', embeddings_regularizer= None,
     activity_regularizer=None, embeddings_constraint=None, mask_zero=False )

    target = embedding(input_target)
    target = Reshape((vector_dim, 1), name='target')(target)
    context = embedding(input_context)
    context = Reshape((vector_dim, 1) , name='context' )(context)
    similarity = dot([target, context], axes=0 , normalize = True )
    # now perform the dot product operation to get a similarity measure
    dot_product = dot([target, context] , axes=1)
    dot_product = Reshape((1,))(dot_product)

    # add the sigmoid output layer
    output = Dense(1, activation='sigmoid' , name = 'out')(dot_product)

    # create the primary training model
    #o = Adagrad(lr=0.001)
    #o = Adadelta(lr=1.0, rho=0.95)

    o = RMSprop(lr=0.025, rho=0.9)


    #o = Adagrad(lr=0.000075)

    model = Model(inputs=[input_target,input_context], outputs=[output])
    model.compile(loss='binary_crossentropy', optimizer=o , metrics = [ 'binary_accuracy'])
    embedder = Model( inputs=[input_target], output=target )
    validation_model = Model(input=[input_target, input_context], output=similarity)

    class SimilarityCallback:
        def run_sim(self):
            for i in range(valid_size):
                valid_word = reverse_dictionary[valid_examples[i]]
                top_k = 8  # number of nearest neighbors
                sim = self._get_sim(valid_examples[i])
                nearest = (-sim).argsort()[1:top_k + 1]
                log_str = 'Nearest to %s:' % valid_word
                for k in range(top_k):
                    close_word = reverse_dictionary[nearest[k]]
                    log_str = '%s %s,' % (log_str, close_word)
                print(log_str)

        @staticmethod
        def _get_sim(valid_word_idx):
            sim = np.zeros((vocab_size,))
            in_arr1 = np.zeros((1,))
            in_arr2 = np.zeros((1,))
            for i in range(vocab_size):
                in_arr1[0,] = valid_word_idx
                in_arr2[0,] = i
                out = validation_model.predict_on_batch([in_arr1, in_arr2])
                sim[i] = out
            return sim
    sim_cb = SimilarityCallback()
    ###modify this
    batchiter = 10000
    epochs = 100


if retrain == True:
    model = print('Load the model..')
    model = load_model(modelfile)
    #o = RMSprop(lr=0.0001, rho=0.9)
    o = Adadelta(lr=0.0001)
    #o = Nadam(lr=0.002, beta_1=0.9, beta_2=0.999)
    model.compile(loss='binary_crossentropy', optimizer=o , metrics = [ 'binary_accuracy'])

mc = ModelCheckpoint(modelfile, monitor = 'loss', mode = 'min', verbose = 1, save_best_only = False)
lr = ReduceLROnPlateau(monitor='loss', factor=0.5, patience= 20 , min_lr=0.000001 , verbose = 1)
tb = TensorBoard(log_dir='./logs', histogram_freq=0, batch_size=32, write_graph=True, write_grads=False, write_images=False, embeddings_freq=0, embeddings_layer_names=None, embeddings_metadata=None, embeddings_data=None, update_freq='epoch')

history = model.fit_generator(makesamples( gaf, obo ,sampling , c, index , Yancestors ), steps_per_epoch=10000, epochs=10000, verbose=1,
callbacks=[ mc , lr , tb  ], max_queue_size=10, workers=1, use_multiprocessing=False, shuffle=True, initial_epoch=0)

model.save(modelfile)
embedder.save(modelfile+'embedder')
