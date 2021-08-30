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
    #ar1 = np.array([ random.uniform(0, 1) > p  for p in [ sampling[ s[0] ] for s in samples ] ] , dtype = np.bool_ )
    select = np.bitwise_and(ar1,ar2)
    samples = np.array(samples)[select,:]
    return samples


def makesamples( mat , sampling , pow):
    #generator function to loop through gaf generating samples...
    
    terms = list(index.keys())
    negatives = []
    print(sum(c.values()))
    pow = 1.5

    while len(negatives)< 10000000:
        neg1 = [ random.choice(terms) for i in range(100000) ]
        neg1 = [ n for n in neg1 if count[n]>0 and random.uniform(0, 1) < sampling[index[n]]**pow and count[n] > 50 ]
        negatives +=neg1
    #at least 100 annot in corpus
    thresh = 1

    infinite_gaf = itertools.cycle(gafreader)

    for i,dat in enumerate(infinite_gaf):
        try:
            #if i == 0:
            #    last =  [  index[s] for s in dat['GO'] if sampling[index[s]]**pow  < random.uniform(0,1) and count[index[s]] > thresh ]
            if len(dat['GO'])>1 and i > 0:
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
                nega = np.array([ [  random.choice(negatives) , random.choice(negatives) ] + [0]  for i in range(posi.shape[0])  ]  )

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

def blurmat( mat , iterations, set_const):


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
