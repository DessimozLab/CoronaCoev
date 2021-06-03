import numpy as np
from collections import Counter
import multiprocessing as mp
import pandas as pd
from Bio import SeqIO
import pandas as pd
import time
import h5py
import dendropy
from Bio import AlignIO , SeqIO , Seq
from scipy import sparse
import sparse as sparseND
import pickle
import copy
import sys
import os
import itertools


from dask.distributed import Client, Variable , Queue , Lock
from dask.distributed import  utils_perf
import gc

import dask.bag as db
import dask.array as da
import dask.dataframe as dd

from dask.delayed import delayed


#throttled = utils_perf.ThrottledGC()
#gc.collect = throttled.collect

if __name__ == '__main__':

    sys.setrecursionlimit( 10 **8 )
    runName = 'sparsemat_AAtransition'
    #number of cores to use

    distributed_computation = True

    if distributed_computation == True:
        NCORE = 100
        ncpu = 10

        print('deploying cluster')
        from dask_jobqueue import SLURMCluster

        cluster = SLURMCluster(
            n_workers = 25,
            cores=ncpu,
            memory="25 GB"
        )

        #cluster.adapt(minimum = 50,  maximum=NCORE)
        print(cluster.dashboard_link)
        client = Client(cluster)
        print('cluster deploy sleep')
        #wait for cluster deploy
        time.sleep( 5 )
        print('done')
    else:
        NCORE = 10
        ncpu = 1
        print('testing')
        client = Client()

    #fraction of genomes to remove if jackknifing
    bootstrap = .2
    #number of replicates
    bootstrap_replicates = 50
    restart = None
    nucleotides_only = False
    #keep track of transitions and not just events as binary
    transition_matrices = True
    future_clean_trigger = 100
    start_worker_trigger = 100

    treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
    alnfile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln'

    #treefile = '../validation_data/16s/16s_salaminWstruct_aln.fasta.treefile'
    #alnfile = '../validation_data/16s/16s_salaminWstruct_aln.fasta'

    #treefile = '../validation_data/dengue/dengue_all.aln.fasta.treefile'
    #alnfile = '../validation_data/dengue/dengue_all.aln.fasta'


    #create distributed cluster


    #use blast based annotation to assign codons to column ranges
    allowed_symbols = [ b'A', b'C', b'G' , b'T' ]
    allowed_transitions = [ c1+c2 for c1 in allowed_symbols for c2 in allowed_symbols  if c1!= c2]
    print('allowed transitions',allowed_transitions)


    transition_dict = {  c : i  for i,c in enumerate( allowed_transitions )  }
    rev_transition_dict= dict( zip(transition_dict.values(), transition_dict.keys()))
    allowed_symbols = set(allowed_symbols)

    print('transition dict', transition_dict)
    ProteinAlphabet = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
    allowed_AA_transitions = [ c1+c2 for c1 in ProteinAlphabet for c2 in ProteinAlphabet  if c1!= c2]
    print(allowed_AA_transitions[0:100] , '...etc...')

    transitiondict_AA = {  c : i  for i,c in enumerate( allowed_AA_transitions )  }
    rev_transitiondict_AA = dict( zip(transitiondict_AA.values(), transitiondict_AA.keys()))


    tree = dendropy.Tree.get(
        path=treefile,
        schema='newick')
    msa = AlignIO.read(alnfile , format = 'fasta')
    if os.path.exists(alnfile +'.h5'):
        pass#implement w np unique could be faster
    else:
        print('aln2numpy ')
        align_array = np.array([ list(rec.upper())  for rec in msa], np.character)
        print('done')
        print('dumping to hdf5')
        with h5py.File(alnfile +'.h5', 'w' ) as hf:
            hf.create_dataset("MSA2array",  data=align_array)
        print('done')



    if nucleotides_only == False:
        #use blast based annotation
        annotation = pd.read_csv( alnfile +'annotation.csv' )
    else:
        #just seperate sequence into dummy codons
        #indexing starts at 1 for blast
        print('using dummy annot')
        with h5py.File(alnfile +'.h5', 'r') as hf:
            align_array = hf['MSA2array']
            print('array shape' ,align_array.shape)
            dummy_annot = {'dummy_gene': { 'qstart':1 , 'qend':align_array.shape[1]-1 , 'evalue':0  }}
            annotation = pd.DataFrame.from_dict( dummy_annot , orient = 'index')
    print('selecting informative sites')


    def retcounts( row ):
        return  np.unique( row , return_counts=True)

    with h5py.File(alnfile +'.h5', 'r') as hf:
        align_array = hf['MSA2array'][:]
        array = da.from_array(align_array)
        array = da.transpose(array)
        #create a df of the columns
        daskdf = dd.from_dask_array(array)
        daskdf['unique']= daskdf.apply( retcounts  , axis =1 , meta = (None, 'object') )


        res = list( daskdf['unique'].compute() )

        print('compiling sites')
        sites= { col : dict(zip(list(unique[0]), list(unique[1]))) for col,unique in enumerate(res) }
        informativesites = set([ s for s in sites if len( set( sites[s].keys()) -set([b'-',b'N']) ) > 1  ] )

    print('done')
    print('informative columns:' , len(informativesites))

    #associate informative sites to a codon
    codon_dict = {}
    print( 'grouping codons')
    for i,r in annotation.iterrows():
        #indexing starts at 1 for blast
        for j,codon in enumerate(range(r.qstart-1, r.qend-1 , 3 )):
            for nt in [codon,codon+ 1, codon+2]:
                if nt in informativesites:
                    if (codon,codon+2) not in codon_dict:
                        codon_dict[(codon,codon+2)] = (nt,)
                    else:
                        codon_dict[(codon,codon+2)]+= (nt,)
    print('done')
    codon_dict_rev = dict(zip ( codon_dict.values() , codon_dict.keys( ) ) )
    def clipID(ID):
        return ID.replace('|',' ').replace('_',' ').replace('/',' ').strip()

    print('preparing tree IDs')
    IDs = {i:clipID(rec.id) for i,rec in enumerate(msa)}
    #IDs = {i:rec.id for i,rec in enumerate(msa)}
    IDindex = dict(zip( IDs.values() , IDs.keys() ) )
    print( [(t,IDindex[t]) for t in list(IDindex.keys())[0:10]] )


    ####small parsimony functions ##########
    def process_node_smallpars_1(node):
        #go from leaves up and generate character sets
        if node.symbols is None:
            for child in node.child_nodes():
                if child.symbols is None:
                    process_node_smallpars_1(child)
            node.symbols = { }
            node.scores = { }
            for pos in [0,1,2]:
                symbols = set.intersection( * [ child.symbols[pos] for child in node.child_nodes( ) ] )
                if len(symbols) == 0:
                    symbols = set.union( * [ child.symbols[pos] for child in node.child_nodes( ) ] )
                node.symbols[pos] = symbols
                node.scores[pos] = { }
                for c in allowed_symbols:
                    if c not in node.symbols[pos]:
                        #add trnasition mat here if needed
                        score = min(  [ child.scores[pos][c] for child in node.child_nodes()])+1
                    else:
                        score = min(  [ child.scores[pos][c] for child in node.child_nodes() ] )
                    node.scores[pos][c] = score

    def process_node_smallpars_2(node , verbose = False):
        #assign the most parsimonious char from children
        if node.char is None:
            if node.parent_node:
                #node has parent
                node.char = {}
                node.event = {}
                node.eventype= {}
                node.AAevent = None
                for pos in [0,1,2]:
                    node.char[pos] = min(node.scores[pos], key=node.scores[pos].get)
                    if node.parent_node.char[pos] == node.char[pos]:
                        node.event[pos] = 0
                    else:
                        if node.scores[pos][node.parent_node.char[pos]] == node.scores[pos][node.char[pos]] :
                            node.char[pos] = node.parent_node.char[pos]
                            node.event[pos] = 0
                        else:
                            node.event[pos] = 1
                            node.eventype[pos] = transition_dict[node.parent_node.char[pos]+node.char[pos]]
                node.AA = str(Seq.Seq(b''.join([ node.char[pos] for pos in [0,1,2] ]).decode() ).translate())

                if node.AA != node.parent_node.AA and nucleotides_only == False:
                    if node.parent_node.AA+node.AA in transitiondict_AA:
                        node.AAevent = transitiondict_AA[node.parent_node.AA+node.AA]
                        if verbose == True:
                            print( node.parent_node.AA , ' -> ' ,  node.AA)
            else:
                #root node
                node.char = {}
                node.event= {}
                node.eventype = {}
                node.AAevent = 0
                for pos in [0,1,2]:
                    node.char[pos] = min(node.scores[pos], key=node.scores[pos].get)
                    node.event[pos] = 0
                if nucleotides_only == False:
                    node.AA = Seq.Seq(b''.join([ node.char[pos] for pos in [0,1,2] ]).decode() ).translate()
                else:
                    node.AA = 'G'
            #down one level
            for child in node.child_nodes():
                if child.char is None:
                    process_node_smallpars_2(child)

    def calculate_small_parsimony( df, tree , row_index , bootstrap , array_size , verbose  = False ):

        #df is 3 columns of a codons
        #setup the tree and matrix for each worker
        missing = 0
        sys.setrecursionlimit( 10 **8 )
        t = copy.deepcopy(pickle.loads(tree))
        #assign leaf values
        #repeat here for bootstrap
        if bootstrap is not None:
            #select portion of random genomes to take out
            del_genomes = set(np.random.randint( array_size , size= int( array_size *bootstrap) ) )
        else:
            del_genomes = set([])


        #codon position
        pos =0
        #change a subset of leaves to ambiguous characters
        for idx,row in df.iterrows():
            for l in t.leaf_nodes():
                #setup for small_pars1
                l.calc[pos] = True
                l.event[pos] = 0
                l.scores[pos] = { c:10**10 for c in allowed_symbols }
                l.symbols[pos] =  allowed_symbols

                if str(l.taxon).replace("'", '') in row_index:
                    char = row[ row_index[str(l.taxon).replace("'", '')] ]
                    if char.upper() in allowed_symbols:
                        l.symbols[pos] = { char }
                        l.scores[pos][char] = 0
                    elif row_index[str(l.taxon).replace("'", '')] in del_genomes:
                        #eliminate for bootstrap
                        l.symbols[pos] =  allowed_symbols
                else:
                    char = None
                    l.symbols[pos] =  allowed_symbols
                l.char[pos] = min(l.scores[pos], key=l.scores[pos].get)
            pos+=1
        #done tree init
        #up
        process_node_smallpars_1(t.seed_node)
        #down
        process_node_smallpars_2(t.seed_node)
        #collect events
        eventdict = {}
        AAeventindex = [ n.matrow for n in t.nodes() if n.AAevent  ]
        AAeventypes = [ n.AAevent for n in t.nodes() if n.AAevent  ]
        for pos in [0,1,2]:
            eventindex = [ n.matrow for n in t.nodes() if n.event[pos] > 0 ]
            eventtypes = [ n.eventype[pos] for n in t.nodes() if n.event[pos] > 0 ]
            if pos ==0:
                eventdict[pos] = { 'type': eventtypes , 'index' : eventindex , 'AAeventindex':AAeventindex , 'AAeventypes': AAeventypes }
            else:
                eventdict[pos] = { 'type': eventtypes , 'index' : eventindex , 'AAeventindex':[] , 'AAeventypes': [] }

        retdf = pd.DataFrame.from_dict(eventdict, orient = 'index' )
        return retdf

    def save_mats(count, runName, AA_mutation,nucleotide_mutation):
        print('saving')
        with open( runName + '_coevmat_AAmutations.pkl' , 'wb') as coevout:
            coevout.write(pickle.dumps((count,AA_mutation)))
        with open( runName + '_coevmat_nucleotidemutations.pkl' , 'wb') as coevout:
            coevout.write(pickle.dumps((count,nucleotide_mutation)))
        print('done saving')


    def collect_futures(  runName  , check_interval= 10 , save_interval = 60, nucleotides_only =False  ):
        AA_mutation = None
        nucleotide_mutation = None
        t0 = time.time()
        runtime = time.time()
        count = 0
        outq = Queue('outq')
        stopiter = Variable('stopiter')

        while stopiter == False:
            #wait a little while
            result = outq.get()
            #get next job completed
            result = future.result()
            column, eventdict , AAeventindex , AAeventypes= result
            #save each position to event mats
            for pos in [0,1,2]:
                col = column+pos
                eventindex = eventdict[pos]['index']
                eventtypes = eventdict[pos]['type']
                if len(eventindex)>0:
                    if nucleotide_mutation:
                        nucleotide_mutation  += sparseND.COO( coords =  ( eventindex  , np.ones(len(selectrows )) * col   , eventtypes ) , data = np.ones(len(eventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transition_dict) ),  dtype = np.int32   )
                    else:
                        nucleotide_mutation  =  sparseND.COO( coords = ( eventindex ,  np.ones(len(selectrows )) * col   , eventtypes ) , data = np.ones(len(eventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transition_dict) ),  dtype = np.int32   )
            if nucleotides_only == False:
                if AA_mutation:
                    AA_mutation  += sparseND.COO( coords =  (AAeventindex , np.ones(len(AAeventindex)) * column , AAeventypes ) , data = np.ones(len(AAeventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transitiondict_AA ) ) ,  dtype = np.int32  )
                else:
                    AA_mutation  = sparseND.COO( coords =  (AAeventindex , np.ones(len(AAeventindex)) * column , AAeventypes ) , data = np.ones(len(AAeventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transitiondict_AA ) )   ,  dtype = np.int32 )
            count +=1
            if time.time() - runtime > save_interval:
                print('saving', time.time()-t0)
                runtime = time.time()
                save_mats(count, runName, AA_mutation,nucleotide_mutation)
        #finish up
        print('FINAL SAVE !')

        for result in  queue.get( timeout=None, batch=True):
                #get next job completed
                column, eventdict , AAeventindex , AAeventypes= result
                #save each position to event mats
                for pos in [0,1,2]:
                    col = column+pos
                    eventindex = eventdict[pos]['index']
                    eventtypes = eventdict[pos]['type']
                    if len(eventindex)>0:
                        if nucleotide_mutation:
                            nucleotide_mutation  += sparseND.COO( coords =  ( eventindex  , np.ones(len(selectrows )) * col   , eventtypes ) , data = np.ones(len(eventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transition_dict) ),  dtype = np.int32    )
                        else:
                            nucleotide_mutation  =  sparseND.COO( coords = ( eventindex ,  np.ones(len(selectrows )) * col   , eventtypes ) , data = np.ones(len(eventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transition_dict) ),  dtype = np.int32   )
                if nucleotides_only == False:
                    if AA_mutation:
                        AA_mutation  += sparseND.COO( coords =  (AAeventindex , np.ones(len(AAeventindex)) * column , AAeventypes ) , data = np.ones(len(AAeventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transitiondict_AA ) ) ,  dtype = np.int32  )
                    else:
                        AA_mutation  = sparseND.COO( coords =  (AAeventindex , np.ones(len(AAeventindex)) * column , AAeventypes ) , data = np.ones(len(AAeventindex)  ,  ) , shape = (matsize[0] , matsize[1] ,len(transitiondict_AA ) )   ,  dtype = np.int32 )
                count +=1

        save_mats(count, runName, AA_mutation,nucleotide_mutation)
        finished = Variable('saverDone')
        finished.set(True)

        print('DONE ! ')
        return None



    #######start the sankof algo here #######################
    print('starting sankof')
    with h5py.File(alnfile +'.h5', 'r', libver='latest', swmr=True ) as hf:
        print('loading array')
        align_array = hf['MSA2array'][:]

        array = da.from_array(align_array)
        array = da.transpose(array)

        #create a df of the columns
        daskdf = dd.from_dask_array(array)
        print('done')
        row_index = IDindex
        keep_codons = []
        keep_positions = []
        count =0



        #init the blank tree
        for i,n in enumerate(tree.nodes()):
            n.matrow = i
            n.symbols = None
            n.scores = None
            n.event = None
            n.char = None
            n.eventype = None
            n.AAevent = 0

        for i,l in enumerate(tree.leaf_nodes()):
            l.event = {}
            l.scores = {}
            l.symbols = {}
            l.char= {}
            l.calc = {}
        remote_tree = client.scatter( pickle.dumps(tree) )



        for annot_index,annot_row in annotation.iterrows():
            #indexing starts at 1 for blast
            #####switch to sending the coordinates and masking for the matrix
            for j,codon in enumerate(range(annot_row.qstart-1, annot_row.qend-1 , 3 )):
                positions = []
                keep_codon = False
                for col in [codon, codon+1 , codon+2]:
                    if col in informativesites:
                        keep_codon = True
                        #just add the alignment character if it doesnt change.
                    positions.append( col )
                if keep_codon == True:
                    keep_codons += [count,count,count]
                    keep_positions += positions
                count+=1

        keep_positions.reverse()
        keep_codons.reverse()

        print(keep_positions[0:100])
        print( keep_codons[0:100] )

        print('selecting')

        daskdf= daskdf.loc[keep_positions]
        daskdf['codon']= keep_codon
        #daskdf = client.persist(daskdf)
        print('done')
        #apply sankof to keep_codons
        res =[]




        for k in range(bootstrap_replicates):
            print('replicate:' , k )
            res.append( daskdf.groupby('codon').apply( calculate_small_parsimony , tree= remote_tree , row_index = row_index , bootstrap = bootstrap , array_size= align_array.shape[0] ,
            meta={ 'type': 'object' , 'index' : 'object' , 'AAeventindex':'object' , 'AAeventypes': 'object' }  ).compute() )
            print('done')
        #make sparse matrices from delayed
