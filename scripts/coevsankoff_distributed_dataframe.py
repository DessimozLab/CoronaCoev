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



import dask
from dask.distributed import Client, Variable , Queue , Lock , LocalCluster
from dask.distributed import  utils_perf
import gc

import dask.bag as db
import dask.array as da
import dask.dataframe as dd

from dask.delayed import delayed

####small parsimony functions ##########
def process_node_smallpars_1(node):
    #go from leaves up and generate character sets
    if node.symbols is None:
        for child in node.child_nodes():
            if child.symbols is None:
                process_node_smallpars_1(child)
        node.symbols = { }
        node.scores = { }
        for pos in child.symbols:
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
            for pos in node.scores:
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

            if len(node.char) == 3:
                node.AA = str(Seq.Seq(b''.join([ node.char[pos] for pos in node.char ]).decode() ).translate())
            else:
                node.AA = None

            if node.AA and node.AA != node.parent_node.AA:
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
            for pos in node.scores:
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

def calculate_small_parsimony( df, tree , row_index , bootstrap , array_size , replicate , verbose  = False ):
    replicate = str(replicate)

    if df.iloc[0,0] =='foo':

        #eventdict[pos] = { replicate+'_type': eventtypes , 'index' : eventindex , 'AAeventindex':AAeventindex , 'AAeventypes': AAeventypes }

        if replicate == '0':
            eventdict = { k:{ replicate+'type': [] , replicate+'index' : [] , replicate+'AAeventindex':[] , replicate+'AAeventypes': [] , 'pos':0  } for k in range(len(df))  }
        else:
            eventdict = { k:{ replicate+'type': [] , replicate+'index' : [] , replicate+'AAeventindex':[] , replicate+'AAeventypes': [] } for k in range(len(df))  }
        ret = pd.DataFrame.from_dict( eventdict , orient = 'index' )

        return ret
    if len(df) == 0:
        if replicate == '0':
            eventdict = { 0:{ replicate+'type': [] , replicate+'index' : [] , replicate+'AAeventindex':[] , replicate+'AAeventypes': [] , 'pos':0  } }
        else:
            eventdict = { 0:{ replicate+'type': [] , replicate+'index' : [] , replicate+'AAeventindex':[] , replicate+'AAeventypes': [] }   }
        ret = pd.DataFrame.from_dict( eventdict , orient = 'index' )

        return ret
    if len(df) > 0 and len(df) != 3:
        print('err:' , df)
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
    #change a subset of leaves to ambiguous characters
    pos = 0
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
    for pos,i in enumerate(df.index):
        eventindex = [ n.matrow for n in t.nodes() if n.event[pos] > 0 ]
        eventtypes = [ n.eventype[pos] for n in t.nodes() if n.event[pos] > 0 ]
        if pos==0:
            eventdict[i] = { replicate+'type': eventtypes , replicate+'index' : eventindex , replicate+'AAeventindex':AAeventindex , replicate+'AAeventypes': AAeventypes  }
        else:
            eventdict[i] = { replicate+'type': eventtypes , replicate+'index' : eventindex , replicate+'AAeventindex':[] , replicate+'AAeventypes': [] }
        if replicate == '0':
            eventdict[i]['pos'] = pos

    retdf = pd.DataFrame.from_dict(eventdict, orient = 'index' )
    return retdf

def retsum(df):
    #testing functions
    print(df)
    if df.iloc[0,0] =='foo':
        print('dummy data')
        ret = pd.DataFrame.from_dict( { k:{ 'sum':0 } for k in range(len(df))  } , orient = 'index' )

    else:
        s = df['codon'].sum()
        ret = pd.DataFrame.from_dict( { k:{ 'sum':s } for k in range(len(df))  } , orient = 'index' )
    ret.index = df.index
    print(ret)
    return ret

def compute_matrices(  resdf  , k , matsize , transitionsNT = 12 , transitionsAA = 380 ):

    AA_mutation = None
    nucleotide_mutation = None
    count = 0
    for idx,row in resdf.iterrows():

        for replicate in range(int(k)+1):

            replicate = str(replicate)
            #get next job completed

            eventtypes , eventindex , AAeventindex , AAeventypes= row[[replicate+'type' , replicate+'index' , replicate+'AAeventindex' , replicate+'AAeventypes']]
            eventtypes , eventindex , AAeventindex , AAeventypes = [ list(a)  for a in [eventtypes , eventindex , AAeventindex , AAeventypes ] ]
            #save each position to event mats
            col = int(idx[1])
            if len(eventindex)>0:

                if nucleotide_mutation:
                    nucleotide_mutation  += sparseND.COO( coords =  ( eventindex  , [ col for i in range(len(eventindex)) ]  , eventtypes ) , data = np.ones( len(eventindex) ) , shape = (matsize[0] , matsize[1] , transitionsNT)  )
                else:
                    nucleotide_mutation  =  sparseND.COO( coords = ( eventindex ,  [ col for i in range(len(eventindex)) ]    , eventtypes ) , data = np.ones( len(eventindex) ) , shape = (matsize[0] , matsize[1] , transitionsNT )  )

            if len(AAeventindex)>0:
                if AA_mutation:

                    AA_mutation  += sparseND.COO( coords =  (AAeventindex ,  [ col for i in range(len(AAeventindex)) ]   , AAeventypes ) , data = np.ones(len(AAeventindex)  ) , shape = (matsize[0] , matsize[1] , transitionsAA )   )
                else:
                    AA_mutation  = sparseND.COO( coords =  (AAeventindex ,  [ col for i in range(len(AAeventindex)) ]  , AAeventypes ) , data = np.ones(len(AAeventindex)    ) , shape = (matsize[0] , matsize[1] , transitionsAA )   )

    return (nucleotide_mutation, AA_mutation)

def add_sparsemats(args1, args2):
    res=[]
    for i in [0,1]:
        if args1[i] is not None and args2[i]is not None:
            res.append(args1[i]+args2[i])
        elif args1[i] is not None:
            res.append(args1[i])
        elif args2[i] is not None:
            res.append(args2[i])
        else:
            res.append(None)
    return  res


#throttled = utils_perf.ThrottledGC()
#gc.collect = throttled.collect

if __name__ == '__main__':

    sys.setrecursionlimit( 10 **8 )
    runName = 'sparsemat_AAtransition'
    #number of cores to use

    distributed_computation = True

    if distributed_computation == True:
        NCORE = 100
        ncpu = 1

        print('deploying cluster')
        from dask_jobqueue import SLURMCluster
        cluster = SLURMCluster(
            walltime='12:00:00',
            n_workers = NCORE,
            cores=ncpu,
            #cores=20,
            processes=1,
            memory="8GB",
            env_extra=[
            'source /scratch/dmoi/miniconda/etc/profile.d/conda.sh',
            'conda activate ML'
            ]
        )

        #cluster.adapt(minimum=40, maximum=NCORE)
        #cluster.adapt(minimum = 50,  maximum=NCORE)
        print(cluster.dashboard_link)
        client = Client(cluster , timeout='45s')
        client.restart()
    else:
        NCORE = 10
        ncpu = 1
        print('testing')
        cluster = LocalCluster(n_workers = NCORE )
        #cluster.adapt(minimum = 50,  maximum=NCORE)
        print(cluster.dashboard_link)
        client = Client(cluster)



    #fraction of genomes to remove if jackknifing
    bootstrap = .2
    #number of replicates
    bootstrap_replicates = 50

    overwrite = False
    restart = None

    nucleotides_only = False

    blastpath = '/scratch/dmoi/software/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin/'


    #treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
    #alnfile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln'

    #treefile = '../validation_data/16s/16s_salaminWstruct_aln.fasta.treefile'
    #alnfile = '../validation_data/16s/16s_salaminWstruct_aln.fasta'

    #treefile = '../validation_data/dengue/dengue_all.aln.fasta.treefile'
    #alnfile = '../validation_data/dengue/dengue_all.aln.fasta'


    treefile = '/scratch/dmoi/datasets/covid_data/30_may/mmsa_2021-06-01/global.tree'
    alnfile = '/scratch/dmoi/datasets/covid_data/30_may/mmsa_2021-06-01/2021-06-01_masked.fa'

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
    print('allowed_AA_transitions:', allowed_AA_transitions[0:100] , '...etc...')

    transitiondict_AA = {  c : i  for i,c in enumerate( allowed_AA_transitions )  }
    rev_transitiondict_AA = dict( zip(transitiondict_AA.values(), transitiondict_AA.keys()))

    #### creat hdf5 aln ###
    if os.path.exists(alnfile +'.h5') and overwrite == False:
        pass
    else:
        print('aln2numpy ')


        def ret_msa(alnfile):
            with open( alnfile ,'r') as msa:
                ret = []
                s = None
                for l in msa:
                    if '>' not in l:
                        if s:
                            s+=l.strip()
                        else:
                            s = l.strip()
                    else:
                        if s:
                            ret.append(list(s.upper()))
                            s = None
                            if len(ret) > 1000:
                                yield np.array(ret,  np.dtype('S1') ).T
                                ret = []
                                s = None
                if len( ret ) > 0:
                    yield np.array(ret, np.dtype('S1')).T

        gen = ret_msa(alnfile)
        chunk = next(gen)
        dtype = chunk.dtype
        row_count = chunk.shape[1]
        print(chunk)
        maxshape = ( 30000, 1000000)
        with h5py.File(alnfile +'.h5', 'w' ) as f:
            # Initialize a resizable dataset to hold the output
            #maxshape = (None,None) + chunk.shape[1:]
            dset = f.create_dataset('MSA2array', shape=chunk.shape, maxshape=maxshape,
                                    chunks=chunk.shape, dtype=chunk.dtype)

            dset[0:chunk.shape[0], 0:chunk.shape[1]] = chunk

            for i,chunk in enumerate(gen):
                if i % 10 == 0:
                    print(dset.shape)

                # Resize the dataset to accommodate the next chunk of rows
                dset.resize(row_count + chunk.shape[1], axis=1)
                # Write the next chunk
                dset[0:chunk.shape[0],row_count:] = chunk
                row_count += chunk.shape[1]
        print('done')

    ###get nonstatic sites
    def retcounts( row ):
        return  np.unique( row , return_counts=True)

    print('finding informative aln columns')
    with h5py.File(alnfile +'.h5', 'r') as hf:
        align_array = hf['MSA2array']
        array = da.from_array(align_array , chunks = {0: 1000})
        #array = da.transpose(array)
        #create a df of the columns
        daskdf = dd.from_dask_array(array)
        daskdf = daskdf.repartition(npartitions= NCORE  )
        daskdf['unique']= daskdf.apply( retcounts  , axis =1 , meta = (None, 'object') )
        res = list( daskdf['unique'].compute())
        #daskdf = client.persist(daskdf)
        print('compiling sites')
        sites= { col : dict(zip(list(unique[0]), list(unique[1]))) for col,unique in enumerate(res) }
        informativesites = set([ s for s in sites if len( set( sites[s].keys()) -set([b'-',b'N']) ) > 1  ] )
        #daskdf = daskdf.repartition(npartitions= NCORE * ncpu )
    print('done')
    print('informative columns:' , len(informativesites))

    def clipID(ID):
        return ID.replace('|',' ').replace('_',' ').replace('/',' ').strip()

    print('reading tree')
    tree = dendropy.Tree.get(
        path=treefile,
        schema='newick')
    print('done')


    print('preparing tree IDs')
    msa = SeqIO.parse(alnfile , format = 'fasta')
    IDs = {i:clipID(rec.id) for i,rec in enumerate(msa)}
    #IDs = {i:rec.id for i,rec in enumerate(msa)}
    IDindex = dict(zip( IDs.values() , IDs.keys() ) )

    print( [(t,IDindex[t]) for t in list(IDindex.keys())[0:10]] )
    print('done')

    seq = IDs[0]
    print('reference seq chosen in aln for codons: ' , seq)

    if nucleotides_only == False:
        #use blast based annotation
        if os.path.exists(alnfile +'annotation.csv' ):
            annotation = pd.read_csv( alnfile +'annotation.csv' )
        else:
            with h5py.File(alnfile + '.h5', 'r') as hf:
                align_array = hf['MSA2array']
                aln_row = align_array[:,IDindex[seq]]

                nongap_cols = [ i for i,c in enumerate(list(aln_row)) if c != b'-' ]
                print('nongap:' , nongap_cols[0:100] , '...')
                submat_aln = align_array[0,nongap_cols ]
                qseq = b''.join(aln_row[nongap_cols])
                qfile = alnfile+'codon_geno.fasta'
                print('qseq' , qseq[0:100] , '...')
                with open(qfile , 'w') as geno_out:
                    geno_out.write((b'>testgeno\n'+qseq).decode())
                import subprocess
                import shlex

                def runblastx( qseq , blastpath = '', outannot = 'outannot.txt' , outfmt = None , db=None):
                    if outfmt is None:
                        outfmt = [ 'qseqid' , 'sseqid' , 'qlen' ,  'slen' , 'qstart' , 'qend' ,  'qframe' , 'evalue' ]
                        outfmt =  ' "10 ' + ''.join([fmt+ ' ' for fmt in outfmt]) + ' " '
                        print(outfmt)
                    args = blastpath + 'blastx -query '+ qfile + ' -db ' +db+' -outfmt' + outfmt + ' -out ' + outannot
                    p = subprocess.run( shlex.split(args) )
                    return p , outannot

                #prepare a ref proteome ahead of time...
                p,annot = runblastx(qfile , blastpath = blastpath , outannot= alnfile+'outannot.csv' , db = '/scratch/dmoi/projects/covid/validation_data/covid19/covidref_Geno.txt' )
                annotation = pd.read_csv( annot , header = None )
                annotation.columns = [ 'qseqid' , 'sseqid' , 'qlen' ,  'slen' , 'qstart' , 'qend' ,  'qframe' , 'evalue' ]
                annotation = annotation[ annotation['evalue'] < 10**-3 ]
                print(annotation)
                print(len(annotation), ' orfs detected')
                #select longest nice hit
                rows = []
                for ID in annotation.sseqid.unique():
                    #print(annotation[annotation.sseqid == ID ].iloc[0])
                    sub = annotation[annotation.sseqid == ID ]
                    rows.append(sub.index[0] )
                annotation = annotation.loc[rows]
                rows = []
                for ID in annotation.qstart.unique():
                    #print(annotation[annotation.sseqid == ID ].iloc[0])
                    sub = annotation[annotation.qstart == ID ]
                    rows.append(sub.index[0] )
                annotation = annotation.loc[rows]
                genes =  {}
                prots = {}
                for i,r in annotation.iterrows():
                    genes[i] = qseq[r.qstart-1:r.qend-1].decode()
                    #print(genes[i])
                    prots[i] = str(Seq.Seq( genes[i]).translate( ) )
                annotation = annotation.sort_values( ['qstart'] )
                annotation['prots'] = annotation.index.map(prots)
                annotation['genes'] = annotation.index.map(genes)
                annotation = annotation[annotation.slen < 2000]

                aln_regions = np.array(list(zip(list(annotation.qstart),list(annotation.qend))))
                aln_regions= aln_regions[1:,:]
                aln_len = np.array(list( annotation.qend - annotation.qstart))

                annotation = pd.DataFrame.sort_values(annotation, by='qstart')
                print(annotation)
                annotation.to_csv(alnfile + 'annotation.csv')

    else:
        #just seperate sequence into dummy codons
        #indexing starts at 1 for blast
        print('using dummy annot')
        with h5py.File(alnfile +'.h5', 'r') as hf:
            align_array = hf['MSA2array']
            print('array shape' ,align_array.shape)
            dummy_annot = {'dummy_gene': { 'qstart':1 , 'qend':align_array.shape[1]-1 , 'evalue':0  }}
            annotation = pd.DataFrame.from_dict( dummy_annot , orient = 'index')

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
    print('selecting informative sites')

    #######start the sankof algo here #######################
    print('starting sankof')
    #daskdf = pd.DataFrame(data = align_array.T)
    row_index = IDindex
    keep_codons = []
    keep_positions = []
    count =0
    daskdf = daskdf.dropna()

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
    row_index = client.scatter(row_index )


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

    print('selecting')
    daskdf= daskdf.loc[keep_positions]
    print('positions to analyze:', len(daskdf))
    daskdf['codon'] = dd.from_array(np.array(keep_codons))
    #daskdf['codon'] = np.array(keep_codons)
    print('done')

    #apply sankof to keep_codons
    daskdf.columns = daskdf.columns.map(lambda x: str(x) )
    ret = None
    for k in range(bootstrap_replicates):
        print('replicate:' , k )
        k = str(k)
        if ret is None:
            ret = daskdf.groupby('codon' ).apply( calculate_small_parsimony , tree= remote_tree , row_index = row_index ,  bootstrap = bootstrap , array_size= align_array.shape[0] , replicate = k  )
        else:
            m = daskdf.groupby('codon' ).apply( calculate_small_parsimony , tree= remote_tree , row_index = row_index ,  bootstrap = bootstrap , array_size= align_array.shape[0] , replicate = k  )
            ret = ret.merge(m , left_index = True, right_index = True )
        dfs = ret.to_delayed()
        matrices = [ delayed(compute_matrices)(df,k, retmatsize) for df in dfs ]
        #do binary sum

        L = matrices
        while len(L) > 1:
            new_L = []
            for i in range(0, len(L), 2):
                try:
                    lazy = delayed(add_sparsemats)(L[i], L[i + 1])  # add neighbors
                except:
                    lazy = L[i]
                new_L.append(lazy)
            L = new_L                       # swap old list for new
        sparsemats = dask.compute(L)
        with open( alnfile +'_' + k +'_coevmats.pkl' , 'wb') as coevout:
            print('saving',sparsemats)
            coevout.write(pickle.dumps(sparsemats))
            print('done')
    print('DONE!')
    #make sparse matrices from delayed
