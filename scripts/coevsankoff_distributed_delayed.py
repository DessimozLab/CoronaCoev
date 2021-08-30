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
import pdb
import lzma

from dask.distributed import fire_and_forget
from dask.distributed import Client, Variable , Queue , Lock ,LocalCluster
from dask_jobqueue import SLURMCluster
from dask.distributed import  utils_perf
import gc
import dask
import dask.bag as db
import dask.array as da
import dask.dataframe as dd
from dask.delayed import delayed
from dask import delayed, compute

from datetime import datetime, timedelta


import warnings


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

@dask.delayed(pure=False)
def delayed_send(data):
    return lzma.compress(pickle.dumps(data))

@dask.delayed(pure=False)
def delayed_receive(data):
    return pickle.loads( lzma.decompress(data) )

@dask.delayed(pure=False)
def calculate_small_parsimony(tree , df  ,  bootstrap = None , position = 0 ):
    #df is 3 columns of a codons
    #setup the tree and matrix for each worker
    sys.setrecursionlimit( 10 **9 )
    t = pickle.loads(lzma.decompress(tree))
    df = pickle.loads(lzma.decompress(df))
    #assign leaf values
    #repeat here for bootstrap

    if bootstrap is not None :
        #select portion of random genomes to take out
        del_genomes = set(np.random.randint( len(t.leaf_nodes()) , size= int( len(t.leaf_nodes()) *bootstrap) ) )
    else:
        del_genomes = set([])

    #change a subset of leaves to ambiguous characters
    pos = 0
    for idx,alncol in df.iterrows():
        for l in t.leaf_nodes():
            #setup for small_pars1
            l.calc[pos] = True
            l.event[pos] = 0
            l.scores[pos] = { c:10**10 for c in allowed_symbols }
            l.symbols[pos] =  allowed_symbols
            if l.aln_row and l.aln_row not in del_genomes and l.aln_row in alncol:
                char = alncol[ l.aln_row ]
                if char.upper() in allowed_symbols:
                    l.symbols[pos] = { char }
                    l.scores[pos][char] = 0
                else:
                    char = None
                    l.symbols[pos] =  allowed_symbols
            else:
                char = None
                l.symbols[pos] =  allowed_symbols
            l.char[pos] = min(l.scores[pos], key=l.scores[pos].get)
        pos+=1
    #up
    process_node_smallpars_1(t.seed_node)
    #down
    process_node_smallpars_2(t.seed_node)
    #collect events
    eventdict = {}
    AAeventindex = [ n.matrow for n in t.nodes() if n.AAevent  ]
    AAeventypes = [ n.AAevent for n in t.nodes() if n.AAevent  ]
    for pos,i in enumerate(t.seed_node.char):
        eventindex = [ n.matrow for n in t.nodes() if n.event[pos] > 0 ]
        eventtypes = [ n.eventype[pos] for n in t.nodes() if n.event[pos] > 0 ]


        if pos==0:
            eventdict[i] = { 'type': eventtypes , 'index' : eventindex , 'AAeventindex':AAeventindex , 'AAeventypes': AAeventypes  }
        
        else:
            eventdict[i] = { 'type': eventtypes , 'index' : eventindex , 'AAeventindex':[] , 'AAeventypes': [] }


        eventdict[i]['codon_pos'] = pos
        eventdict[i]['column'] = position + pos
    retdf = pd.DataFrame.from_dict(eventdict, orient = 'index' )
    return retdf


###compute spares matrics from results #######################################################################

@dask.delayed(pure=False)
def compute_matrices(  resdf  ,  matsize , transitionsNT = 12 , transitionsAA = 380      ):
    count = 0
    AA_mutation = None
    nucleotide_mutation = None

    for idx,row in resdf.iterrows():
        #get next job completed
        eventtypes , eventindex , AAeventindex , AAeventypes= row[['type' , 'index' , 'AAeventindex' , 'AAeventypes']]
        eventtypes , eventindex , AAeventindex , AAeventypes = [ list(a)  for a in [eventtypes , eventindex , AAeventindex , AAeventypes ] ]
        #save each position to event mats
        col = row.column
        if nucleotide_mutation is not None:
            nucleotide_mutation  += sparseND.COO( coords =  ( eventindex  , [ col for i in range(len(eventindex)) ]  , eventtypes ) , data = np.ones( len(eventindex) ) , shape = (matsize[0] , matsize[1] , transitionsNT)  )
        else:
            nucleotide_mutation  =  sparseND.COO( coords = ( eventindex ,  [ col for i in range(len(eventindex)) ]    , eventtypes ) , data = np.ones( len(eventindex) ) , shape = (matsize[0] , matsize[1] , transitionsNT )  )
        if AA_mutation  is not None:
            AA_mutation  += sparseND.COO( coords =  ( AAeventindex ,  [ col for i in range(len(AAeventindex)) ]   , AAeventypes ) , data = np.ones(len(AAeventindex)  ) , shape = (matsize[0] , matsize[1] , transitionsAA )   )
        else:
            AA_mutation  = sparseND.COO( coords =  ( AAeventindex ,  [ col for i in range(len(AAeventindex)) ]  , AAeventypes ) , data = np.ones(len(AAeventindex)    ) , shape = (matsize[0] , matsize[1] , transitionsAA )   )
    return  nucleotide_mutation, AA_mutation


@dask.delayed(pure=False)
def retloc(df,loc):
    return df.loc[loc]


if __name__ == '__main__':

    sys.setrecursionlimit( 10 **8 )
    ts = datetime.utcnow().isoformat()
    print('timestamp',ts)


    tag = 'small_test'
    #number of cores to use

    distributed_computation = True

    #fraction of genomes to remove if jackknifing
    bootstrap = .25
    overwrite = False
    overwrite_annot = True
    overwrite_index = False
    #number of replicates
    bootstrap_replicates = 20
    restart = None
    nucleotides_only = False



    blastpath = '/scratch/dmoi/software/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin/'


    refproteodb = '/scratch/dmoi/datasets/covid_data/refproteome/covidrefproteome.fasta'
    #refproteodb = '/scratch/dmoi/datasets/covid_data/structs/covid_structs.fasta'

    alnfile = '/scratch/dmoi/datasets/covid_data/msa_0730/msa_0730.fasta'
    treefile = '/scratch/dmoi/datasets/covid_data/msa_0730/global.tree'


    #treefile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.aln.EPIID.treefile'
    #alnfile = '../validation_data/covid19/gisaid_hcov-2020_08_25.QC.NSoutlier.filter.deMaiomask.EPIID.aln'

    #treefile = '../validation_data/16s/16s_salaminWstruct_aln.fasta.treefile'
    #alnfile = '../validation_data/16s/16s_salaminWstruct_aln.fasta'

    #create distributed cluster
    #use blast based annotation to assign codons to column ranges


    allowed_symbols = set([ b'A', b'C', b'G' , b'T' ])
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

    
    ###########################################aln to hdf5 

    print('preparing tree IDs')
    if os.path.exists( alnfile + '_IDs.pkl') and overwrite_index==False:
        with open( alnfile + '_IDs.pkl' , 'rb') as idxin:
            IDindex = pickle.loads(idxin.read())
        IDs = dict(zip( IDindex.values() , IDindex.keys() ) )
        IDindex = dict(zip( IDs.values() , IDs.keys() ) )
    else:
        msa = SeqIO.parse(alnfile , format = 'fasta')
        #def clipID(ID):
        #    return ''.join( [ s +'|' for s in str(ID).split('|')[:-1] ])[:-1].replace('_',' ')
        def clipID(ID):
            if '|' in ID:
                ID =  ID.split('|')[1]
            if '_' in ID:
                ID = ID.replace( '_' , ' ')
            return ID
        IDs = {i:clipID(rec.id) for i,rec in enumerate(msa)}
        IDindex = dict(zip( IDs.values() , IDs.keys() ) )
        print( [(t,IDindex[t]) for t in list(IDindex.keys())[0:10]] )
        with open( alnfile + '_IDs.pkl' , 'wb') as idxout:
            idxout.write(pickle.dumps(IDindex))
    
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
                            if clipID(l) in found:
                                ret.append(list(s.upper()))
                                if len(ret) > 1000:
                                    yield np.array(ret,  np.dtype('S1') ).T
                                    ret = []
                                    s = None
                        s = None
                if len( ret ) > 0:
                    yield np.array(ret, np.dtype('S1')).T
        gen = ret_msa(alnfile)
        chunk = next(gen)
        dtype = chunk.dtype
        row_count = chunk.shape[1]
        maxshape = ( 60000, None)
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


    seq = IDs[100]
    print('reference seq chosen in aln for codons: ' , seq)

    print('loading aln')
    with h5py.File(alnfile +'.h5', 'r') as hf:
        align_array = hf['MSA2array'][:]
        print(align_array.shape)
        aln_row = align_array[:,IDindex[seq]]
        nongap_cols = [ i for i,c in enumerate(list(aln_row)) if c != b'-' ]        
        print(aln_row[0:100])
        print('nongap:' , nongap_cols[0:100] , '...')
        print('selecting non gap columns of aln from ref geno')
        align_array = align_array[nongap_cols, :] 
        print( align_array.shape )
        daskdf = pd.DataFrame( data = align_array )
        del align_array


    ###########################################init tree
    print('reading tree')
    tree = dendropy.Tree.get(
        path=treefile,
        schema='newick')
    print('done')

    print('init blank tree for sankof')
    #init the blank tree
    missing = []
    found = []
    for i,n in enumerate(tree.nodes()):
        n.matrow = i
        n.symbols = None
        n.scores = None
        n.event = None
        n.char = None
        n.eventype = None
        n.AAevent = 0
    for i,l in enumerate(tree.leaf_nodes()):        
        tax = str(l.taxon).replace("'" , '')
        try:
            l.aln_row = IDindex[tax]
            found.append(tax)
        except:
            l.aln_row = None
            missing.append(tax)
        l.event = {}
        l.scores = {}
        l.symbols = {}
        l.char= {}
        l.calc = {}
    print('missing leaves:' , len(missing))

    found = set(found)
    print(len(tree.leaf_nodes()))
    print('done')


    print('flashing up a dask cluster')
    if distributed_computation == True:


        NCORE = 10
        njobs = 20
        print('deploying cluster')


        cluster = SLURMCluster(
            walltime='24:00:00',
            n_workers = NCORE,
            cores=NCORE,
            processes = NCORE,
            interface='ib0',
            memory="150GB",
            env_extra=[
            'source /scratch/dmoi/miniconda/etc/profile.d/conda.sh',
            'conda activate ML'
            ],
            scheduler_options={'interface': 'ens2f0' }
        )
        print(cluster.job_script())
        #cluster.adapt(minimum=10, maximum=30)
        cluster.scale(jobs=20)
        time.sleep(5)
        print(cluster)
        print(cluster.dashboard_link)
        client = Client(cluster , timeout='450s' , set_as_default=True )
    else:
        NCORE = 50
        njobs = 1
        print('testing')
        cluster = LocalCluster(n_workers = NCORE )
        #cluster.adapt(minimum = 50,  maximum=NCORE)
        print(cluster.dashboard_link)
        client = Client(cluster)
    print('done')

    #keep non gap positions for one referene geno
    print('done')
    if nucleotides_only == False:
        #use blast based annotation
        if os.path.exists(alnfile +'annotation.csv' ) and overwrite_annot == False:
            annotation = pd.read_csv( alnfile +'annotation.csv' )
        else:
            qseq = b''.join(aln_row[nongap_cols] )
            qfile = alnfile+'codon_geno.fasta'
            print('qseq' , qseq[0:500] , '...')
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
            p,annot = runblastx(qfile , blastpath = blastpath , outannot= alnfile+'outannot.csv' , db = refproteodb )
            annotation = pd.read_csv( annot , header = None )
            annotation.columns = [ 'qseqid' , 'sseqid' , 'qlen' ,  'slen' , 'qstart' , 'qend' ,  'qframe' , 'evalue' ]
            annotation = annotation[ annotation['evalue'] < 10**-6 ]
            genes =  {}
            prots = {}
            positions = []            
            for i,r in annotation.iterrows():
                genes[i] = qseq[r.qstart-1:r.qend-1].decode()
                prots[i] = str(Seq.Seq( genes[i]).translate( ) )
                positions+=[ i for i in range(r.qstart-1,r.qend-1,3) ]            
            positions = set(positions)
            annotation = annotation.sort_values( ['qstart'] )
            annotation['prots'] = annotation.index.map(prots)
            annotation['genes'] = annotation.index.map(genes)
            annotation = pd.DataFrame.sort_values(annotation, by='qstart')
            print(len(annotation), ' prots detected')
            print(annotation)
            annotation.to_csv(alnfile + 'annotation.csv')
    else:
        #just seperate sequence into dummy codons
        #indexing starts at 1 for blast
        print('using dummy annot')
        dummy_annot = {'dummy_gene': { 'qstart':1 , 'qend':align_array.shape[0]-1 , 'evalue':0  }}
        annotation = pd.DataFrame.from_dict( dummy_annot , orient = 'index')

    #######start the sankof algo here #######################
    print('starting sankof')
    row_index = IDindex
    count =0
    remote_tree = client.scatter( lzma.compress( pickle.dumps(tree) ) , broadcast=True )
    retmatsize = ( len(tree.nodes()) ,daskdf.shape[0]  )
    coordinates = []


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        for k in range(bootstrap_replicates):
            inlist = []
            count = 0 
            #indexing starts at 1 for blast
            #####switch to sending the coordinates and masking for the matrix
            for i,codon in enumerate(positions):
                pos = [codon, codon+1 , codon+2]
                res = calculate_small_parsimony( remote_tree , delayed_send( daskdf.loc[pos] )  , bootstrap , codon )
                inlist.append( res )
                if len(inlist)== NCORE*njobs:
                    print('codon positions left to calclulate' , len(positions) - i )
                    delayed_mats = [ compute_matrices(df, retmatsize) for df in inlist ]
                    delayed_mats = dask.compute( * delayed_mats )

                    AAbag = dask.bag.from_sequence([ m[1] for m in delayed_mats ]) 
                    NTbag = dask.bag.from_sequence([ m[0] for m in delayed_mats ])
                    
                    if count ==0 :
                        matricesAA = dask.compute(AAbag.sum())[0]
                        matricesNT = dask.compute(NTbag.sum())[0]

                    else:
                        matricesAA += dask.compute(AAbag.sum())[0]
                        matricesNT += dask.compute(NTbag.sum())[0]

                    print(sparseND.argwhere(matricesAA))


                    count += 1
                    inlist =[]        
            
            #last batch
            print('codon positions to calclulate' , len(inlist) )
            delayed_mats = [ compute_matrices(df, retmatsize) for df in inlist ]
            delayed_mats = dask.compute( *delayed_mats )
            AAbag = dask.bag.from_sequence([  m[1] for m in delayed_mats ]) 
            NTbag = dask.bag.from_sequence([ m[0] for m in delayed_mats ])

            matricesAA += dask.compute(AAbag.sum())[0]
            matricesNT += dask.compute(NTbag.sum())[0]
            
            print(matricesAA)
            print(matricesNT)
            print(sparseND.argwhere(matricesAA))
            
            print('done')
            if bootstrap is None:
                filename = alnfile +'_' + str(k) +ts+ tag + '_coevmats.pkl' 
            else:
                filename = alnfile +'_' + str(k) +ts+ tag + '_BS_coevmats.pkl' 

            with open( filename , 'wb') as coevout:
                print(filename)
                print('saving',(matricesAA, matricesNT ) )
                coevout.write(pickle.dumps((matricesAA, matricesNT )))
                print('done')