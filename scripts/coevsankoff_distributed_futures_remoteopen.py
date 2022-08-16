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

from dask.distributed import fire_and_forget , secede, rejoin, as_completed

from dask.distributed import Client, Variable , Queue , Lock ,LocalCluster

from dask.distributed import worker_client

from dask_jobqueue import SLURMCluster
from dask.distributed import  utils_perf
import gc
import dask
import dask.bag as db
import dask.array as da
import dask.dataframe as dd
from dask.delayed import delayed
from dask import delayed, compute
import argparse
from datetime import datetime, timedelta
import warnings
from functools import wraps


sys.setrecursionlimit( 10 **9 )
parsbatch = 10

####small parsimony functions ##########
def process_node_smallpars_1(node, count = 0 ):
    sys.setrecursionlimit( 10 **9 )

    #go from leaves up and generate character sets
    if node.symbols is None:
        for child in node.child_nodes():
            if child.symbols is None:
                if count < parsbatch:
                    child = process_node_smallpars_1( child , count+1 )
                else:
                    with worker_client() as client:
                        f = client.submit( process_node_smallpars_1, child , 0 )
                        child = client.gather( f )
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
        return node


def process_node_smallpars_2(node , count =0 , verbose = False):

    sys.setrecursionlimit( 10 **9 )

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
            if count < parsbatch:
                child = process_node_smallpars_2( child , count + 1)
            else:
                with worker_client() as client:
                    if child.char is None:
                        f = client.submit( process_node_smallpars_2, child)
                        child = client.gather(f)
        return node


def calculate_small_parsimony(tree ,  posvec  ,  bootstrap = None , position = 0 , alnfile = None ):
    #df is 3 columns of a codons
    #setup the tree and matrix for each worker
    sys.setrecursionlimit( 10 **9 )
    with open(tree , 'rb') as treein:
        t = pickle.loads(treein.read())    
    
    #assign leaf values
    #repeat here for bootstrap

    if bootstrap is not None :
        #select portion of random genomes to take out
        del_genomes = set(np.random.randint( len(t.leaf_nodes()) , size= int( len(t.leaf_nodes()) *bootstrap) ) )
    else:
        del_genomes = set([])
    retdfs = []
    #change a subset of leaves to ambiguous characters
    for idx in posvec:
        with h5py.File(alnfile +'.h5', 'r') as hf:
            align_array = hf['MSA2array'] 
            df = align_array[idx  , : ]
        alnmax = df.shape[1]
        for codonpos in [0,1,2]:
            alncol = df[codonpos , : ]
            for l in t.leaf_nodes():
                #setup for small_pars1
                l.calc[codonpos] = True
                l.event[codonpos] = 0
                l.scores[codonpos] = { c:10**10 for c in allowed_symbols }
                l.symbols[codonpos] =  allowed_symbols
                if l.aln_row and l.aln_row < alnmax and l.aln_row not in del_genomes:
                    char = alncol[ l.aln_row ]
                    if char.upper() in allowed_symbols:
                        l.symbols[codonpos] = { char }
                        l.scores[codonpos][char] = 0
                    else:
                        char = None
                        l.symbols[codonpos] =  allowed_symbols
                else:
                    char = None
                    l.symbols[codonpos] =  allowed_symbols
                l.char[codonpos] = min(l.scores[codonpos], key=l.scores[codonpos].get)
        
        #calculate smallpars updown
        with worker_client() as client:
            t = client.submit( process_node_smallpars_1, t.seed_node ) 
            t = client.compute(t).result()
            t = client.submit( process_node_smallpars_2, t )
            t = client.compute(t).result()
        eventdict = {}
        AAeventindex = [ n.matrow for n in t.nodes() if n.AAevent  ]
        AAeventypes = [ n.AAevent for n in t.nodes() if n.AAevent  ]
        for codonpos,i in enumerate(t.seed_node.char):
            eventindex = [ n.matrow for n in t.nodes() if n.event[codonpos] > 0 ]
            eventtypes = [ n.eventype[codonpos] for n in t.nodes() if n.event[codonpos] > 0 ]
            if codonpos==0:
                eventdict[i] = { 'type': eventtypes , 'index' : eventindex , 'AAeventindex':AAeventindex , 'AAeventypes': AAeventypes  }
            else:
                eventdict[i] = { 'type': eventtypes , 'index' : eventindex , 'AAeventindex':[] , 'AAeventypes': [] }
            eventdict[i]['codon_pos'] = codonpos
            eventdict[i]['column'] = idx[0] + codonpos
        retdf = pd.DataFrame.from_dict(eventdict, orient = 'index' )
        retdfs.append(retdf)

    del t
    gc.collect()
    return retdfs

###compute spares matrics from results #######################################################################
@dask.delayed(pure=False)
def compute_matrices(  retdfs  ,  matsize , transitionsNT = 12 , transitionsAA = 380      ):
    count = 0
    AA_mutation = None
    nucleotide_mutation = None
    for resdf in retdfs:
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


    '''
    parser=argparse.ArgumentParser()
    parser.add_argument('--tree', help='tree corresponding to the alignment',type = str)
    parser.add_argument('--aln', help='alignment corresponding to the tree',type = str)
    parser.add_argument('--distributed', help='run in distributed mode' , type = str)
    parser.add_argument('--blastpath', help='blastpath',type = str)
    parser.add_argument('--refproteome', help='reference proteome' , type = str)
    parser.add_argument('--overwrite', help='overwrite HDF5 of aln' , type = str)
    parser.add_argument('--verbose', help='overwrite HDF5 of aln' , type = str)

    parser.add_argument('--reload', help='overwrite HDF5 of aln' , type = str)
    parser.add_argument('--njobs', help='overwrite HDF5 of aln' , type = str)


    args = vars(parser.parse_args(sys.argv[1:]))




    if args['tree']:
        treefile = args['tree']
    if args['aln']:
        alnfile = args['aln']
    if args['blastpath']:
        taxfilter = args['blastpath']
    if args['refproteome']:
        taxmask = args['refproteome']

    if args['reload']:
        reload_file = args['reload']

    if args['reload']:
        reload_file = args['reload']
    if args['reload']:
        taxmask = args['refproteome']
    if args['distributed']:
        if args['distributed'] == 'True':
            distributed_computation = True
        else:
            distributed_computation = False

    if args['distributed']:
        if args['distributed'] == 'True':
            distributed_computation = True
            njobs = 200
        else:
            distributed_computation = False
            njobs = 10

    if args['njobs']:
            njobs = int(args['njobs'])

    print(args)

    '''

    print('timestamp',ts)
    tag = 'small_test'

    #fraction of genomes to remove if jackknifing
    #bootstraps = [ .001 , .01, .1 , .5  ,  .9  , .99    ]
    bootstraps = [ None , .25  ]
    bootstrap_replicates = [1 , 5]

    #defaults

    verbose = True
    
    distributed_computation = True
    njobs = 200 


    overwrite = False 
    overwrite_annot = False
    overwrite_index = False
    overwrite_tree = False


    #restart computation after failure
    #reload_file = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/apr_4_2022/mmsa_2022-04-04/2022-04-04_masked.fa_02022-06-13T22:03:53.971781small_test_coevmats.pkl'
    #startposition = 15708

    reload_file = None
    startposition = 0


    #number of replicates
    restart = None
    nucleotides_only = False

    blastpath = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/software/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin/'
    refproteodb = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/refproteome/covidrefproteome.fasta'
    #refproteodb = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/structs/covid_structs.fasta'


    alnfile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/apr_4_2022/mmsa_2022-04-04/2022-04-04_masked.fa'
    treefile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/apr_4_2022/GISAID-hCoV-19-phylogeny-2022-02-28/global.tree'
    
    #treefile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/feb_2021/GISAID-hCoV-19-phylogeny-2021-02-21/global.tree'

    #alnfile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/dec_7/2021-12-07_masked.fa'
    #treefile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/dec_7/global.tree'

    #alnfile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/mmsa_2021-11-02/2021-11-02_masked.fa'
    #treefile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/mmsa_2021-11-02/global.tree'

    #alnfile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/msa_0730/msa_0730.fasta'
    #treefile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/msa_0730/global.tree'

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
    def clipID(ID):
        if '|' in ID:
            ID =  ID.split('|')[1]
        if '_' in ID:
            ID = ID.replace( '_' , ' ')
        return ID


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
        IDs = {i:clipID(rec.id) for i,rec in enumerate(msa)}
        IDindex = dict(zip( IDs.values() , IDs.keys() ) )
        print( [(t,IDindex[t]) for t in list(IDindex.keys())[0:10]] )
        with open( alnfile + '_IDs.pkl' , 'wb') as idxout:
            idxout.write(pickle.dumps(IDindex))
    
    ###########################################init tree
    if os.path.exists( alnfile +'preptree.pkl') and overwrite_tree==False:
        print('reloading tree')
        with open( alnfile +'preptree.pkl' , 'rb') as treein:
            tree = pickle.loads(treein.read()) 
        print('done')
    else:
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
        print('found leaves:', len(found))
        print('out of n nodes:', len(tree.leaf_nodes()))
        print('done')

        print('saving tree')
        with open( alnfile +'preptree.pkl' ,  'wb') as treeout:
            treeout.write( pickle.dumps(tree) )

    seq = IDs[100]
    with h5py.File(alnfile +'.h5', 'r') as hf:     
        align_array = hf['MSA2array']
        print(align_array.shape)
        aln_row = align_array[:,IDindex[seq]]
        nongap_cols = [ i for i,c in enumerate(list(aln_row)) if c != b'-' ]        
        print(aln_row[0:100])
        print('nongap:' , nongap_cols[0:100] , '...')    

    #make map from selected ref geno positions to aln columns
    column_map = { i:nongap_cols[i-1] for i in range(1,len(nongap_cols)+1)}
    ##annotate sequence with ref proteome
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
        
        positions = []
        for i,row in annotation.iterrows():
            positions = positions + [ i for i in range(row.qstart,row.qend, 3 ) ]

    else:
        #just seperate sequence into dummy codons
        #indexing starts at 1 for blast
        print('using dummy annot')
        dummy_annot = {'dummy_gene': { 'qstart':1 , 'qend':align_array.shape[0]-1 , 'evalue':0  }}
        annotation = pd.DataFrame.from_dict( dummy_annot , orient = 'index')
        positions = [ i for i in range(1,align_array.shape[0], 3 ) ]

    positions = sorted(list(set(positions)))
    print(annotation)
    print( 'calcilating sankof on positions', positions[0:500] , '...' )
    print('codons to calclulate:' , len(positions))
    print('flashing up a dask cluster')
    if distributed_computation == True:


        NCORE = 1
        print('deploying cluster')
        cluster = SLURMCluster(
            walltime='1:00:00',
            n_workers = njobs,
            cores=NCORE,
            processes = NCORE,
            interface='ib0',
            memory="50GB" ,
            env_extra=[
            'source /work/FAC/FBM/DBC/cdessim2/default/dmoi/condaenvs/etc/profile.d/conda.sh',
            'conda activate ML2'
            ],
            scheduler_options={'interface': 'ens2f0' }
        )
        print(cluster.job_script())
        cluster.scale(jobs=10)
        cluster.adapt(minimum=njobs, maximum=1000)
        time.sleep(5)
        print(cluster)
        print(cluster.dashboard_link)
        client = Client(cluster , timeout='1000s' , set_as_default=True )
        
    else:
        NCORE = 50
        print('testing')
        cluster = LocalCluster(n_workers = njobs )
        #cluster.adapt(minimum = 50,  maximum=NCORE)
        print(cluster.dashboard_link)
        client = Client(cluster)
    print('done')
    
    #######start the sankof algo here #######################
    print('starting sankof')
    row_index = IDindex


    print('done')
    coordinates = []
    while len(client.scheduler_info()['workers']) < 5:
        time.sleep(5)
        cluster.scale(jobs=njobs)
        cluster.adapt(minimum=njobs, maximum=1000)
        print('waiting for workers')

    with h5py.File(alnfile +'.h5', 'r') as hf:
        align_array = hf['MSA2array'] 
        retmatsize = ( len(tree.nodes()) ,align_array.shape[0]  )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for bs , bootstrap in enumerate(bootstraps):
            this_bootstrap_replicates = bootstrap_replicates[bs]
            for k in range(this_bootstrap_replicates):
                if reload_file == None:
                    if bootstrap is None:
                        filename = alnfile +'_' + str(k) +ts+ tag  + '_coevmats.pkl' 
                    else:
                        filename = alnfile +'_' + str(k) +ts+ tag + str(bootstrap ) + '_BS_coevmats.pkl' 
                else:
                    filename = reload_file
                    with open( filename , 'rb') as coevout:
                        print('reloading' , filename)
                        matricesAA, matricesNT  = pickle.loads( coevout.read() )
                        print('loaded',(matricesAA, matricesNT ) )
                print('working on ', filename)
                fitch_inlist = []
                positions_batch = []

               

                with open( filename+'.logfile.txt' , 'w') as logout:
                    logout.write('cluster link')
                    logout.write(cluster.dashboard_link)
                    logout.flush()
                    computebatch = 1
                    codonbatch = 5
                    #indexing starts at 1 for blast
                    #####switch to sending the coordinates and masking for the matrix
                    results = []
                    workbatch = []

                    for i,codon in enumerate(positions):
                        if ( reload_file is not None and len(positions) - i < startposition ) or reload_file is None:
                            if nucleotides_only == False:
                                pos = [column_map[codon], column_map[codon+1] , column_map[codon+2]]
                            else:
                                pos = [codon + i for i in range(3) if codon +i < align_array.shape[0]]
                            positions_batch.append(pos)
                        future = client.submit( calculate_small_parsimony , alnfile +'preptree.pkl' ,  positions_batch  , bootstrap , codon , alnfile = alnfile)
                        workbatch.append(future)
                        #loop here to compute results if q is too long
                        if len ( workbatch )> computebatch:
                            print('computing')

                            ac = as_completed(workbatch)     

                            fitch_inlist =  [f.result() for f in ac]

                            print(fitch_inlist)

                            delayed_mats = [ compute_matrices(df, retmatsize) for df in fitch_inlist ]
                            delayed_mats = dask.compute( * fitch_inlist , retries = 100)
                            AAbag = dask.bag.from_sequence([ m[1] for m in delayed_mats ]) 
                            NTbag = dask.bag.from_sequence([ m[0] for m in delayed_mats ])
                            if count ==0 :
                                matricesAA = dask.compute(AAbag.sum())[0]
                                matricesNT = dask.compute(NTbag.sum())[0]
                            else:
                                matricesAA += dask.compute(AAbag.sum())[0]
                                matricesNT += dask.compute(NTbag.sum())[0]
                            print(sparseND.argwhere(matricesAA))
                            
                            results = []
                            workbatch = []

                        if i % 50 == 0 and i > 0:
                            with open( filename , 'wb') as coevout:
                                print(filename)
                                print(count, 'intermediate saving',(matricesAA, matricesNT , i ) )
                                coevout.write(pickle.dumps((matricesAA, matricesNT , i )))
                                print('done')
                    #last batch
                    while len( workbatch ) > 0:                        
                            time.sleep(1)
                            for f in workbatch:
                                if f.done():
                                    results.append(f)
                                    workbatch.remove(f)
                    fitch_inlist =  [f.result() for f in results]
                    delayed_mats = [ compute_matrices(df, retmatsize) for df in fitch_inlist ]
                    AAbag = dask.bag.from_delayed([  m[1] for m in delayed_mats ]) 
                    NTbag = dask.bag.from_delayed([ m[0] for m in delayed_mats ])
                    matricesAA += dask.compute(AAbag.sum())[0]
                    matricesNT += dask.compute(NTbag.sum())[0]
                    print(matricesAA)
                    print(matricesNT)
                    print(sparseND.argwhere(matricesAA))
                    print('done')

                    with open( filename , 'wb') as coevout:
                        print(filename)
                        print('saving',(matricesAA, matricesNT ) )
                        coevout.write(pickle.dumps((matricesAA, matricesNT )))
                        print('done')
                    print('DONE!')