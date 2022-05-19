import os
import h5py
import numpy as np

overwrite = True

#alnfile = '/scratch/dmoi/datasets/covid_data/dec_7/2021-12-07_masked.fa'
#alnfile = '../validation_data/16s/16s_salaminWstruct_aln.fasta'
alnfile = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/covid_data/apr_4_2022/mmsa_2022-04-04/2022-04-04_masked.fa'

print('aln2hdf5 start')

if not (os.path.exists(alnfile +'.h5') and overwrite == False):
    def ret_msa(alnfile):
        with open( alnfile ,'r') as msa:
            ret = []
            s = None
            count = 0
            for l in msa:
                if '>' not in l:
                    if s:
                        s+=l.strip()
                    else:
                        s = l.strip()
                else:
                    if s:
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
        print(dset)
    print('done')