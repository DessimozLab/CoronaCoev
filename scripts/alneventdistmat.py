import pickle
import numpy as np

with open( 'blurmat_codons.pkl' , 'rb') as blurout:
    blurmat = pickle.loads(blurout.read())

blurmat= blurmat.todense()

import h5py
import multiprocessing as mp
#calculate a distmat between cols w events
def distfun(stuff):
    i,j,v1,v2 = stuff
    #distance function between two continuous vectors
    #euclidean dist
    return (i,j,np.linalg.norm( v1-v2 ))

def process(q,retq, iolock):
    #calculate compression distances
    with iolock:
        print('init worker')
    from time import sleep
    while True:
        stuff = q.get()
        if stuff is None:
            break
        retq.put(distfun(stuff))
    print('done')

def mat_creator(retq,matsize,iolock):
    with iolock:
        print('init matcreator')
    #collect distances and create final mat
    calculations = (matsize**2 - matsize) / 2
    distmat = np.zeros((matsize,matsize))
    count = 0
    with h5py.File('./gisaid/alnEventdistmat.h5', 'w') as hf:
        try:
            hf.create_dataset("alnEventdists",  data=distmat)
        except:
            pass
        while True:
            r = retq.get()
            count+=1
            if r is None:
                break
            row,col,cdist = r
            hf['alnEventdists'][row,col] = cdist
            if count% 10000 == 0 :
                with iolock:
                    print(count/calculations)
                    print((row,col))
                hf.flush()

    print('done saver')

startk = 0
startl = 0

if __name__ == '__main__':

    NCORE = 40
    q = mp.Queue(maxsize=NCORE*5000)
    retq = mp.Queue(maxsize=NCORE*5000)
    iolock = mp.Lock()
    pool = mp.Pool(NCORE, initializer=process, initargs=(q,retq, iolock))
    p = mp.Process(target=mat_creator, args=(retq,blurmat.shape[1], iolock))
    p.start()

    for i in range(blurmat.shape[1]):
        for j in range(blurmat.shape[1]):
            if i < j :
                s1 = blurmat[:,i].ravel()
                s2 = blurmat[:,j].ravel()
                q.put( (i,j,s1,s2) )
                if i % 500 == 0 and j % 500 == 0 :
                    print(i)
                    print(j)

    for _ in range(NCORE):  # tell workers we're done
        q.put(None)
    retq.put(None)
    pool.close()
    pool.join()
    p.join()
