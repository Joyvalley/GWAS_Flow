import pandas as pd 
import numpy as np
import sys
from scipy.stats import f
import tensorflow as tf
import limix
import herit
import h5py
import multiprocessing as mlt
from pandas_plink import read_plink
import pickle

def kinship(M):
        mafs = np.sum(M,axis=0)/M.shape[0]
        P = np.repeat(mafs,M.shape[0]).reshape(M.shape[0],M.shape[1],order="F")
        Z = M - P
        K = (np.matmul(Z , Z.T))  / (2 * np.sum(mafs*(1-mafs)))
        return K

def load_and_prepare_data(X_file,Y_file,K_file,m,cof_file):  
    if K_file != 'not_prov':
        type_K = K_file.split(".")[-1]
    type_X = X_file.split(".")[-1]

        ## load and preprocess genotype matrix 
    Y = pd.read_csv(Y_file,engine='python').sort_values(['accession_id']).groupby('accession_id').mean()
    Y = pd.DataFrame({'accession_id' :  Y.index, 'phenotype_value' : Y[m]})
    if type_X == 'hdf5' or type_X == 'h5py'  :
        SNP = h5py.File(X_file,'r')
        markers= np.asarray(SNP['positions'])
        acc_X =  np.asarray(SNP['accessions'][:],dtype=np.int)
    elif type_X == 'csv' :
        X = pd.read_csv(X_file,index_col=0)
        markers = X.columns.values
        acc_X = X.index
        X = np.asarray(X,dtype=np.float32)/2
    elif type_X.lower() == 'plink':
        my_prefix = X_file.split(".")[0]
        (bim,fam,bed) = read_plink(my_prefix)
        acc_X = np.array(fam[['fid']],dtype=np.int).flatten()
        markers = np.array(bim[['snp']]).flatten()       
    else :
        sys.exit("Only hdf5, h5py, plink and csv files are supported")
    if K_file != 'not_prov':  
        if type_K == 'hdf5' or type_K == 'h5py':
            k = h5py.File(K_file,'r')
            acc_K = np.asarray(k['accessions'][:],dtype=np.int)
        elif type_K == 'csv':
            k = pd.read_csv(K_file,index_col=0)
            acc_K = k.index
            k = np.array(k, dtype=np.float32)

    acc_Y =  np.asarray(Y[['accession_id']]).flatten()
    acc_isec = [isec for isec in acc_X if isec in acc_Y]
           
    idx_acc = list(map(lambda x: x in acc_isec, acc_X))
    idy_acc = list(map(lambda x: x in acc_isec, acc_Y))
    if K_file != 'not_prov':
        idk_acc = list(map(lambda x: x in acc_isec, acc_K))
    else:
        idk_acc = idx_acc
       
        
    cof = 0
    if cof_file != 0 :
        cof = pd.read_csv(cof_file,index_col=0)
        idc = cof.index
        cof= np.array(cof['cof'])
        acc_isec = [isec for isec in idc if isec in acc_Y]
        idc_acc = list(map(lambda x: x in acc_isec, idc))
        if not all(idx_acc):
            print("accessions ids in the covariate file must be identical to the ones in the phenotype file")
            quit()
    Y_ = np.asarray(Y.drop('accession_id',1),dtype=np.float32)[idy_acc,:]

    if type_X == 'hdf5' or type_X == 'h5py' :
        X = np.asarray(SNP['snps'][0:(len(SNP['snps'])+1),],dtype=np.float32)[:,idx_acc].T
        X = X[np.argsort(acc_X[idx_acc]),:]
        if K_file != 'not_prov':
            k1 = np.asarray(k['kinship'][:])[idk_acc,:]
            K  = k1[:,idk_acc]
            K = K[np.argsort(acc_X[idx_acc]),:]
            K = K[:,np.argsort(acc_X[idx_acc])]
        else: 
            K = kinship(X)          
    elif type_X.lower() == 'plink':
        X = np.asarray(bed.compute()/2,dtype=np.float32)[:,idx_acc].T
        if K_file != 'not_prov' : 
            k1 = np.asarray(k['kinship'][:])[idk_acc,:]
            K  = k1[:,idk_acc]
            K = K[np.argsort(acc_X[idx_acc]),:]
            K = K[:,np.argsort(acc_X[idx_acc])]
        else:
            K = kinship(X)
    else:
        X  = X[idx_acc,:]
        if K_file != 'not_prov':
            k1 = k[idk_acc,:]
            K  = k1[:,idk_acc]
        else:
            K = kinship(X)
       
    print("data has been imported")
    return X,K,Y_,markers,cof


def mac_filter(mac_min, X, markers):
    ac1 = np.sum(X,axis=0)
    ac0 = X.shape[0] - ac1
    macs = np.minimum(ac1,ac0)
    markers_used  = markers[macs >= mac_min]
    X = X[:,macs >= mac_min]
    return markers_used, X, macs


## calculate betas and se of betas 
def stderr(a,M,Y_t2d,int_t):
        n = len(int_t)
        x = tf.stack((int_t,tf.squeeze(tf.matmul(M.T,tf.reshape(a,(n,-1))))),axis=1)
        coeff = tf.matmul(tf.matmul(tf.linalg.inv(tf.matmul(tf.transpose(x),x)),tf.transpose(x)),Y_t2d)
        SSE = tf.reduce_sum(tf.math.square(tf.math.subtract(Y_t2d,tf.math.add(tf.math.multiply(x[:,1],coeff[0,0]),tf.math.multiply(x[:,1],coeff[1,0])))))
        SE = tf.math.sqrt(SSE/(471-(1+2)))
        StdERR = tf.sqrt(tf.linalg.diag_part(tf.math.multiply(SE , tf.linalg.inv(tf.matmul(tf.transpose(x),x)))))[1]
        return tf.stack((coeff[1,0],StdERR))

## calculate residual sum squares 
def rss(a,M,Y_t2d,int_t):
        x_t = tf.reduce_sum(tf.math.multiply(M.T,a),axis=1)
        lm_res = tf.linalg.lstsq(tf.transpose(tf.stack((int_t,x_t),axis=0)),Y_t2d)
        lm_x = tf.concat((tf.squeeze(lm_res),x_t),axis=0)
        return tf.reduce_sum(tf.math.square(tf.math.subtract(tf.squeeze(Y_t2d),tf.math.add(tf.math.multiply(lm_x[1],lm_x[2:]), tf.multiply(lm_x[0],int_t)))))

# calculate residual sum squares with co-variates
def rss_cof(a,M,Y_t2d,int_t,cof_t):
        x_t = tf.reduce_sum(tf.math.multiply(M.T,a),axis=1)
        lm_res = tf.linalg.lstsq(tf.transpose(tf.stack((int_t,x_t,cof_t),axis=0)),Y_t2d)
        return tf.math.reduce_sum(tf.math.square(Y_t2d - (lm_res[1] * x_t + lm_res[0] * int_t + lm_res[2] * cof_t)))
    
def get_K_stand(K):
    n = K.shape[0]
    return (n-1)/np.sum((np.identity(n) - np.ones((n,n))/n) * K) * K

def get_herit(Y,K_stand):
    return herit.estimate(Y,"normal",K_stand,verbose = False)

def transform_kinship(vg, K_stand, ve): 
    n = K_stand.shape[0]
    return np.transpose(np.linalg.inv(np.linalg.cholesky(vg * K_stand + ve  * np.identity(n)))).astype(np.float32)

def transform_Y(M,Y):
    return  np.sum(np.multiply(np.transpose(M),Y),axis=1).astype(np.float32)

def transform_int(M):
    n = M.shape[0]
    return np.sum(np.multiply(np.transpose(M),np.ones(n)),axis=1).astype(np.float32)

def emmax(int_t,Y_t):
    n = len(int_t)
    return (np.linalg.lstsq(np.reshape(int_t,(n,-1)) , np.reshape(Y_t,(n,-1)),rcond=None)[1]).astype(np.float32)

def transform_cof(M,cof):
    return np.sum(np.multiply(np.transpose(M),cof),axis=1).astype(np.float32)

def  getOutput(F_1,X_sub,StdERR):
    return tf.concat([tf.reshape(F_1,(X_sub.shape[1],-1)),StdERR],axis=1)

def getSTDERR(M,Y_t2d,int_t,X_sub):
    return tf.map_fn(lambda a : stderr(a,M,Y_t2d,int_t), X_sub.T)

def getF1(RSS_env,R1_full,n):
    return tf.divide(tf.subtract(RSS_env, R1_full),tf.divide(R1_full,(n-3)))

def getPval(F_dist,n):
    return 1 - f.cdf(F_dist,1,n-3)

def getR1Full(M,Y_t2d,int_t,X_sub):
     return tf.map_fn(lambda a: rss(a,M,Y_t2d,int_t), X_sub.T)

def gwas(X,K,Y,batch_size,cof):
    Y = Y.flatten()
    n_marker = X.shape[1]
    n = len(Y)
    ## REML   
    K_stand = get_K_stand(K)
    vg, delta, ve  = get_herit(Y,K_stand)
    print(" Pseudo-heritability is " , vg / (ve + vg + delta))
    print(" Performing GWAS on ", n , " phenotypes and ", n_marker ,"markers")
    ## Transform kinship-matrix, phenotypes and estimate intercpt
   #  Xo = np.ones(K.shape[0]).flatten()
    M = transform_kinship(vg, K_stand, ve)
    Y_t =   transform_Y(M,Y)
    int_t = transform_int(M)
    #transform  co-factor
    if isinstance(cof,int) == False:
        cof_t = transform_cof(M,cof)
    ## EMMAX Scan
    RSS_env = emmax(int_t,Y_t)
    ## loop over th batches 
    for i in range(int(np.ceil(n_marker/batch_size))):
        tf.compat.v1.reset_default_graph()
        if n_marker < batch_size:
            X_sub = X
        else:
            lower_limit = batch_size * i 
            upper_limit = batch_size * i + batch_size
            if upper_limit <= n_marker :
                X_sub = X[:,lower_limit:upper_limit]
                print("Working on markers ", lower_limit , " to ", upper_limit, " of ", n_marker )    
            else:
                X_sub = X[:,lower_limit:]
                print("Working on markers ", lower_limit , " to ", n_marker, " of ", n_marker )    
        config = tf.compat.v1.ConfigProto()
        n_cores = mlt.cpu_count()
        config.intra_op_parallelism_threads = n_cores
        config.inter_op_parallelism_threads = n_cores
        sess = tf.compat.v1.Session(config=config)                                             
        Y_t2d = tf.cast(tf.reshape(Y_t,(n,-1)),dtype=tf.float32)                     
      #  y_tensor =  tf.convert_to_tensor(Y_t,dtype = tf.float32)                                      
        StdERR = getSTDERR(M,Y_t2d,int_t,X_sub)              
        if isinstance(cof,int) == False :
            R1_full = tf.map_fn(lambda a: rss_cof(a,M,Y_t2d,int_t,cof_t), X_sub.T)
        else:
            R1_full = getR1Full(M,Y_t2d,int_t,X_sub)
        F_1 = getF1(RSS_env,R1_full,n)
        if i == 0 :
            output = sess.run(getOutput(F_1,X_sub,StdERR))
        else :
            tmp = sess.run(getOutput(F_1,X_sub,StdERR))
            output = np.append(output,tmp,axis=0)
        sess.close()
        F_dist = output[:,0]
    pval  = getPval(F_dist,n)
    output[:,0] = pval
    return output 


