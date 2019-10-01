import pandas
import numpy as np
import pandas as pd
import h5py
import tensorflow as tf

X = h5py.File("./gwas_sample_data/AT_geno.hdf5","r")
SNPs = np.asarray(X['snps'],dtype='int8').T



X = tf.convert_to_tensor(SNPs,dtype=tf.float64)
X = X - tf.math.reduce_mean(X,axis=0)/tf.math.reduce_std(X,axis=0)
K = 1/X.shape[1] * tf.tensordot(X,tf.transpose(X),axes=1)


#K = np.zeros(shape=(SNPs.shape[0],SNPs.shape[0]))
X = (SNPs-np.mean(SNPs,axis=0))/np.std(SNPs,axis=0)
kernel = 1.0/X.shape[1] * np.dot(X,X.T)


#for i in range(SNPs.shape[0]):
#    for j in list(range(i,SNPs.shape[0])):
#        K[i,j] = sum(SNPs[i,:] == SNPs[i,:]) / SNPs.shape[1]


print("bla")
