import pandas
import numpy as np
import pandas as pd
import h5py

X = h5py.File("./gwas_sample_data/AT_geno.hdf5","r")
SNPs = np.asarray(X['snps'],dtype='int8').T

SNPs = SNPs[:,5000:5220]

K = np.zeros(shape=(SNPs.shape[0],SNPs.shape[0]))

for i in range(SNPs.shape[0]):
    for j in list(range(i,SNPs.shape[0])):
        K[i,j] = sum(SNPs[i,:] == SNPs[i,:]) / SNPs.shape[1]


print("bla")
