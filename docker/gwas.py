import os
import sys
import time
import numpy as np
import pandas as pd
import main
import h5py

# set defaults 
mac_min = 1
batch_size =  500000 
out_file = "results.csv"
m = 'phenotype_value'
perm = 1
mac_min= 6

X_file = 'gwas_sample_data/AT_geno.hdf5'
Y_file = 'gwas_sample_data/phenotype.csv'
K_file = 'gwas_sample_data/kinship_ibs_binary_mac5.h5py'



for i in range (1,len(sys.argv),2):
    if sys.argv[i] == "-x" or sys.argv[i] == "--genotype":
        X_file = sys.argv[i+1]
    elif sys.argv[i] == "-y" or sys.argv[i] == "--phenotype":
        Y_file = sys.argv[i+1]
    elif sys.argv[i] == "-k" or sys.argv[i] == "--kinship":
        K_file = sys.argv[i+1]
    elif sys.argv[i] == "-m":
        m = sys.argv[i+1]
    elif sys.argv[i] == "-a" or sys.argv[i] == "--mac_min":
        mac_min = int(sys.argv[i+1])
    elif sys.argv[i] == "-bs" or sys.argv[i] == "--batch-size":
        batch_size = int(sys.argv[i+1])
    elif sys.argv[i] == "-p" or sys.argv[i] == "--perm":
        perm  = int(sys.argv[i+1])
    elif sys.argv[i] == "-o" or sys.argv[i] == "--out":
        out_file = sys.argv[i+1]
    elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
        print("-x , --genotype :file containing marker information in csv or hdf5 format of size")
        print("-y , --phenotype: file container phenotype information in csv format"  )
        print("-k , --kinship : file containing kinship matrix of size k X k in csv or hdf5 format")
        print("-m : name of columnn containing the phenotype : default m = phenotype_value")
        print("-a , --mac_min : integer specifying the minimum minor allele count necessary for a marker to be included. Default a = 1" )
        print("-bs, --batch-size : integer specifying the number of markers processed at once. Default -bs 500000" )
        print("-p , --perm : single integer specifying the number of permutations. Default 1 == no perm ")
        print("-o , --out : name of output file. Default -o results.csv  ")
        print("-h , --help : prints help and command line options")
        quit()
    else:
        print('unknown option ' + str(sys.argv[i]))
        quit()



print("parsed commandline args")

start = time.time()

X,K,Y_,markers = main.load_and_prepare_data(X_file,Y_file,K_file,m)


## MAF filterin
markers_used , X , macs = main.mac_filter(mac_min,X,markers)

## prepare
print("Begin performing GWAS on ", Y_file)

if perm == 1:
    output = main.gwas(X,K,Y_,batch_size)   
    if( X_file.split(".")[-1] == 'csv'):
        chr_pos = np.array(list(map(lambda x : x.split("- "),markers_used)))
    else: 
        chr_reg = h5py.File(X_file,'r')['positions'].attrs['chr_regions']
        mk_index= np.array(range(len(markers)),dtype=int)[macs >= mac_min]
        chr_pos = np.array([list(map(lambda x: sum(x > chr_reg[:,1]) + 1, mk_index)), markers_used]).T
        my_time = np.repeat((time.time()-start),len(chr_pos))
    pd.DataFrame({
        'chr' : chr_pos[:,0] ,
        'pos' : chr_pos[:,1] , 
        'pval': output[:,0] ,
        'mac' : np.array(macs[macs >= mac_min],dtype=np.int) ,
        'eff_size': output[:,1] ,
        'SE' : output[:,2]}).to_csv(out_file,index=False)
elif perm > 1:
    min_pval = []
    perm_seeds = []
    my_time = []
    for i in range(perm):
        start_perm = time.time()
        print("Running permutation ", i+1, " of ",perm)
        my_seed  = np.asscalar(np.random.randint(9999,size=1))
        perm_seeds.append(my_seed)
        np.random.seed(my_seed)
        Y_perm = np.random.permutation(Y_)
        output = main.gwas(X,K,Y_perm,batch_size)
        min_pval.append(np.min(output[:,0]))
        print("Elapsed time for permuatation",i+1 ," with p_min", min_pval[i]," is",": ", round(time.time() - start_perm,2))
        my_time.append(time.time()-start_perm)
    pd.DataFrame({
        'time': my_time ,
        'seed': perm_seeds ,
        'min_p': min_pval }).to_csv(out_file,index=False)

print("done")
 
end = time.time()
eltime = np.round(end -start,2)

if eltime <= 59:
    print("Total time elapsed",  eltime, "seconds")
elif eltime > 59 and eltime <= 3600:
    print("Total time elapsed",  np.round(eltime / 60,2) , "minutes")
elif eltime > 3600 :
    print("Total time elapsed",  np.round(eltime / 60 / 60,2), "hours")
