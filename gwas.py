''' main script for gwas '''
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
from gwas_flow import main
import h5py as h5
import pandas as pd
from pandas_plink import read_plink
import numpy as np
import time
import sys

# set defaults
BATCH_SIZE = 500000
OUT_FILE = "results.csv"
M_PHE = 'phenotype_value'
PERM = 1
MAC_MIN = 1
COF_FILE = 0
COF = "nan"
PLOT = False
K_FILE = 'not_prov'
OUT_PERM = False


for i in range(1, len(sys.argv), 2):
    if sys.argv[i] == "-x" or sys.argv[i] == "--genotype":
        X_FILE = sys.argv[i + 1]
    elif sys.argv[i] == "--cof":
        COF_FILE = sys.argv[i + 1]
    elif sys.argv[i] == "-y" or sys.argv[i] == "--phenotype":
        Y_FILE = sys.argv[i + 1]
    elif sys.argv[i] == "-k" or sys.argv[i] == "--kinship":
        K_FILE = sys.argv[i + 1]
    elif sys.argv[i] == "-m":
        M_PHE = sys.argv[i + 1]
    elif sys.argv[i] == "-a" or sys.argv[i] == "--mac_min":
        MAC_MIN = int(sys.argv[i + 1])
    elif sys.argv[i] == "-bs" or sys.argv[i] == "--batch-size":
        BATCH_SIZE = int(sys.argv[i + 1])
    elif sys.argv[i] == "-p" or sys.argv[i] == "--perm":
        PERM = int(sys.argv[i + 1])
    elif sys.argv[i] == "-o" or sys.argv[i] == "--out":
        OUT_FILE = sys.argv[i + 1]
    elif sys.argv[i] == "--plot":
        PLOT = bool(sys.argv[i + 1])
    elif sys.argv[i] == "--out_perm":
        OUT_PERM = bool(sys.argv[i + 1])
    elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
        print("-x , --genotype :file containing marker information in csv or hdf5 format of size")
        print("-y , --phenotype: file container phenotype information in csv format")
        print("-k , --kinship : file containing kinship matrix of size k X k in csv or hdf5 format")
        print("-m : name of columnn containing the phenotype : default m = phenotype_value")
        print('''
        -a , --mac_min : integer specifying the minimum minor allele count necessary 
        for a marker to be included. Default a = 1
        ''')
        print('''
        -bs, --batch-size : integer specifying the number of markers processed at once.
        Default -bs 500000
        ''')
        print('''
        -p , --perm : single integer specifying the number of permutations. 
        Default 1 == no perm 
        ''')
        print('''
        --out_perm : output the results of the individual permuations.
        Default False
        ''')
        print("-o , --out : name of output file. Default -o results.csv  ")
        print("-h , --help : prints help and command line options")
        print("--plot: creates manhattan plot")
        sys.exit()
    else:
        print('unknown option ' + str(sys.argv[i]))
        sys.exit()


print("parsed commandline arguments")

start = time.time()

X, K, Y_, markers, COF = main.load_and_prepare_data(
    X_FILE, Y_FILE, K_FILE, M_PHE, COF_FILE)

# MAF filterin
markers_used, X, macs = main.mac_filter(MAC_MIN, X, markers)

# prepare
print("Begin performing GWAS on ", Y_FILE)

output = main.gwas(X, K, Y_, BATCH_SIZE, COF)
if X_FILE.split(".")[-1] == 'csv':
    CHR_POS = np.array(list(map(lambda x: x.split("- "), markers_used)))
elif X_FILE.split(".")[-1].lower() == 'plink':
    my_prefix = X_FILE.split(".")[0]
    (bim, fam, bed) = read_plink(my_prefix)
    bim.set_index("snp", inplace=True)
    my_chr = bim.loc[markers_used, "chrom"]
    my_pos = bim.loc[markers_used, "pos"]
    CHR_POS = np.vstack((my_chr, my_pos)).T
else:
    chr_reg = h5.File(X_FILE, 'r')['positions'].attrs['chr_regions']
    mk_index = np.array(range(len(markers)), dtype=int)[macs >= MAC_MIN]
    CHR_POS = np.array(
        [list(map(lambda x: sum(x > chr_reg[:, 1]) + 1, mk_index)), markers_used]).T
    my_time = np.repeat((time.time() - start), len(CHR_POS))
res = pd.DataFrame({
    'chr': CHR_POS[:, 0],
    'pos': CHR_POS[:, 1],
    'pval': output[:, 0],
    'mac': np.array(macs[macs >= MAC_MIN], dtype=np.int),
    'eff_size': output[:, 1],
    'SE': output[:, 2]})
res.to_csv(OUT_FILE, index=False)
if PERM > 1:
    min_pval = []
    perm_seeds = []
    my_time = []
    for i in range(PERM):
        perm_out = 'perm_' + OUT_FILE
        start_perm = time.time()
        print("Running permutation ", i + 1, " of ", PERM)
        my_seed = np.asscalar(np.random.randint(9999, size=1))
        perm_seeds.append(my_seed)
        np.random.seed(my_seed)
        Y_perm = np.random.permutation(Y_)
        output = main.gwas(X, K, Y_perm, BATCH_SIZE, COF)
        min_pval.append(np.min(output[:, 0]))
        if OUT_PERM:
            res = pd.DataFrame({
                'chr': CHR_POS[:, 0],
                'pos': CHR_POS[:, 1],
                'pval': output[:, 0],
                'mac': np.array(macs[macs >= MAC_MIN], dtype=np.int),
                'eff_size': output[:, 1],
                'SE': output[:, 2]})
            res.to_csv(OUT_FILE.replace(
                ".csv", "_" + str(i + 1) + ".csv"), index=False)

        print(
            "Elapsed time for permuatation",
            i + 1, " with p_min", min_pval[i],
            " is", ": ", round(
                time.time() -
                start_perm,
                2))
        my_time.append(time.time() - start_perm)
    res_perm = pd.DataFrame({
        'time': my_time,
        'seed': perm_seeds,
        'min_p': min_pval}).sort_values('min_p')
    res_perm.to_csv(perm_out, index=False)


if PLOT:
    import src.plot as plot
    plot.manhattan(OUT_FILE, PERM)


print("Finished performing GWAS")


end = time.time()
eltime = np.round(end - start, 2)

if eltime <= 59:
    print("Total time elapsed", eltime, "seconds")
elif 59 < eltime <= 3600:
    print("Total time elapsed", np.round(eltime / 60, 2), "minutes")
elif eltime > 3600:
    print("Total time elapsed", np.round(eltime / 60 / 60, 2), "hours")
end = time.time()
