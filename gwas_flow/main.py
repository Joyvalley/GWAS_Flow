import sys
import pandas as pd
import numpy as np
from scipy.stats import f
import tensorflow as tf
from pandas_plink import read_plink
import h5py
from .herit import estimate_variance_components


def kinship(marker):
    ''' returns kinship matrix after vanRaden '''
    n_phe = marker.shape[0]
    n_mar = marker.shape[1]
    mafs = np.sum(marker, axis=0) / n_phe
    p_mat = np.repeat(mafs, n_phe).reshape(n_phe, n_mar, order="F")
    z_mat = marker - p_mat
    kin_vr = (np.matmul(z_mat, z_mat.T)) / (2 * np.sum(mafs * (1 - mafs)))
    return kin_vr


def load_and_prepare_data(x_file, y_file, k_file, m_phe, cof_file):
    ''' etl the data '''
    if k_file != 'not_prov':
        type_k = k_file.split(".")[-1]
    type_x = x_file.split(".")[-1]
    y_phe = pd.read_csv(y_file, engine='python').sort_values(
        ['accession_id']).groupby('accession_id').mean()
    y_phe = pd.DataFrame({'accession_id': y_phe.index,
                         'phenotype_value': y_phe[m_phe]})
    if type_x in ('hdf5', 'h5py'):
        snp = h5py.File(x_file, 'r')
        markers = np.asarray(snp['positions'])
        acc_x = np.asarray(snp['accessions'][:], dtype=np.int)
    elif type_x == 'csv':
        x_gen = pd.read_csv(x_file, index_col=0)
        markers = x_gen.columns.values
        acc_x = x_gen.index
        x_gen = np.asarray(x_gen, dtype=np.float64) / 2
    elif type_x.lower() == 'plink':
        my_prefix = x_file.split(".")[0]
        (bim, fam, bed) = read_plink(my_prefix)
        acc_x = np.array(fam[['fid']], dtype=np.str).flatten()
        markers = np.array(bim[['snp']]).flatten()
    else:
        sys.exit("Only hdf5, h5py, plink and csv files are supported")
    if k_file != 'not_prov':
        if type_k in ('hdf5', 'h5py'):
            k = h5py.File(k_file, 'r')
            acc_k = np.asarray(k['accessions'][:], dtype=np.int)
        elif type_k == 'csv':
            k = pd.read_csv(k_file, index_col=0)
            acc_k = k.index
            k = np.array(k, dtype=np.float64)

    acc_y = np.asarray(y_phe[['accession_id']]).flatten()
    acc_isec = [isec for isec in acc_x if isec in acc_y]
    if(len(acc_isec) == 0):
        print("WARNING: accessions in X do not overlap with accessions in Y")
        print("Accessions X:")
        print(acc_x)
        print("Accessions Y:")
        print(acc_y)

    idx_acc = list(map(lambda itt: itt in acc_isec, acc_x))
    idy_acc = list(map(lambda itt: itt in acc_isec, acc_y))
    if k_file != 'not_prov':
        idk_acc = list(map(lambda itt: itt in acc_isec, acc_k))
        if len(idk_acc) != len(acc_isec):
            print("WARNING: not all accessions are in the kinship matrix")
            print("Accessions X/Y:")
            print(acc_isec)
            print("Accessions K:")
            print(acc_k)
    if cof_file != 0:
        cof = pd.read_csv(cof_file, index_col=0)
        idc = cof.index
        cof = np.array(cof['cof'])
        acc_isec = [isec for isec in idc if isec in acc_y]
        #idc_acc = list(map(lambda x: x in acc_isec, idc))
        if not all(idx_acc):
            print('''
            accessions ids in the covariate file must be 
            identical to the ones in the phenotype file
            ''')
            sys.exit()
    else:
        cof = 0

    y_phe_ = np.asarray(y_phe.drop('accession_id', 1),
                        dtype=np.float64)[idy_acc, :]
    if type_x in ('hdf5', 'h5py'):
        x_gen = np.asarray(snp['snps'][:, np.where(idx_acc)[0]], np.float64).T
        x_gen = x_gen[np.argsort(acc_x[idx_acc]), :]
        if k_file != 'not_prov':
            k_1 = np.asarray(k['kinship'][:])[idk_acc, :]
            kin_vr = k_1[:, idk_acc]
            kin_vr = kin_vr[np.argsort(acc_x[idx_acc]), :]
            kin_vr = kin_vr[:, np.argsort(acc_x[idx_acc])]
        else:
            kin_vr = kinship(x_gen)
    elif type_x.lower() == 'plink':
        x_gen = np.asarray(bed.compute() / 2, dtype=np.float64)[:, idx_acc].T
        if k_file != 'not_prov':
            k_1 = k
            if 'kinship' in k:
                k_1 = k['kinship'][:]
            k_1 = np.asarray(k_1)[idk_acc, :]
            kin_vr = k_1[:, idk_acc]
            kin_vr = kin_vr[np.argsort(acc_x[idx_acc]), :]
            kin_vr = kin_vr[:, np.argsort(acc_x[idx_acc])]
        else:
            kin_vr = kinship(x_gen)
    else:
        x_gen = x_gen[idx_acc, :]
        if k_file != 'not_prov':
            k_1 = k[idk_acc, :]
            kin_vr = k_1[:, idk_acc]
        else:
            kin_vr = kinship(x_gen)

    print("data has been imported")
    return x_gen, kin_vr, y_phe_, markers, cof


def mac_filter(mac_min, x_gen, markers):
    ''' filter for minor allele frequencies'''
    ac1 = np.sum(x_gen, axis=0)
    ac0 = x_gen.shape[0] - ac1
    macs = np.minimum(ac1, ac0)
    markers_used = markers[macs >= mac_min]
    x_gen = x_gen[:, macs >= mac_min]
    return markers_used, x_gen, macs

# calculate betas and se of betas


def stderr_func(itt, marker, y_t2d, int_t):
    '''returns standard errors for induvidual markers '''
    n_phe = len(int_t)
    x = tf.stack(
        (int_t, tf.squeeze(
            tf.matmul(
                marker.T, tf.reshape(
                    itt, (n_phe, -1))))), axis=1)
    coeff = tf.matmul(
        tf.matmul(
            tf.linalg.inv(
                tf.matmul(
                    tf.transpose(x, perm=None),
                    x)),
            tf.transpose(x, perm=None)),
        y_t2d)
    sum_sq_e = tf.reduce_sum(tf.math.square(tf.math.subtract(y_t2d, tf.math.add(
        tf.math.multiply(x[:, 1], coeff[0, 0]), tf.math.multiply(x[:, 1], coeff[1, 0])))))
    stand_err = tf.math.sqrt(sum_sq_e / (n_phe - (1 + 2)))
    stdr_glob = tf.sqrt(
        tf.linalg.diag_part(
            tf.math.multiply(
                stand_err,
                tf.linalg.inv(
                    tf.matmul(
                        tf.transpose(x, perm=None),
                        x)))))[1]
    return tf.stack((coeff[1, 0], stdr_glob))

# calculate residual sum squares


def rss(itt, marker, y_t2d, int_t):
    ''' calculates the residual sum of squares '''
    x_t = tf.reduce_sum(tf.math.multiply(marker.T, itt), axis=1)
    lm_res = tf.linalg.lstsq(
        tf.transpose(
            tf.stack(
                (int_t, x_t), axis=0), perm=None), y_t2d, l2_regularizer=0.0)
    lm_x = tf.concat((tf.squeeze(lm_res), x_t), axis=0)
    return tf.reduce_sum(tf.math.square(tf.math.subtract(tf.squeeze(y_t2d), tf.math.add(
        tf.math.multiply(lm_x[1], lm_x[2:]), tf.multiply(lm_x[0], int_t)))))

# calculate residual sum squares with co-variates


def rss_cof(itt, marker, y_t2d, int_t, cof_t):
    ''' calculates the residual sum of squares when cof is included '''
    x_t = tf.reduce_sum(tf.math.multiply(marker.T, itt), axis=1)
    lm_res = tf.linalg.lstsq(
        tf.transpose(
            tf.stack(
                (int_t, x_t, cof_t), axis=0), perm=None), y_t2d, l2_regularizer=0.0)
    return tf.math.reduce_sum(tf.math.square(
        y_t2d - (lm_res[1] * x_t + lm_res[0] * int_t + lm_res[2] * cof_t)))


def get_k_stand(kin_vr):
    ''' obtains the standardized kinship matrix'''
    n_phe = kin_vr.shape[0]
    return (n_phe - 1) / np.sum((np.identity(n_phe) -
                                np.ones((n_phe, n_phe)) / n_phe)
                                * kin_vr) * kin_vr


def get_herit(y_phe, k_stand):
    ''' calculates the heritabilty'''
    return estimate_variance_components(y_phe, k_stand, verbose=False)


def transform_kinship(v_g, k_stand, v_e):
    ''' transform the kinship matrix with cholesky transformation '''
    n_phe = k_stand.shape[0]
    return np.transpose(
        np.linalg.inv(
            np.linalg.cholesky(
                v_g *
                k_stand +
                v_e *
                np.identity(n_phe)))).astype(
        np.float64)


def transform_y(marker, y_phe):
    ''' transform phenotypes '''
    return np.sum(np.multiply(np.transpose(marker), y_phe), axis=1).astype(np.float64)


def transform_int(marker):
    ''' transform the intercept'''
    n_phe = marker.shape[0]
    return np.sum(
        np.multiply(
            np.transpose(marker),
            np.ones(n_phe)),
        axis=1).astype(
            np.float64)


def emmax(int_t, y_trans):
    ''' run emmax according to Kang et al 2010'''
    n_phe = len(int_t)
    return (np.linalg.lstsq(np.reshape(int_t, (n_phe, -1)),
                            np.reshape(y_trans, (n_phe, -1)), rcond=None)[1]).astype(np.float64)


def transform_cof(marker, cof):
    ''' transform the coefficients '''
    return np.sum(np.multiply(np.transpose(marker), cof), axis=1).astype(np.float64)


def get_output(f_1, x_sub, stdr_glob):
    ''' get the F1 values'''
    return tf.concat([tf.reshape(f_1, (x_sub.shape[1], -1)), stdr_glob], axis=1)


def get_stderr(marker, y_t2d, int_t, x_sub):
    ''' build tensor loping of all markers to obtain all standerros'''
    return tf.map_fn(lambda mar: stderr_func(mar, marker, y_t2d, int_t), x_sub.T)


def get_f1(rss_env, r1_full, n_phe):
    '''calculate the f1 scors for all markers'''
    return tf.divide(
        tf.subtract(
            rss_env, r1_full), tf.divide(
            r1_full, (n_phe - 3)))


def get_pval(f_dist, n_phe):
    '''get p values from f1 scores'''
    return f.logsf(f_dist, 1, n_phe - 3)


def get_r1_full(marker, y_t2d, int_t, x_sub):
    ''' build tensor for full model '''
    return tf.map_fn(lambda mar: rss(mar, marker, y_t2d, int_t), x_sub.T)


def gwas(x_gen, kin_vr, y_phe, batch_size, cof):
    ''' get gwas results, calls all the subfunctions '''
    y_phe = y_phe.flatten()
    n_marker = x_gen.shape[1]
    n_phe = len(y_phe)
    # REML
    k_stand = get_k_stand(kin_vr)
    v_g, v_e = get_herit(y_phe, k_stand)
    print(" Pseudo-heritability is ", v_g / (v_e + v_g))
    print(" Performing GWAS on ", n_phe,
          " phenotypes and ", n_marker, "markers")
    # Transform kinship-matrix, phenotypes and estimate intercpt
   #  Xo = np.ones(K.shape[0]).flatten()
    marker = transform_kinship(v_g, k_stand, v_e)
    y_trans = transform_y(marker, y_phe)
    int_t = transform_int(marker)
    # transform  co-factor
    if isinstance(cof, int) == False:
        cof_t = transform_cof(marker, cof)
    # EMMAX Scan
    rss_env = emmax(int_t, y_trans)
    # loop over th batches
    for i in range(int(np.ceil(n_marker / batch_size))):
        if n_marker < batch_size:
            x_sub = x_gen
        else:
            lower_limit = batch_size * i
            upper_limit = batch_size * i + batch_size
            if upper_limit <= n_marker:
                x_sub = x_gen[:, lower_limit:upper_limit]
                print(
                    "Working on markers ",
                    lower_limit,
                    " to ",
                    upper_limit,
                    " of ",
                    n_marker)
            else:
                x_sub = x_gen[:, lower_limit:]
                print(
                    "Working on markers ",
                    lower_limit,
                    " to ",
                    n_marker,
                    " of ",
                    n_marker)
        y_t2d = tf.cast(tf.reshape(y_trans, (n_phe, -1)), dtype=tf.float64)
     
        stdr_glob = get_stderr(marker, y_t2d, int_t, x_sub)
        if isinstance(cof, int) == False:
            r1_full = tf.map_fn(
                lambda mar: rss_cof(
                    mar, marker, y_t2d, int_t, cof_t), x_sub.T)
        else:
            r1_full = get_r1_full(marker, y_t2d, int_t, x_sub)
        f_1 = get_f1(rss_env, r1_full, n_phe)
        if i == 0:
            output = get_output(f_1, x_sub, stdr_glob)
        else:
            tmp = get_output(f_1, x_sub, stdr_glob)
            output = np.append(output, tmp, axis=0)
        f_dist = output[:, 0]
    pval = np.exp(get_pval(f_dist, n_phe))
    pval = tf.expand_dims(pval,-1)
    output = tf.concat([output,pval],1)
    return output
