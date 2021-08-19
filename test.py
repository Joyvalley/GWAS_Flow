''' Test functions for gwas '''
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import unittest
import pickle
import numpy as np
import main


class TestGwas(unittest.TestCase):
    ''' class containing tests for gwas'''
    def test_kinship(self):
        ''' testing the kinship function '''
        with open('test_data/K_test', 'rb') as k_file:
            kin_mat = pickle.load(k_file)
        with open('test_data/M_test', 'rb') as m_file:
            markers = pickle.load(m_file)
        self.assertIsNone(np.testing.assert_array_equal(main.kinship(markers), kin_mat))

    def test_gwas(self):
        ''' tests for the gwas function '''
        with open('test_data/X_test', 'rb') as g_file:
            x_gen = pickle.load(g_file)
        with open('test_data/K_gwas_test', 'rb') as k_file:
            kin_mat = pickle.load(k_file)
        with open('test_data/Y_test', 'rb') as phe_file:
            y_pheno = pickle.load(phe_file)
        with open('test_data/Out_test', 'rb') as output:
            output = pickle.load(output)
        with open('test_data/cof_test', 'rb') as coutt:
            cof = pickle.load(coutt)
        with open('test_data/cof_output', 'rb') as coutput:
            cof_output = pickle.load(coutput)
        batch_size = 500000
        print(output)
        self.assertIsNone(
            np.testing.assert_array_almost_equal(
                main.gwas(
                    x_gen,
                    kin_mat,
                    y_pheno,
                    batch_size,
                    cof=0),
                output))
        print(output)
        self.assertIsNone(
            np.testing.assert_array_almost_equal(
                main.gwas(
                    x_gen,
                    kin_mat,
                    y_pheno,
                    batch_size,
                    cof),
                cof_output))
