import main 
import unittest
import pickle 
import numpy as np


class TestGwas(unittest.TestCase):
    def test_K(self):          
        with open('test_data/K_test','rb') as f: K = pickle.load(f)
        with open('test_data/M_test','rb') as g: M = pickle.load(g)
        self.assertIsNone(np.testing.assert_array_equal(main.kinship(M),K))
        
    def test_gwas(self):
        with open('test_data/X_test','rb') as f: X = pickle.load(f)
        with open('test_data/K_gwas_test','rb') as f: K = pickle.load(f)
        with open('test_data/Y_test','rb') as f: Y = pickle.load(f)
        with open('test_data/Out_test','rb') as f: output = pickle.load(f)
        batch_size=500000        
        cof = 0
        self.assertIsNone(np.testing.assert_array_equal(main.gwas(X,K,Y,batch_size,cof),output))

