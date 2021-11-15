''' custom function from limix package that calculates the variance components estiamtes '''
import numpy as np
from numpy_sugar.linalg import economic_qs
from numpy import pi, var
from numpy_sugar import is_all_finite
from glimix_core.lmm import LMM



def estimate_variance_components(y_phe, kin, marker_mat=None, verbose=True):
    if not is_all_finite(y_phe):
        raise ValueError("Outcome must have finite values only.")
    if not is_all_finite(kin):
        raise ValueError("Outcome must have finite values only.")
    marker_mat = np.full ((y_phe.shape[0],1),1.) 
    if kin is not None:
           # K = K / diag(K).mean()
            q_s = economic_qs(kin)
    else:
            q_s = None
       
    method = LMM(y_phe, marker_mat, q_s, restricted=True)
    method.fit(verbose=verbose)
        
    v_g = method.scale * (1 - method.delta)
    v_e = method.scale * method.delta
    
    return v_g, v_e
