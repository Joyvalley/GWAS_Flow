from numpy_sugar.linalg import economic_qs
from numpy import pi, var
from glimix_core.glmm import GLMMExpFam
from glimix_core.lmm import LMM
from limix._data import normalize_likelihood, conform_dataset
from limix.qtl._assert import assert_finite
from limix._display import session_block, session_line
""" custom function from limix package that calculates the variance components estiamtes """


def estimate(y_phe, lik, kin, marker_mat=None, verbose=True):
    """ estimate variance components """
    lik = normalize_likelihood(lik)
    lik_name = lik[0]
    with session_block("Heritability analysis", disable=not verbose):
        with session_line("Normalising input...", disable=not verbose):
            data = conform_dataset(y_phe, M=marker_mat, K=kin)
        y_phe = data["y"]
        marker_mat = data["M"]
        kin = data["K"]
        assert_finite(y_phe, marker_mat, kin)
        if kin is not None:
           # K = K / diag(K).mean()
            q_s = economic_qs(kin)
        else:
            q_s = None
        if lik_name == "normal":
            method = LMM(y_phe.values, marker_mat.values, q_s, restricted=True)
            method.fit(verbose=verbose)
        else:
            method = GLMMExpFam(y_phe, lik, marker_mat.values, q_s, n_int=500)
            method.fit(verbose=verbose, factr=1e6, pgtol=1e-3)
        v_g = method.scale * (1 - method.delta)
        v_e = method.scale * method.delta
        if lik_name == "bernoulli":
            v_e += pi * pi / 3
        v_v = var(method.mean())
        return v_g, v_v, v_e
