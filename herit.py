
def estimate(y, lik, K, M=None, verbose=True):
    from numpy_sugar.linalg import economic_qs
    from numpy import pi, var, diag
    from glimix_core.glmm import GLMMExpFam
    from glimix_core.lmm import LMM
    from limix._data._assert import assert_likelihood
    from limix._data import normalize_likelihood, conform_dataset 
    from limix.qtl._assert import assert_finite
    from limix._display import session_block, session_line
    lik = normalize_likelihood(lik)
    lik_name = lik[0]
    with session_block("Heritability analysis", disable=not verbose):
        with session_line("Normalising input...", disable=not verbose):
            data = conform_dataset(y, M=M, K=K)
        y = data["y"]
        M = data["M"]
        K = data["K"]
        assert_finite(y, M, K)
        if K is not None:
           # K = K / diag(K).mean()
            QS = economic_qs(K)
        else:
            QS = None
        if lik_name == "normal":
            method = LMM(y.values, M.values, QS, restricted=True)
            method.fit(verbose=verbose)
        else:
            method = GLMMExpFam(y, lik, M.values, QS, n_int=500)
            method.fit(verbose=verbose, factr=1e6, pgtol=1e-3)
        g = method.scale * (1 - method.delta)
        e = method.scale * method.delta
        if lik_name == "bernoulli":
            e += pi * pi / 3
        v = var(method.mean())
        return g , v , e 

