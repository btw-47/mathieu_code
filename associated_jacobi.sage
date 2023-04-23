## Make sure the imports are correct!
from weilrep.weilrep.weilrep_modular_forms_class import WeilRepModularForm



def associated_jacobi(f, num = False):
    r"""
    Compute the sequence of Jacobi forms of lattice index "L"
    attached to a vector-valued modular form for the Weil representation for L+U(N).

    NOTE: The result is actually a sequence of vector-valued functions.
    The Jacobi forms can be recovered from them by typing
    [x.jacobi_form() for x in X]
    if 'X' is the output of this function.
    """
    w = f.weilrep()
    f_comp = f.components()
    N = w._N()
    if num:
        zeta = exp(2 * pi * I / N).n()
    else:
        zeta, = CyclotomicField(N).gens()
    w0 = WeilRep(w.gram_matrix()[1:-1, 1:-1])
    nl = w0.norm_list()
    ds = w0.ds()
    L = []
    for d in divisors(N):
        X = [None for _ in ds]
        for i, g in enumerate(ds):
            h = 0
            for c in range(N):
                h += (zeta ** (-d * c)) * f_comp[tuple([0] + list(g) + [c / N])]
            X[i] = g, nl[i], h
        L.append(WeilRepModularForm(f.weight(), w0.gram_matrix(), X, weilrep = w0))
    return L