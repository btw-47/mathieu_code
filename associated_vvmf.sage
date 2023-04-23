## Make sure the imports are correct!
from weilrep import *
from weilrep.weilrep.weilrep_modular_forms_class import WeilRepModularForm

two_pi_i = complex(0.0, 2 * math.pi)
exp = cmath.exp


## numerical cutoff in case we want to do that
eps = 0.00001

@cached_function
def cached_exp(N):
    return exp(two_pi_i / N)

def trace(f, d):
    r"""
    The "trace map"  M_*(Gamma_0(N)) --> M_*(Gamma0(d)) for divisors d|N.
    """
    N = f.level()
    if N % d:
        raise ValueError('Invalid level')
    return sum(f.slash(matrix(x)) for x in Gamma0(N).coset_reps() if x[1][0] % d == 0)


class LinearCombo:
    r"""
    Linear combinations of eta products and Eisenstein series.
    """
    def __init__(self, forms, multiples, lvl = None):
        self.__forms = forms
        self.__multiples = multiples
        self.__qexp = sum(n * forms[i].qexp() for i, n in enumerate(multiples))
        if lvl is None:
            lvl = lcm([x.level() for x in forms])
        self.__level = lvl

    def __repr__(self):
        return str(self.__qexp)

    def forms(self):
        return self.__forms

    def level(self):
        return self.__level

    def multiples(self):
        return self.__multiples

    def qexp(self):
        return self.__qexp

    def weight(self):
        return self.forms()[0].weight()

    def __add__(self, other):
        try:
            forms = self.forms() + other.forms()
            multiples = self.multiples() + other.multiples()
        except AttributeError:
            forms = self.forms() + [other]
            multiples = self.multiples() + [1]
        return LinearCombo(forms, multiples)

    def __sub__(self, other):
        try:
            forms = self.forms() + other.forms()
            multiples = self.multiples() + [-x for x in other.multiples()]
        except AttributeError:
            forms = self.forms() + [other]
            multiples = self.multiples() + [-1]
        return LinearCombo(forms, multiples)

    def __mul__(self, N):
        return LinearCombo(self.forms(), [x * N for x in self.multiples()])

    def __div__(self, N):
        return LinearCombo(self.forms(), [x / N for x in self.multiples()])

    def __neg__(self):
        return LinearCombo(self.forms(), [-x for x in self.multiples()])

    __radd__ = __add__
    __rmul__ = __mul__
    __truediv__ = __div__

    def slash(self, M, lvl = None, num = False):
        r"""
        Apply the Petersson slash operator.
        """
        multiples = self.multiples()
        X = [x.slash(M, lvl=lvl, num = num) for x in self.forms()]
        if num:
            return sum(multiples[i] * x for i, x in enumerate(X))
        N = lcm([x.base_ring().gens()[0].multiplicative_order() for x in X])
        K0.<zetaN> = CyclotomicField(N)
        return sum(multiples[i] * x.base_extend(K0) for i, x in enumerate(X))

class LambdaN:
    r"""
    This is the modular form (N/24) * (N*E2(Nz) - E2(z)) of weight 2 and level Gamma0(N).
    """
    def __init__(self, N, prec, lvl = None):
        e2 = eisenstein_series_qexp(2, prec) * (-24)
        qexp = (N / 24) * ( N * e2.V(N) - e2 )
        if lvl is None:
            lvl = N
        self.__qexp = qexp
        self.__N = N
        self.__lvl = lvl
        self.__prec = prec

    def __repr__(self):
        return str(self.__qexp)

    def level(self):
        return self.__lvl

    def qexp(self):
        return self.__qexp

    def slash(self, M, lvl = None, num = False):
        r"""
        Apply the Petersson slash operator.
        """
        a, b, c, d = M.list()
        N = self.__N
        if lvl is None:
            lvl = self.level()
        prec = self.__prec
        if c % N:
            if num:
                K = CC
                zeta = cached_exp(N)
            else:
                K.<zeta> = CyclotomicField(lvl)
                zeta = zeta**(lvl/N)
            R.<q> = PuiseuxSeriesRing(K)
            f = R(0)
            g = gcd(c, N)
            for n in range(1, N * prec):
                s = 0
                for delt in divisors(n):
                    n_delt = n // delt
                    for j in range(1, N):
                        if n_delt % N == j*c % N:
                            s += delt * zeta**(j*d*delt % N)
                f -= s * q^(n / N)
            f = f + (N - 1) * eisenstein_series_qexp(2, prec)
            return f + (g^2 - 1) / 24
        return self.__qexp

    def weight(self):
        return 2

    def __mul__(self, n):
        return LinearCombo([self], [n])

    def __add__(self, other):
        return LinearCombo([self], [1]) + other

    def __sub__(self, other):
        return LinearCombo([self], [1]) - other

    def __div__(self, n):
        return LinearCombo([self], [1/n])

    def __neg__(self):
        return LinearCombo([self], [-1])

    __radd__ = __add__
    __rmul__ = __mul__
    __truediv__ = __div__

def eta_product(*args, **kwargs):
    r"""
    Redefine "eta_product".
    """
    h = EtaProduct(*args, **kwargs)
    return LinearCombo([h], [1], **kwargs)

class ModularConstant:
    r"""
    Constants as modular forms of weight 0.
    """
    def __init__(self, u, N, prec=None):
        if prec is not None:
            R.<q> = PowerSeriesRing(QQ)
            u = R(u) + O(q ** prec)
        self.__u = u
        self.__N = N
    def __repr__(self):
        return str(self.__u)
    def level(self):
        return self.__N
    def qexp(self):
        return self.__u
    def slash(self, M, **kwargs):
        return self.__u
    def weight(self):
        return 0

def associated_vvmf(X, lvl, t0 = None, verbose = False, num=False):
    if num:
        return associated_vvmf_num(X, lvl, t0=t0, verbose=verbose)
    weight = X[0].weight()
    N = lvl
    w = II(N)
    div = divisors(N)
    if N > 1:
    	K.<zeta> = CyclotomicField(N)
    else:
        K = QQ
        zeta = 1
    forms_dict = {}
    F, q = PuiseuxSeriesRing(K, 'q').objgen()
    for i in range(N):
        if verbose:
            print('i:', i)
            print('current zetas:', zeta.parent(), zeta0.parent())
        for j in range(N):
            m = gcd([i, j, N])
            if m == N:
                forms_dict[(i, j)] = F(X[-1].qexp())
            else:
                c, d = i/m, j/m
                while gcd(c, d) != 1:
                    if not d:
                        d += N/m
                    else:
                        c += N/m
                    if c > 10000:
                        raise RuntimeError
                _, b, a = xgcd(c, d)
                M = matrix([[a, -b], [c, d]])
                if M.determinant() != 1:
                    print(M)
                    assert False
                h = X[div.index(m)].slash(M, lvl = t0)
                mu, = h.base_ring().gens()
                lvl = lcm(lvl, mu.multiplicative_order())
                K0.<zeta0> = CyclotomicField(lvl)
                forms_dict[(i, j)] = h.base_extend(K0)
    X = []
    nl = w.norm_list()
    for ii, g in enumerate(w.ds()):
        i, k = N * g
        f = 0 * zeta
        for j in range(N):
            u = zeta**(j*k)
            h = forms_dict[(i, j)]
            try:
                f += u * h
            except TypeError:
                lvl = lcm(lvl, h.base_ring().gens()[0].multiplicative_order())
                K0.<zeta0> = CyclotomicField(lvl)
                print('Extended cyclotomy level to %s'%lvl)
                f += u * h.base_extend(K0)
        f = f / N
        if f:
            q, = f.parent().gens()
            f = q**(-nl[ii]) * f
            f = f + O(q ** floor(f.prec()))
        try:
            X.append([g, nl[ii], f.laurent_part()])
        except AttributeError:
            X.append([g, nl[ii], f])
    return WeilRepModularForm(weight, w.gram_matrix(), X, weilrep = w)

def associated_vvmf_num(X, lvl, t0 = None, verbose = False):
    weight = X[0].weight()
    N = lvl
    w = II(N)
    div = divisors(N)
    K = CC
    zeta_N = cached_exp(N)
    forms_dict = {}
    F, q = PuiseuxSeriesRing(K, 'q').objgen()
    for i in range(N):
        if verbose:
            print(i)
        for j in range(N):
            m = gcd([i, j, N])
            if m == N:
                forms_dict[(i, j)] = F(X[-1].qexp())
            else:
                c, d = i/m, j/m
                while gcd(c, d) != 1:
                    if not d:
                        d += N/m
                    else:
                        c += N/m
                    if c > 10000:
                        raise RuntimeError
                _, b, a = xgcd(c, d)
                M = matrix([[a, -b], [c, d]])
                if M.determinant() != 1:
                    print(M)
                    assert False
                h = X[div.index(m)].slash(M, lvl = t0, num=True)
                forms_dict[(i, j)] = h
    X = []
    nl = w.norm_list()
    for ii, g in enumerate(w.ds()):
        i, k = N * g
        f = 0.0
        for j in range(N):
            u = zeta_N**(j*k)
            h = forms_dict[(i, j)]
            f += u * h
        f = f / N
        if all(abs(x) < eps for x in f.list()):
            f = f * 0.0
        if f:
            q, = f.parent().gens()
            f = q**(-nl[ii]) * f
            f = f + O(q ** floor(f.prec()))
        try:
            X.append([g, nl[ii], f.laurent_part()])
        except AttributeError:
            X.append([g, nl[ii], f])
    return WeilRepModularForm(weight, w.gram_matrix(), X, weilrep = w)