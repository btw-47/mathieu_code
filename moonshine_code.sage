#############
## A typical use:
## for x in conj_classes:
##    compute_forms(x)
#############


##########
## Imports
from weilrep import *
from eta_products import *
load('associated_vvmf.sage')
load('associated_jacobi.sage')
##
##########


##########
## Global settings
prec = 10
##
##########


##########
## The basic jacobi forms phi0 and phi2 (i.e. phi_{-2})
w0 = WeilRep([[2]])
phi0, = w0.nearly_holomorphic_modular_forms_basis(-1/2, 1/4, prec)
phi2, = w0.nearly_holomorphic_modular_forms_basis(-5/2, 1/4, prec)
##
##########




##########
## Data. We construct 5 lists and/or dictionaries:
## 1) a list 'conj_classes' of names of the conjugacy classes in the Mathieu group. (Some classes such as 7AB that yield the same information are combined)
## 2) a dictionary 'levels' assigning an integer to each conjugacy class
## 3) a dictionary 'powers' assigning a list of conjugacy classes to each conjugacy class. These are the classes of g^d,
##    where 'g' is any element in the given class and 'd' runs through the divisors of g's level.
## 4) a dictionary 'chi'
## 5) a dictionary 'Tg'
##########



conj_classes = ['1A', '2A', '2B', '3A', '3B', '4A', '4B', '4C', '5A', '6A',
'6B', '7AB', '8A', '10A', '11A', '12A', '12B', '14AB', '15AB', '21AB', '23AB']

levels = {
'1A': 1,
'2A': 2,
'2B': 4,
'3A': 3,
'3B': 9,
'4A': 8,
'4B': 4,
'4C': 16,
'5A': 5,
'6A': 6,
'6B': 36,
'7AB': 7,
'8A': 8,
'10A': 20,
'11A': 11,
'12A': 24,
'12B': 144,
'14AB': 14,
'15AB': 15,
'21AB': 63,
'23AB': 23,
}

powers = {
'1A':  ['1A'],
'2A':  ['2A', '1A'],
'2B':  ['2B', '1A', '1A'],
'3A':  ['3A', '1A'],
'3B':  ['3B', '1A', '1A'],
'4A':  ['4A', '2A', '1A', '1A'],
'4B':  ['4B', '2A', '1A'],
'4C':  ['4C', '2B', '1A', '1A', '1A'],
'5A':  ['5A', '1A'],
'6A':  ['6A', '3A', '2A', '1A'],
'6B':  ['6B', '3B', '2B', '3B', '1A', '2B', '1A', '1A', '1A'],
'7AB': ['7AB', '1A'],
'8A':  ['8A', '4B', '2A', '1A'],
'10A': ['10A', '5A', '5A', '2B', '1A', '1A'],
'11A': ['11A', '1A'],
'12A': ['12A', '6A', '4A', '3A', '2A', '3A', '1A', '1A'],
'12B': ['12B', '6B', '4C', '3B', '2B', '3B', '4C', '1A', '3B', '2B', '1A', '1A', '1A', '1A', '1A'],
'14AB':['14AB', '7AB', '2A', '1A'],
'15AB':['15AB', '5A', '3A', '1A'],
'21AB':['21AB', '7AB', '3B', '7AB', '1A', '1A'],
'23AB':['23AB', '1A'],
}

chi = {
'1A':  24,
'2A':  8,
'2B':  0,
'3A':  6,
'3B':  0,
'4A':  0,
'4B':  4,
'4C':  0,
'5A':  4,
'6A':  2,
'6B':  0,
'7AB': 3,
'8A':  2,
'10A': 0,
'11A': 2,
'12A': 0,
'12B': 0,
'14AB':1,
'15AB':1,
'21AB':0,
'23AB':1,
}

Tg = {
'1A':  0 * LambdaN(1, prec),
'2A':  16 * LambdaN(2, prec),
'2B':  -24 * LambdaN(2, prec) + 8 * LambdaN(4, prec),
'3A':  6 * LambdaN(3, prec),
'3B':  2 * eta_product([[1, 6], [3, -2]], prec, lvl = levels['3B']),
'4A':  2 * eta_product([[2, 8], [4, -4]], prec, lvl = levels['4A']),
'4B':  -4 * LambdaN(2, prec) + 4 * LambdaN(4, prec),
'4C':  2 * eta_product([[1, 4], [2, 2], [4, -2]], prec, lvl = levels['4C']),
'5A':  2 * LambdaN(5, prec),
'6A':  -2 * LambdaN(2, prec) - 2 * LambdaN(3, prec) + 2 * LambdaN(6, prec),
'6B':  2 * eta_product([[1, 2], [2, 2], [3, 2], [6, -2]], prec, lvl = levels['6B']),
'7AB': 1*LambdaN(7, prec),
'8A':  -1*LambdaN(4, prec) + LambdaN(8, prec),
'10A': 2 * eta_product([[1, 3], [2, 1], [5, 1], [10, -1]], prec, lvl = levels['10A']),
'11A': (2/5) * LambdaN(11, prec) - (22/5) * eta_product([[1, 2], [11, 2]], prec, lvl = levels['11A']),
'12A': 2 * eta_product([[1, 3], [2, -1], [3, -1], [4, 2], [6, 3], [12, -2]], prec, lvl = levels['12A']),
'12B': 2 * eta_product([[1, 4], [2, -1], [4, 1], [6, 1], [12, -1]], prec, lvl = levels['12B']),
'14AB':(-LambdaN(2, prec) - LambdaN(7, prec) + LambdaN(14, prec) - 14 * eta_product([[1, 1], [2, 1], [7, 1], [14, 1]], prec, lvl = levels['14AB']) ) / 3,
'15AB':(-LambdaN(3, prec) - LambdaN(5, prec) + LambdaN(15, prec) - 15 * eta_product([[1, 1], [3, 1], [5, 1], [15, 1]], prec, lvl = levels['15AB']) ) / 4,
'21AB':(7 * eta_product([[1, 3], [3, -1], [7, 3], [21, -1]], prec, lvl = levels['21AB']) - eta_product([[1, 6], [3, -2]], prec, lvl = levels['21AB']) ) /3,
'23AB':(LambdaN(23, prec) - 138 * eta_product([[1, 2], [23, 2]], prec, lvl = levels['23AB']) - 23 * eta_product([[1, 3], [2, -1], [23, 3], [46, -1]]) - 23 * 4 * eta_product([[1, 1], [2, 1], [23, 1], [46, 1]], prec, lvl = levels['23AB']) - 23 * 4 * eta_product([[2, 2], [46, 2]], prec, lvl = levels['23AB']) )/11,
}



##########
## The main script.
## For the conjugacy class "x", we compute the vector-valued modular form attached to the input function,
## to the precision specified above (default 10)
## its associated sequence of Jacobi forms
## and its principal part
## The forms F0 and F2 are attached to the Gram matrix [[0, N], [N, 0]] and
## the forms phi0, phi2 are attached to the Gram matrix [[2]].
## Line 7 F.conjugate(...) passes from Gram matrix [[0, N, 0], [N, 0, 0], [0, 0, 2]] to [[0, 0, N], [0, 2, 0], [N, 0, 0]].
## The result is saved in a txt file with filename "jacobi_forms_(...).txt".
##########


def compute_forms(x, **kwargs):
    lvl = levels[x]
    H2 = [Tg[y] for y in powers[x]]
    H0 = [ModularConstant(chi[y] / 12, lvl, prec=prec) for y in powers[x]]
    F2 = associated_vvmf(H2, lvl, **kwargs)
    F0 = associated_vvmf(H0, lvl, **kwargs)
    F = F0 * phi0 + F2 * phi2
    F = F.conjugate(matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]]))
    F._WeilRepModularForm__weilrep = w0 + II(lvl)
    j = associated_jacobi(F)
    with open('jacobi_forms_%s.txt'%str(x), 'a') as f:
        s = str(x) + '\n\nAssociated Jacobi forms:\n' + '\n\n'.join('phi_%s: %s'%(d, str(j[i].jacobi_form())) for i, d in enumerate(divisors(lvl))) + '\n\nPrincipal part: %s\n\n'%F.principal_part()+'-'*80+'\n\n'
        print(s)
        f.write(s)
    return F