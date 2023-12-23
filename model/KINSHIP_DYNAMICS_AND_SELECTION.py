# EVALUATE BOTH KINSHIP DYNAMICS & SELECTIVE PRESSURES ON HELPING AND HARMING (FOR MAIN & SUPPLEMENTARY RESULTS)

import os
from sympy import *
from mpi4py import MPI
import SELECTION as SEL
from pathlib import Path

comm = MPI.COMM_WORLD
if comm.size != 3:
    if comm.rank == 0:
        print("\n---------- SPAWN *THREE* MPI PROCESSES INSTEAD ----------\n")
    comm.Abort()

wd = os.getcwd() + '/'
os.chdir(wd)
dir = './res/'
Path(dir).mkdir(parents = True, exist_ok = True)

s = 5  # the size of smaller groups: 2 * s
l = 20 # the size of larger groups: 2 * l

pars = ([
    [Rational(.15), Rational(.85)], # df --- 0
    [Rational(.15), Rational(.85)], # dm --- 1
    [Rational(.1)], # muf --- 2
    [Rational(.1)], # mum --- 3
    [s, l], # nfa --- 4 ------------------------------------------------- ALPHA F
    [s, l], # nfb --- 5 ---------------------------------------------------------------- BETA F
    [s, l], # nma --- 6 ------------------------------------------------- ALPHA M
    [s, l], # nmb --- 7 ---------------------------------------------------------------- BETA M
    [Rational(.1), Rational(.5), Rational(.9)], # u --- 8
    [Rational(0), Rational(.82)], # m --- 9
    ['F', 'M'], # sexes --- 10
    ['ALPHA', 'BETA'], # patch_type --- 11
    ['FECUNDITY', 'MORTALITY'], # fec_or_mor --- 12
    [30] # lifespan --- 13
    ])
import itertools
pars = list(itertools.product(*pars))

##### male-biased dispersal with local mating
pars_ML = [p for p in pars if 
        (p[0] == .15 and p[1] == .85 and p[9] == .82 and p[4] == s and p[5] == s and p[6] == s and p[7] == s and p[8] == .5 and p[11] == 'ALPHA') or
        (p[0] == .15 and p[1] == .85 and p[9] == .82 and p[4] == s and p[5] == l and p[6] == s and p[7] == l) or
        (p[0] == .15 and p[1] == .85 and p[9] == .82 and p[4] == l and p[5] == l and p[6] == l and p[7] == l and p[8] == .5 and p[11] == 'BETA')]
o_fn_ML = dir + 'ML.txt'

##### female-biased dispersal with local mating
pars_FL = [p for p in pars if 
        (p[0] == .85 and p[1] == .15 and p[9] == .82 and p[4] == s and p[5] == s and p[6] == s and p[7] == s and p[8] == .5 and p[11] == 'ALPHA') or
        (p[0] == .85 and p[1] == .15 and p[9] == .82 and p[4] == s and p[5] == l and p[6] == s and p[7] == l) or
        (p[0] == .85 and p[1] == .15 and p[9] == .82 and p[4] == l and p[5] == l and p[6] == l and p[7] == l and p[8] == .5 and p[11] == 'BETA')]
o_fn_FL = dir + 'FL.txt'

##### bisexual philopatry with non-local mating
pars_BN = [p for p in pars if 
        (p[0] == .15 and p[1] == .15 and p[9] == 0 and p[4] == s and p[5] == s and p[6] == s and p[7] == s and p[8] == .5 and p[11] == 'ALPHA') or
        (p[0] == .15 and p[1] == .15 and p[9] == 0 and p[4] == s and p[5] == l and p[6] == s and p[7] == l) or
        (p[0] == .15 and p[1] == .15 and p[9] == 0 and p[4] == l and p[5] == l and p[6] == l and p[7] == l and p[8] == .5 and p[11] == 'BETA')]
o_fn_BN = dir + 'BN.txt'

par = [pars_ML, pars_FL, pars_BN]
o_fn = [o_fn_ML, o_fn_FL, o_fn_BN]
SEL.EVALUATE(o_fn[comm.rank], par[comm.rank])