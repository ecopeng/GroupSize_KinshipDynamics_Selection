# CALCULATE INDIVIDUALS' SEX- AND AGE-SPECIFIC AVERAGE LOCAL GENETIC RELATEDNESS TO OTHERS

import os
from sympy import *
import RELATEDNESS as R
from pathlib import Path

wd = os.getcwd() + '/'
os.chdir(wd)
dir = './res/'
Path(dir).mkdir(parents = True, exist_ok = True)
fn = dir + 'r.txt'

# parameter values
nfa = 5        # N females in alpha-type groups
nma = 5        # N males in alpha-type groups
nfb = 20       # N females in beta-type groups
nmb = 20       # N males in beta-type groups
u = .5         # % of alpha-type groups
d_f = .15      # female dispersal rate
d_m = .85      # male dispersal rate
m = .82        # rate of local mating
mu_f = .1      # mortality rate females
mu_m = .1      # mortality rate males
lifespan = 30  # lifespan

# predict the patterns of kinship dynamics for both sexes
R.CALCULATE(fn, lifespan, u, nfa, nma, nfb, nmb, d_f, d_m, m, mu_f, mu_m)