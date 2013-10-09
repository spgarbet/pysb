from ode import model as base_model
import numpy
from pylab import *
from hybrid import run_hybrid
import itertools

t = numpy.linspace(0, 5, 1001)


smoldyn_species = [
    ('S(b=1) % E(b=1)',100),
    ('P()',100)
]

result = run_hybrid(base_model, t, smoldyn_species)

ls={}
ls['_s0'] = 'r-'
ls['_s1'] = 'b-'
ls['_s2'] = 'g-'
ls['_s3'] = 'c-'
ls['enzyme_total'] = 'k-'

for i in result.dtype.names:
    plot(t,result[i],ls[i])
    
show()