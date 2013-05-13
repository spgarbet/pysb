# Basic Libs
from pysb import *

# Integrate BNGL
from pysb.integrate import odesolve
from pylab import *

# Model
from pysb.examples.michaelis_menten import model

# Smoldyn
from pysb.smoldynlib import *
from pysb.generator.smoldynlib import *
from pysb.generator.smoldyn import *

# ODE Solution

t = linspace(0, 1*5, 200*5+1)  
model.parameters['k1'].value  = 0.05
model.parameters['kb1'].value = 5
model.parameters['k2'].value  = 1
y = odesolve(model, t)


# SmolDyn

def stepsim(t):
    smolRunSimUntil(g.sim,t)
    return [smolGetMoleculeCount(g.sim, "S", MolecState.SOLN),
            smolGetMoleculeCount(g.sim, "E", MolecState.SOLN),
            smolGetMoleculeCount(g.sim, "C", MolecState.SOLN),
            smolGetMoleculeCount(g.sim, "P", MolecState.SOLN)]

adjustment = 1660   # Works with difc = 100
model.parameters['k1' ].value = 0.05 * adjustment
model.parameters['kb1'].value = 5
model.parameters['k2' ].value = 1


g = SmoldynlibGenerator(model)
g.generate_sim()


result = map(stepsim,  t)
result
r = array(result).T
plot(t, y['__s2']) # Blue
plot(t, y['__s1']) # Green
plot(t, y['__s0']) # Red
plot(t, y['__s3']) # Cyan
plot(t, r[0], 'r-')
plot(t, r[1], 'g-')
plot(t, r[2], 'b-')
plot(t, r[3], 'c-')
show()