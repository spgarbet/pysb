# Basic Libs
from pysb import *

# Integrate BNGL
from pysb.integrate import odesolve
from pylab import *

# Model
from pysb.examples.mass_action import model

# Smoldyn
from pysb.smoldynlib import *
from pysb.generator.smoldynlib import *
from pysb.generator.smoldyn import *

model.parameters['kf'].value = 0.0001

t = linspace(0, 1*60, 1*60+1)  
for a in [5, 50, 500, 5000]:
    y = odesolve(model, t, y0=[a, 5.0, 0.0])
    plot(t, y['Na']/a)

show()

model.parameters['kf'].value = 10
g = SmoldynlibGenerator(model)
g.generate_sim()

def stepsim(t):
    smolRunSimUntil(g.sim,t)
    return smolGetMoleculeCount(g.sim, "A", MolecState.SOLN)

result = map(stepsim,  t)
plot(t, array(result)/3000.0)
show()


#h = SmoldynGenerator(model)
#h.generate_content()
#print h.get_content()
a=3000.0
model.parameters['kf'].value = 0.0001
y = odesolve(model, t, y0=[a, 1000.0, 0.0])
plot(t, y['Na']/a)
plot(t, array(result)/3000.0)
show()