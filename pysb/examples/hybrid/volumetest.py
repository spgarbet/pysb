from pysb import *
from pysb.geometry import *
from pysb.bng import generate_equations
from pysb.generator.smoldynlib import *
from pysb.smoldynlib import *
import numpy as np
from matplotlib.pyplot import *
from pysb.integrate import odesolve

# The following reaction system:
# S + E <> C > P + E

Model()

Parameter('Sinit',  100)
Parameter('Einit',  50)

Monomer('S')
Monomer('E')
Monomer('P')

# cell radius
radius = 0.5

# Main compartment's side
side = 1000

main =Compartment("Main", parent=None, dimension=3, geometry=SquareSpace(side,[0,0,0]))
mem  =Compartment('Membrane', main, dimension=2, geometry=SphericalSurface(radius, [0,0,0]))
cell =Compartment("Cell", parent=mem, dimension=3, geometry=SphericalSpace(radius,[0,0,0]))

Initial(S() ** cell, Sinit)
Initial(E() ** cell, Einit)

k1_val = 0.05
kb1_val = 5.0

Parameter('k1',  k1_val*cell.geometry.shape.volume) 
Parameter('kb1', kb1_val)

Rule("r1", S() ** cell + E() ** cell <> P() ** cell, k1, kb1)

Observable("enzyme_total",E())

generate_equations(model)

Parameter('time_start',0)
Parameter('time_step',0.1)
Parameter('time_stop',5)

mem.action=[ [repr(spec).replace('*','.'), MolecState.ALL, SurfAction.REFLECT] for spec in model.species]

def run():
    k1.value=k1_val/cell.size.value
    g = SmoldynlibGenerator(model)
    g.generate_sim()
    
    def stepsim(t):
        result = [smolGetMoleculeCount(g.sim, repr(spec).replace('*','.'), MolecState.SOLN) for spec in model.species]
        smolRunSimUntil(g.sim,t)
        return result
    
    t = np.linspace(0,5,(5./0.1)+1)
    result = map(stepsim,t)
    for j in range(len(result[0])):
        plot(t,[result[i][j] for i in range(len(result))])
    #show()
    g.free_sim()
    k1.value = k1_val
    y=odesolve(model,t)
    count=0
    for i in y.dtype.names:
        plot(t,y[i])
    show()
    return y
y=run()