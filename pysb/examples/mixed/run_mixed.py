from pysb.integrate import odesolve
from pylab import *
from pysb.smoldynlib import *

from pysb.examples.mixed.smoldyn import model,g as smoldyn,g
from pysb.examples.mixed.ode     import model   as ode

t = linspace(0, 1*0.025, 200*0.25+1)

# Idea for final interface

#result = smoldyn_mixed_simulation(smoldyn, ode, t, {"toODEfront": [("OP", "IE")], "toODEback": [("OP", "IE")]}, None )

# Main Loop
smolRunSimUntil(g.sim, 0.0)
result = [[100.00, 0.0, 0.0, 0.0,
           smolGetMoleculeCount(g.sim, "OS", MolecState.SOLN),
           smolGetMoleculeCount(g.sim, "OE", MolecState.SOLN),
           smolGetMoleculeCount(g.sim, "OC", MolecState.SOLN),
           smolGetMoleculeCount(g.sim, "OP", MolecState.SOLN)]]
          
for sec in t[1:]:
    y0 = result[-1]
    y = odesolve(ode, [0, 0.005], y0=[y0[0], y0[1], y0[2], y0[3]])
    smolRunSimUntil(g.sim, sec)
    extra =  smolGetPortMolecules(g.sim, 'toODEfront', "OP", MolecState.ALL, 1)
    extra += smolGetPortMolecules(g.sim, 'toODEback',  "OP", MolecState.ALL, 1)
    result.append([  y[1][0],
                     y[1][1] + extra,
                     y[1][2],
                     y[1][3],
                     smolGetMoleculeCount(g.sim, "OS", MolecState.SOLN),
                     smolGetMoleculeCount(g.sim, "OE", MolecState.SOLN),
                     smolGetMoleculeCount(g.sim, "OC", MolecState.SOLN),
                     smolGetMoleculeCount(g.sim, "OP", MolecState.SOLN)])
    print sec



r = array(result).T
for i in range(4):
    plot(t, r[i]+r[i+4]) # Blue

show()