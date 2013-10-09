from pysb import *
from pysb.generator import *
from pysb.geometry import *

from pysb.smoldynlib import *

from pysb.generator.smoldynlib_spec import *
from pysb.bng import generate_equations

Model()

# For SmolDyn
Parameter('time_start', 0)
Parameter('time_stop',  5)
Parameter('time_step',  0.001)

Parameter('OSinit',  10000.0)
Parameter('OEinit',  5000.0)

Parameter('ok1',  83.0) # 0.05 * 1660
Parameter('okb1',  5.0)
Parameter('ok2',   1.0)

Monomer('OS', ['b'], difc=100)
Monomer('OE', difc=100)
Monomer('OC', ['b'], difc=100)
Monomer('OP', difc=100)


# Compute space based on 100 uM concentration for Sinit
# 6.022045e-2 / 100.0 => 1660.565
# Cube root(1660.565) => 11.84183 (must be an int, so using 12)
main=Compartment('Main', None, geometry=SquareSpace(12, [0, 0, 0]))

Initial(OS(b=None)**main, OSinit)
Initial(OE()**main, OEinit)
#Initial(OS(b=None),OSinit)
#Initial(OE(),OEinit)

Rule("r1", OS(b=None) + OE() <> OC(b=None) + OP(), ok1, okb1)
Rule("r2", OC(b=None) + OP() <> OP() + OE(), ok2,ok2)
Rule("bt", OS(b=None) + OC(b=None) >> OS(b=1) % OC(b=1), ok1)

m=Compartment('Membrane', Main,
              geometry=SphericalSurface(2, [0,0,0]),
              action=[ ["OS(b=None)",MolecState.ALL,SurfAction.REFLECT],
                       ["OE()",MolecState.ALL,SurfAction.REFLECT],
                       ["OC(b=None)",MolecState.ALL,SurfAction.REFLECT],
                       ["OP()",MolecState.ALL,SurfAction.PORT]
                     ])


c=Compartment('Cell',        m, geometry=SphericalSpace(  2, [0,0,0]))

generate_equations(model)

g = SmoldynlibGenerator(model)
g.generate_sim()

# This 
err = smolAddPort(g.sim, 'toODEfront', 'Membrane', PanelFace.FRONT)

if not err == 0: 
    raise Exception("ERROR adding front port to smoldyn simulation")

err = smolAddPort(g.sim, 'toODEback', 'Membrane', PanelFace.BACK)

if not err == 0:
    raise Excception("ERROR adding back port to smoldyn simulation")

#g.generate_content()
#print g.get_content()