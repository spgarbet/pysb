from pysb import *
from pysb.generator import *
from pysb.geometry import *

from pysb.smoldynlib import *

from pysb.generator.smoldynlib import *

Model()

Parameter('display_size', 3)

Parameter('time_start', 0)
Parameter('time_stop',  250)
Parameter('time_step',  0.01)

Parameter('aca_0',  30.0)
Parameter('car1_0', 30.0)

Parameter('k1',  0.02)
Parameter('k2',  2.30)
Parameter('k3',  0.03)
Parameter('k4',  2.30)
Parameter('k5',  0.20)
Parameter('k6',  0.20)
Parameter('k7', 20.00)

Monomer('ACA', {'Orient':['up', 'down']})  # Up, Down state (Down is active)
Monomer('ATP')
Monomer('cAMP')
Monomer('cAR1',{ 'Y': ['up','down']}) 

m=Compartment('Main',      None,     geometry=SquareSpace(100, [0, 0]))

cells = [[-20, 20],
         [ 20, 20],
         [ 20,-20],
         [-20,-20]]

# Create the 4 cells
for idx, loc in enumerate(cells):
    m=Compartment('Membrane%02d'%idx, Main, geometry=SphericalSurface(10, loc), action=[ ["ATP",PanelFace.BOTH,SurfAction.REFLECT] ])
    c=Compartment('Cell%02d'%idx,        m, geometry=SphericalSpace(  10, loc))
    Initial(ACA(Orient="down")**m, aca_0)
    Initial(cAR1(Y="up")**m, car1_0)
    Rule("r1%02d"%idx, None >> ATP(), k1, compartment = c)

Rule("r2", ACA(Orient="down") + ATP() >> ACA(Orient="down") + cAMP(), k2)
Rule("r3", cAMP() >> None, k3)
Rule("r4", cAR1(Y="up") + cAMP() >> cAR1(Y="down") + cAMP(), k4)
Rule("r5", cAR1(Y="down") >> cAR1(Y="up"), k5)
Rule("r6", ACA(Orient="down") >> ACA(Orient="up"), k6)
Rule("r7", ACA(Orient="up") +cAR1(Y="down") >> ACA(Orient="down") + cAR1(Y="down"), k7)


g = SmoldynlibGenerator(model)
g.generate_sim()

#g.generate_content()
#print g.get_content()