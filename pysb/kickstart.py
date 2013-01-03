from pysb import *
from pysb.generator import *
from pysb.geometry import *

from pysb.generator.smoldyn import SmoldynGenerator

Model()

Parameter('display_size', 3)
Parameter('time_start', 0)
Parameter('time_stop',  250)
Parameter('time_step',  0.01)

Compartment('Main',      None,     geometry=SquareSpace(100, [0, 0]))
Compartment('Membrane1', Main,     geometry=SphericalSurface(10, [0,0]) )
Compartment('Cell1',     Membrane1, geometry=SphericalSpace(10, [0, 0]))

Monomer('ACA', {'Orient':['up', 'down']})  # Up, Down state (Down is active)
Monomer('ATP')
Monomer('cAMP')
Monomer('cAR1',{ 'Y': ['U','P']})

Parameter('aca_0', 30.0)
Parameter('car1_0', 30.0)
Initial(ACA(Orient="down") ** Membrane1, aca_0 )
Initial(cAR1(Y="U")  ** Membrane1, car1_0)


Parameter('k1', 1)
Parameter('k2', 1)
Rule('test_rule',
     (ACA(Orient="up")   + ATP() <>
      ACA(Orient="down") + ATP())**Cell1,
     k1, k2)


