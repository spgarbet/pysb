from pysb import *
from pysb.generator import *
from pysb.geometry import *



Model()

# For SmolDyn
Parameter('time_start', 0)
Parameter('time_stop',  60)
Parameter('time_step',  0.001)

Parameter('Ainit',  3000.0)
Parameter('Binit',  3000.0)

Parameter('kf',  0.001)

Monomer('A')
Monomer('B')
Monomer('C')

Observable('Na', A())

main=Compartment('Main', None, geometry=SquareSpace(1000, [0, 0, 0]))

Initial(A()**main, Ainit)
Initial(B()**main, Binit)

Rule("r1", A() + B() >> B() + C(), kf)


