from pysb import *
from pysb.generator import *
from pysb.geometry import *

Model()

# For SmolDyn
Parameter('time_start', 0)
Parameter('time_stop',  5)
Parameter('time_step',  0.001)

Parameter('Sinit',  100.0)
Parameter('Einit',  50.0)

Parameter('k1',  0.05) # Type-o in paper 
Parameter('kb1', 5.0)
Parameter('k2',  1.0)

Monomer('S', difc=100)
Monomer('E', difc=100)
Monomer('C', difc=100)
Monomer('P', difc=100)


# Compute space based on 100 uM concentration for Sinit
# 6.022045e-2 / 100.0 => 1660.565
# Cube root(1660.565) => 11.84183 (must be an int, so using 12)
main=Compartment('Main', None, geometry=SquareSpace(12, [0, 0, 0]))

Initial(S()**main, Sinit)
Initial(E()**main, Einit)

Rule("r1", S() + E() <> C(), k1, kb1)
Rule("r2", C() >> P() + E(), k2)
