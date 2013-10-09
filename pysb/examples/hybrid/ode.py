from pysb import *

# The following reaction system:
# S + E <> C > P + E
# .05*10000*5000=2500000=2.5e6*.05=
# Breaks apart into the following individual unidirectional reactions:
# S + E > C
# C > S + E
# C > P + E

Model()

Parameter('Sinit',  100)
Parameter('Einit',  50)

Parameter('k1',  0.05) 
Parameter('kb1', 5.0)
Parameter('k2',  1.0)

Monomer('S', ['b'], difc=100)
Monomer('E', ['b'], difc=100)
#Monomer('C', difc=100)
Monomer('P', difc=100)

Initial(S(b=None), Sinit)
Initial(E(b=None), Einit)

Rule("r1", S(b=None) + E(b=None) <> S(b=1) % E(b=1), k1, kb1)
Rule("r2", S(b=1) % E(b=1) >> P() + E(b=None), k2)

Observable("enzyme_total",E())
