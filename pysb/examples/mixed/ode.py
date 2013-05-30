from pysb import *

Model()

Parameter('ISinit',  100.0)
Parameter('IEinit',  0.0001) # Had to put some in or it's eliminated in the odes

Parameter('ik1',  0.05) 
Parameter('ikb1', 5.0)
Parameter('ik2',  1.0)

Monomer('IS')
Monomer('IE') # This will come from Smoldyn
Monomer('IC')
Monomer('IP')


# Compute space based on 100 uM concentration for Sinit
# 6.022045e-2 / 100.0 => 1660.565
# Cube root(1660.565) => 11.84183 (must be an int, so using 12)

Initial(IS(), ISinit)
Initial(IE(), IEinit)   # Get from Smoldyn Port

Rule("ir1", IS() + IE() <> IC(), ik1, ikb1)
Rule("ir2", IC() >> IP() + IE(), ik2)