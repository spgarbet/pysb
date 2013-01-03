# 1991, Tyson, Modeling the cell division cycle: cdc2 and cyclin interactions

from pysb import *
from pysb.macros import *
from pysb.bng import generate_equations


Model()

Parameter('k1',  0.015)
#Parameter('k2',  0)
Parameter('k3',  200)
Parameter('k4',  180)
Parameter('kp4', 0.018)
Parameter('k5',  0)
Parameter('k6',  1)
Parameter('k7',  0.6)
Parameter('k8',  1e6)
Parameter('k9',  1e3)

Monomer('cyclin', ['Y', 'b'],{ 'Y': ['U','P']})
Monomer('cdc2',   ['Y', 'b'],{ 'Y': ['U','P']})

synthesize(cyclin(Y='U', b=None), k1)
#degrade(cyclin(Y='U', b=None), k2)

# One-way binding and phosophylation in one step (?)
Rule('Step3', cyclin(Y='U',b=None) + cdc2(Y='P',b=None) >> cyclin(Y='P', b=1) % cdc2(Y='P', b=1), k4)
#Rule('Step45', cyclin(Y='P', 1) % cdc2(Y='P', 1) <> cyclin(Y='P', 1) % cdc2(Y='U', 1), kp4, k5)
Rule('Step4', cyclin(Y='P', b=1) % cdc2(Y='P', b=1) >> cyclin(Y='P', b=1) % cdc2(Y='U', b=1), kp4)
Rule('Step4a', cyclin(Y='P', b=1) % cdc2(Y='P', b=1) + cyclin(Y='U', b=2) % cdc2(Y='P', b=2) >> 
               cyclin(Y='P', b=1) % cdc2(Y='U', b=1) + cyclin(Y='P', b=2) % cdc2(Y='U', b=2), k4)

Rule('Step6', cyclin(Y='P', b=1) % cdc2(Y='U', b=1) >> cyclin(Y='P',b=None) + cdc2(Y='U',b=None), k6)

degrade(cdc2(Y='U'), k7)
equilibrate(cdc2(Y='U', b=None), cdc2(Y='P', b=None), [k8, k9])

# This creates model.odes which contains the math
generate_equations(model)