from pysb import *
from pysb.geometry import *
from pysb.bng import generate_equations
import numpy as np
from pysb.integrate import odesolve
from pysb import macros

# This program fails, but I think it should not.  The problem is that the model has compartments,
# but they get associated to the ComplexPattern, not the MonomerPattern.  However,
# MonomerPattern.is_concrete() returns false if the MonomerPattern doesn't have a compartment,
# when the model does.  So for now, in this version, synthesis will not work in PySB with compartments

# From macros example page on pysb.org
Model()
Monomer('A', sites=['x','y'], site_states={'y':['e','f']})

# Cell radius and Main compartment side length
radius = 0.5
side = 1000

# Define compartments (not even giving the membrane a reflecting action for now)
main =      Compartment("Main",     parent=None,     dimension=3, geometry=SquareSpace(side,[0,0,0]))
membrane  = Compartment('Membrane', parent=main,     dimension=2, geometry=SphericalSurface(radius, [0,0,0]))
cell =      Compartment("Cell",     parent=membrane, dimension=3, geometry=SphericalSpace(radius,[0,0,0]))

# Synthesis rule (a la example on pysb.org)
macros.synthesize(A(x=None,y='e') ** cell, 5)

# BNGL quits because A(x=None,y='e') ** cell is "not concrete"
generate_equations(model)