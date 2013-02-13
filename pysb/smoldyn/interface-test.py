from pysb.smoldynlib import *

# dim 2
# boundaries 0 -50 50
# boundaries 1 -50 50
lbound = (c_double * 2)(-50.0, -50.0)
ubound = (c_double * 2)(50.0, 50.0)
s = smolNewSim(2, pointer(lbound), pointer(ubound))

# species   ACA  ATP  cAMP  cAR1
# difc ACA(all) 1
# difc ATP(all) 1
# difc cAMP(all) 1
# difc cAR1(all) 1
species = ["ACA", "ATP", "cAMP", "cAR1"]
for sp in species:
    smolAddSpecies(s, sp, "")
    smolSetSpeciesMobility(s, sp, MolecState.ALL, 1, 0, 0)


# color ACA(all) maroon
# color ATP(all) red
# color cAMP(all) orange
# color cAR1(all) yellow
# display_size all 3
display_size = 3
maroon = (c_double * 3)(195.0/255.0, 33.0/255.0, 72.0/255.0)
smolSetMoleculeStyle(s, "ACA", MolecState.ALL, display_size, maroon)

red = (c_double * 3)(1.0, 0.0, 0.0)
smolSetMoleculeStyle(s, "ATP", MolecState.ALL, display_size, red)

orange = (c_double * 3)(1.0, 165.0/255.0, 0.0)
smolSetMoleculeStyle(s, "cAMP", MolecState.ALL, display_size, orange)

yellow = (c_double * 3)(1.0, 1.0, 51.0/255.0)
smolSetMoleculeStyle(s, "cAR1", MolecState.ALL, display_size, yellow)

# 
# time_start 0.000000
# time_stop 250.000000
# time_step 0.010000
smolSetSimTimes(s, 0.0, 250.0, 0.01)

# frame_thickness 1
# ?????

#start_surface Membrane00
#action ATP both reflect
#panel sph -20 20 10 20 20
#thickness 1
#color both grey 0.5
#end_surface
smolAddSurface(s, "Membrane00")
smolSetSurfaceAction(s, "Membrane00", PanelFace.BOTH, "ATP", MolecState.ALL, SurfAction.REFLECT)
params = (c_double * 5)(-20.0, 20.0, 10.0, 20.0, 20.0)
smolAddPanel(s, "Membrane00", PanelShape.SPHERE, "", "", params)
grey = (c_double * 3)(0.5, 0.5, 0.5)
smolSetSurfaceStyle(s, "Membrane00", PanelFace.BOTH, DrawMode.NONE, 1, grey, -1, -1, 0.5)

#start_compartment Cell00
#surface Membrane00
#point  -20 20
#end_compartment
smolAddCompartment(s, "Cell00")
smolAddCompartmentSurface(s, "Cell00", "Membrane00")
point = (c_double * 2)(-20.0, 20.0)
smolAddCompartmentPoint(s, "Cell00", point)


# start_surface Membrane01
# action ATP both reflect
# panel sph 20 20 10 20 20
# thickness 1
# color both grey 0.5
# end_surface
# 
# start_compartment Cell01
# surface Membrane01
# point  20 20
# end_compartment
# 
# start_surface Membrane02
# action ATP both reflect
# panel sph 20 -20 10 20 20
# thickness 1
# color both grey 0.5
# end_surface
# 
# start_compartment Cell02
# surface Membrane02
# point  20 -20
# end_compartment
# 
# start_surface Membrane03
# action ATP both reflect
# panel sph -20 -20 10 20 20
# thickness 1
# color both grey 0.5
# end_surface
# 
# start_compartment Cell03
# surface Membrane03
# point  -20 -20
# end_compartment

#surface_mol 30 ACA(down) Membrane00 all all
#surface_mol 30 cAR1(up) Membrane00 all all
smolAddSurfaceMolecules(s, "ACA", MolecState.DOWN, 30, "Membrane00", PanelShape.ALL, "all", 0)
smolAddSurfaceMolecules(s, "cAR1", MolecState.UP, 30, "Membrane00", PanelShape.ALL, "all", 0)


#surface_mol 30 ACA(down) Membrane01 all all
#surface_mol 30 cAR1(up) Membrane01 all all
#surface_mol 30 ACA(down) Membrane02 all all
#surface_mol 30 cAR1(up) Membrane02 all all
#surface_mol 30 ACA(down) Membrane03 all all
#surface_mol 30 cAR1(up) Membrane03 all all


#reaction_cmpt Cell00 r100 0 -> ATP 0.02
# HOW TO DEAL WITH POINTER TO POINTER LISTS OF PRODUCTS!!!!!
species = (c_char_p * 1)("ATP")
#states  = pointer(c_int(MolecState.NONE))
states = c_int(MolecState.NONE)

smolAddReaction(s, "r100", "", MolecState.ALL, "", MolecState.ALL, 1,  byref(species), byref(states), 0.02)
smolSetReactionRegion(s, "r100", "Cell00", "")
#reaction_cmpt Cell01 r101 0 -> ATP 0.02
#reaction_cmpt Cell02 r102 0 -> ATP 0.02
#reaction_cmpt Cell03 r103 0 -> ATP 0.02

# reaction r2 ACA(down) + ATP(bsoln) -> ACA(down) + cAMP(bsoln) 2.3
species = (c_char_p * 2)("ACA", "cAMP")
states = (c_int * 2)(MolecState.DOWN, MolecState.BSOLN)

smolAddReaction(s, "r2", "ACA", MolecState.DOWN, "ATP", MolecState.BSOLN, 2, species, states ,2.3)


reaction r3 cAMP -> 0 0.03
reaction r4 cAR1(up) + cAMP(bsoln) -> cAR1(down) + cAMP(bsoln) 2.3
reaction r5 cAR1(down) -> cAR1(up) 0.2
reaction r6 ACA(down) -> ACA(up) 0.2
reaction r7 ACA(up) + cAR1(down) -> ACA(down) + cAR1(down) 20.0



smolDisplaySim(s)
smolFreeSim(s)