from ctypes import *
from ctypes.util import find_library

class MolecState(object):
    (SOLN,FRONT,BACK,UP,DOWN,BSOLN,ALL,NONE,SOME) = range(9)

class PanelFace(object):
    (FRONT, BACK, NONE, BOTH) = range(4)

class PanelShape(object):
    (RECTANGLE, TRIANGLE, SPHERE, CYLINDER, HEMISPHERE, DISK, ALL, NONE) = range(8)

class SurfAction(object):
    (REFLECT, TRANS, ABSORB, JUMP, PORT, MULT, NO, NONE, ADSORB, REVDES, IRREVDES, FLIP) = range(12)

class DrawMode(object):
    (NO, VERT, EDGE, VE, FACE, VF, EF, VEF, NONE) = range(9)

class RevParam(object):
    (NONE,IRREV,CONFSPREAD,BOUNCE,PGEM,PGEMMAX,PGEMMAXW,RATIO,
     UNBINDRAD,PGEM2,PGEMMAX2,RATIO2,OFFSET,FIXED) = range(14);


# Load Smoldyn
path = find_library("smoldyn_shared")
smoldyn = cdll.LoadLibrary(path)

# Define Get Version
smolGetVersion = smoldyn.smolGetVersion
smolGetVersion.restype = c_double

# Error Handling
smolSetDebugMode = smoldyn.smolSetDebugMode
smolSetDebugMode.restype = None
smolSetDebugMode.argtypes = [c_int]

# Simulation Structure
smolNewSim = smoldyn.smolNewSim
smolNewSim.restype = c_void_p
smolNewSim.argtypes = [c_int, c_void_p, c_void_p]

smolUpdateSim = smoldyn.smolUpdateSim
smolUpdateSim.restype = c_int
smolUpdateSim.argtypes = [c_void_p]

smolRunTimeStep = smoldyn.smolRunTimeStep
smolRunTimeStep.restype = c_int
smolRunTimeStep.argtypes = [c_void_p]

smolRunSim = smoldyn.smolRunSim
smolRunSim.restype = c_int
smolRunSim.argtypes = [c_void_p]

smolRunSimUntil = smoldyn.smolRunSimUntil
smolRunSimUntil.restype = c_int
smolRunSimUntil.argtypes = [c_void_p, c_double]

smolFreeSim = smoldyn.smolFreeSim
smolFreeSim.restype = c_int
smolFreeSim.argtypes = [c_void_p]

smolDisplaySim = smoldyn.smolDisplaySim
smolDisplaySim.restype = c_int
smolDisplaySim.argtypes = [c_void_p]

# Simulation Settings
smolSetSimTimes = smoldyn.smolSetSimTimes
smolSetSimTimes.restype = c_int
smolSetSimTimes.argtypes = [c_void_p, c_double, c_double, c_double]

smolSetTimeStart = smoldyn.smolSetTimeStart
smolSetTimeStart.restype = c_int
smolSetTimeStart.argtypes = [c_void_p, c_double]

smolSetTimeStop = smoldyn.smolSetTimeStop
smolSetTimeStop.restype = c_int
smolSetTimeStop.argtypes = [c_void_p, c_double]

smolSetTimeNow = smoldyn.smolSetTimeNow
smolSetTimeNow.restype = c_int
smolSetTimeNow.argtypes = [c_void_p, c_double]

smolSetTimeStep = smoldyn.smolSetTimeStep
smolSetTimeStep.restype = c_int
smolSetTimeStep.argtypes = [c_void_p, c_double]

smolSetRandomSeed = smoldyn.smolSetRandomSeed
smolSetRandomSeed.restype = c_int
smolSetRandomSeed.argtypes = [c_void_p, c_long]

smolSetPartitions = smoldyn.smolSetPartitions
smolSetPartitions.restype = c_int
smolSetPartitions.argtypes = [c_void_p, c_char_p, c_double]

# Runtime commands (Output values to file)
smolSetOutputPath = smoldyn.smolSetOutputPath
smolSetOutputPath.restype = c_int
smolSetOutputPath.argtypes = [c_void_p, c_char_p]

smolAddOutputFile = smoldyn.smolAddOutputFile
smolAddOutputFile.restype = c_int
smolAddOutputFile.argtypes = [c_void_p, c_char_p, c_int, c_int]

smolAddCommand = smoldyn.smolAddCommand
smolAddCommand.restype = c_int
smolAddCommand.argtypes = [c_void_p, c_char, c_double, c_double, c_double, c_double, c_char_p]

smolAddCommandFromString = smoldyn.smolAddCommandFromString
smolAddCommandFromString.restype = c_int
smolAddCommandFromString.argtypes = [c_void_p, c_char_p]


# Molecules
smolAddSpecies = smoldyn.smolAddSpecies
smolAddSpecies.restype = c_int
smolAddSpecies.argtypes = [c_void_p, c_char_p, c_char_p]
#extern "C" int            smolGetSpeciesIndex(simptr sim,const char *species);
#extern "C" int            smolGetSpeciesIndexNT(simptr sim,const char *species);
#extern "C" char*          smolGetSpeciesName(simptr sim,int speciesindex,char *species);

smolSetSpeciesMobility = smoldyn.smolSetSpeciesMobility
smolSetSpeciesMobility.restype = c_int
smolSetSpeciesMobility.argtypes = [c_void_p, c_char_p, c_int, c_double, c_void_p, c_void_p]

#//?? needs function smolSetSpeciesSurfaceDrift
#extern "C" enum ErrorCode smolAddMolList(simptr sim,const char *mollist);
#extern "C" int            smolGetMolListIndex(simptr sim,const char *mollist);
#extern "C" int            smolGetMolListIndexNT(simptr sim,const char *mollist);
#extern "C" char*          smolGetMolListName(simptr sim,int mollistindex,char *mollist);
#extern "C" enum ErrorCode smolSetMolList(simptr sim,const char *species,enum MolecState state,const char *mollist);
#extern "C" enum ErrorCode smolSetMaxMolecules(simptr sim,int maxmolecules);

#extern "C" enum ErrorCode smolAddSolutionMolecules(simptr sim,const char *species,int number,double *lowposition,double *highposition);
smolAddSolutionMolecules = smoldyn.smolAddSolutionMolecules
smolAddSolutionMolecules.restype = c_int
smolAddSolutionMolecules.argtypes = [c_void_p, c_char_p, c_int, c_void_p, c_void_p]

smolAddCompartmentMolecules = smoldyn.smolAddCompartmentMolecules
smolAddCompartmentMolecules.restype = c_int
smolAddCompartmentMolecules.argtypes = [c_void_p, c_char_p, c_int, c_char_p]

smolAddSurfaceMolecules = smoldyn.smolAddSurfaceMolecules
smolAddSurfaceMolecules.restype = c_int
smolAddSurfaceMolecules.argtypes = [c_void_p, c_char_p, c_int, c_int, c_char_p, c_int, c_char_p, c_void_p]

#extern "C" int            smolGetMoleculeCount(simptr sim,const char *species,enum MolecState state);
smolGetMoleculeCount = smoldyn.smolGetMoleculeCount
smolGetMoleculeCount.restype = c_int
smolGetMoleculeCount.argtypes = [c_void_p, c_char_p, c_int]

smolSetMoleculeStyle = smoldyn.smolSetMoleculeStyle
smolSetMoleculeStyle.restype = c_int
smolSetMoleculeStyle.argtypes = [c_void_p, c_char_p, c_int, c_double, c_void_p]

# Surfaces
smolSetBoundaryType = smoldyn.smolSetBoundaryType
smolSetBoundaryType.restype = c_int
smolSetBoundaryType.argtypes = [c_void_p, c_int, c_int, c_char]

smolAddSurface = smoldyn.smolAddSurface
smolAddSurface.restype = c_int
smolAddSurface.argtypes = [c_void_p, c_char_p]

smolGetSurfaceIndex = smoldyn.smolGetSurfaceIndex
smolGetSurfaceIndex.restype = c_int
smolGetSurfaceIndex.argtypes = [c_void_p, c_char_p]

smolGetSurfaceIndexNT = smoldyn.smolGetSurfaceIndexNT
smolGetSurfaceIndexNT.restype = c_int
smolGetSurfaceIndexNT.argtypes = [c_void_p, c_char_p]

smolGetSurfaceName = smoldyn.smolGetSurfaceName
smolGetSurfaceName.restype = c_char_p
smolGetSurfaceName.argtypes = [c_void_p, c_int, c_char_p]

smolSetSurfaceAction = smoldyn.smolSetSurfaceAction
smolSetSurfaceAction.restype = c_int
smolSetSurfaceAction.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_int, c_int]

smolSetSurfaceRate = smoldyn.smolSetSurfaceRate
smolSetSurfaceRate.restype = c_int
smolSetSurfaceRate.argtypes = [c_void_p, c_char_p, c_char_p, c_int, c_int, c_int, c_double, c_char_p, c_int]

smolAddPanel = smoldyn.smolAddPanel
smolAddPanel.restype = c_int
smolAddPanel.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_char_p, c_void_p]

#extern "C" int            smolGetPanelIndex(simptr sim,const char *surface,enum PanelShape *panelshapeptr,const char *panel);
#smolDisplaySim = smoldyn.smolDisplaySim
#smolDisplaySim.restype = c_int
#smolDisplaySim.argtypes = [c_void_p, c_char_p, ???]
#
#extern "C" int            smolGetPanelIndexNT(simptr sim,const char *surface,enum PanelShape *panelshapeptr,const char *panel);
#smolDisplaySim = smoldyn.smolDisplaySim
#smolDisplaySim.restype = c_int
#smolDisplaySim.argtypes = [c_void_p]

smolGetPanelName = smoldyn.smolGetPanelName
smolGetPanelName.restype = c_char_p
smolGetPanelName.argtypes = [c_void_p, c_char_p, c_int, c_int, c_char_p]

smolSetPanelJump = smoldyn.smolSetPanelJump
smolSetPanelJump.restype = c_int
smolSetPanelJump.argtypes = [c_void_p, c_char_p, c_char_p, c_int, c_char_p, c_int, c_int]

smolAddSurfaceUnboundedEmitter = smoldyn.smolAddSurfaceUnboundedEmitter
smolAddSurfaceUnboundedEmitter.restype = c_int
smolAddSurfaceUnboundedEmitter.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_double, c_void_p]

smolSetSurfaceSimParams = smoldyn.smolSetSurfaceSimParams
smolSetSurfaceSimParams.restype = c_int
smolSetSurfaceSimParams.argtypes = [c_void_p, c_char_p, c_double]

smolAddPanelNeighbor = smoldyn.smolAddPanelNeighbor
smolAddPanelNeighbor.restype = c_int
smolAddPanelNeighbor.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p, c_char_p, c_int]

smolSetSurfaceStyle = smoldyn.smolSetSurfaceStyle
smolSetSurfaceStyle.restype = c_int
smolSetSurfaceStyle.argtypes = [c_void_p, c_char_p, c_int, c_int, c_double, c_void_p, c_int, c_int, c_double]


# Compartments
smolAddCompartment = smoldyn.smolAddCompartment
smolAddCompartment.restype = c_int
smolAddCompartment.argtypes = [c_void_p, c_char_p]

smolGetCompartmentIndex = smoldyn.smolGetCompartmentIndex
smolGetCompartmentIndex.restype = c_int
smolGetCompartmentIndex.argtypes = [c_void_p, c_char_p]

smolGetCompartmentIndexNT = smoldyn.smolGetCompartmentIndexNT
smolGetCompartmentIndexNT.restype = c_int
smolGetCompartmentIndexNT.argtypes = [c_void_p, c_char_p]

smolGetCompartmentName = smoldyn.smolGetCompartmentName
smolGetCompartmentName.restype = c_char_p
smolGetCompartmentName.argtypes = [c_void_p, c_int, c_char_p]

smolAddCompartmentSurface = smoldyn.smolAddCompartmentSurface
smolAddCompartmentSurface.restype = c_int
smolAddCompartmentSurface.argtypes = [c_void_p, c_char_p, c_char_p]

smolAddCompartmentPoint = smoldyn.smolAddCompartmentPoint
smolAddCompartmentPoint.restype = c_int
smolAddCompartmentPoint.argtypes = [c_void_p, c_char_p, c_void_p]

smolAddCompartmentLogic = smoldyn.smolAddCompartmentLogic
smolAddCompartmentLogic.restype = c_int
smolAddCompartmentLogic.argtypes = [c_void_p, c_char_p, c_int, c_char_p]

# Reactions 
smolAddReaction = smoldyn.smolAddReaction
smolAddReaction.restype = c_int
smolAddReaction.argtypes = [c_void_p, c_char_p, c_char_p, c_int, c_char_p, c_int, c_int, c_void_p, c_void_p, c_double]

smolGetReactionIndex = smoldyn.smolGetReactionIndex
smolGetReactionIndex.restype = c_int
smolGetReactionIndex.argtypes = [c_void_p, c_void_p, c_char_p]

smolGetReactionIndexNT = smoldyn.smolGetReactionIndexNT
smolGetReactionIndexNT.restype = c_int
smolGetReactionIndexNT.argtypes = [c_void_p, c_void_p, c_char_p]

smolGetReactionName = smoldyn.smolGetReactionName
smolGetReactionName.restype = c_char_p
smolGetReactionName.argtypes = [c_void_p, c_int, c_int, c_char_p]

smolSetReactionRate = smoldyn.smolSetReactionRate
smolSetReactionRate.restype = c_int
smolSetReactionRate.argtypes = [c_void_p, c_char_p, c_double, c_int]

smolSetReactionRegion = smoldyn.smolSetReactionRegion
smolSetReactionRegion.restype = c_int
smolSetReactionRegion.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p]

#enum ErrorCode smolSetReactionProducts(simptr sim,const char *reaction,enum RevParam method,double parameter,const char *product,double *position);
smolSetReactionProducts = smoldyn.smolSetReactionProducts
smolSetReactionProducts.restype = c_int
smolSetReactionProducts.argtypes = [c_void_p, c_char_p, c_int, c_double, c_char_p, c_void_p]

# Ports
#enum ErrorCode smolAddPort(simptr sim,const char *port,const char *surface,enum PanelFace face);
smolAddPort = smoldyn.smolAddPort
smolAddPort.restype = c_int
smolAddPort.argtypes = [c_void_p, c_char_p, c_char_p, c_int]

#int            smolGetPortIndex(simptr sim,const char *port);
smolGetPortIndex = smoldyn.smolGetPortIndex
smolGetPortIndex.restype = c_int
smolGetPortIndex.argtypes = [c_void_p, c_char_p]

#int            smolGetPortIndexNT(simptr sim,const char *port);
smolGetPortIndexNT = smoldyn.smolGetPortIndexNT
smolGetPortIndexNT.restype = c_int
smolGetPortIndexNT.argtypes = [c_void_p, c_char_p]

#char*          smolGetPortName(simptr sim,int portindex,char *port);
smolGetPortName = smoldyn.smolGetPortName
smolGetPortName.restype = c_char_p
smolGetPortName.argtypes = [c_void_p, c_int, c_char_p]  # Returns name in last arg, ugh!

#enum ErrorCode smolAddPortMolecules(simptr sim,const char *port,int nmolec,const char *species,double **positions);
smolAddPortMolecules = smoldyn.smolAddPortMolecules
smolAddPortMolecules.restypes = c_int
smolAddPortMolecules.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_void_p]

#int            smolGetPortMolecules(simptr sim,const char *port,const char *species,enum MolecState state,int remove);
smolGetPortMolecules = smoldyn.smolGetPortMolecules
smolGetPortMolecules.restypes = c_int
smolGetPortMolecules.argtypes = [c_void_p, c_char_p, c_char_p, c_int, c_int]
#


# Bit of test code
#smolGetVersion()
#lbound = (c_double * 2)(0.0, 0.0)
#ubound = (c_double * 2)(100.0, 100.0)
#s = smolNewSim(2, pointer(lbound), pointer(ubound))
#smolDisplaySim(s)
#smolFreeSim(s)