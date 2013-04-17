import pysb
import ctypes
from pysb.smoldynlib import *

from pysb.geometry import *

class SmoldynlibGenerator(object):

    def __init__(self, model):
        self.sim = None
        self.model = model
    
    def __exit__(self, type, value, traceback):
        if(self.sim != None): smolFreeSim(self.sim)

    def generate_sim(self):
        self.open_sim()
        self.generate_parameters()
        self.generate_molecule_types()
        self.generate_compartments()

        self.generate_species()
        self.generate_reaction_rules()

    def run_sim(self):
        smolUpdateSim(self.sim)
        smolDisplaySim(self.sim)
        smolRunSim(self.sim)

    def open_sim(self):
        for c in self.model.compartments:
            # Master frame (universe boundaries)
            # FIXME: For now defaults to reflective, must be square hypersurface
            if c.parent is None:
                if not isinstance(c.geometry, SquareSpace):
                    raise Exception("SmolDyn specification must have a master SquareSpace to work within")
                lbound = (c_double * (c.dimension+1))()
                ubound = (c_double * (c.dimension+1))()
                for i in range(0, c.dimension):
                    lbound[i] = c.geometry.location[i] - c.geometry.shape.side/2
                    ubound[i] = c.geometry.location[i] + c.geometry.shape.side/2
                self.sim = smolNewSim(2, lbound, ubound)
        if self.sim is None:
            raise Exception("No main geometry specified")

    def generate_parameters(self):
        for p in self.model.parameters:
            if p.name in ['molperbox', 'boxsize'] :
                smolSetPartitions(self.sim, p.name, p.value)

        if not 'time_start' in self.model.parameters.keys()  :
            raise Exception("No time_start specified")
        if not 'time_stop' in self.model.parameters.keys() :
            raise Exception("No time_stop specified")
        if not 'time_step' in self.model.parameters.keys() :
            raise Exception("No time_step specified")

        smolSetSimTimes(self.sim, self.model.parameters['time_start'].value,
                                  self.model.parameters['time_stop'].value,
                                  self.model.parameters['time_step'].value)

    def generate_compartments(self):
        if not self.model.compartments:
            return
        for c in self.model.compartments:
            # Master frame (universe boundaries)
            if c.parent is None:
                pass # Dealt with in open_sim()
            # Deal with surfaces
            elif isinstance(c.geometry, SphericalSurface):
                smolAddSurface(self.sim, c.name)
                for action in c.action:
                    smolSetSurfaceAction(self.sim, c.name, PanelFace.BOTH, action[0], action[1], action[2])
                params = (c_double * (len(c.geometry.location)+3))()
                for idx, coord in enumerate(c.geometry.location):
                    params[idx] = coord
                params[-1] = 20
                params[-2] = 20
                params[-3] = c.geometry.shape.radius
                smolAddPanel(self.sim, c.name, PanelShape.SPHERE, "", "", params)
            elif isinstance(c.geometry, SphericalSpace) or isinstance(c.geometry, SquareSpace):
                smolAddCompartment(self.sim, c.name)
                smolAddCompartmentSurface(self.sim, c.name, c.parent.name)
                point = (c_double * len(c.geometry.location))()
                for idx, coord in enumerate(c.geometry.location):
                    point[idx] = coord
                smolAddCompartmentPoint(self.sim, c.name, point)
            else:
                raise Exception(("Unimplemented Geometry for Compartment %s") % (c.name))

    def generate_molecule_types(self):
        for m in self.model.monomers:
            #print "smolAddSpecies(<addr>, %s, '')" % m.name
            smolAddSpecies(self.sim, m.name, "")
# FIXME: How to add diffusion as configurable
            smolSetSpeciesMobility(self.sim, m.name, MolecState.ALL, 1, 0, 0)

    def generate_species(self):
        if not self.model.initial_conditions:
            raise Exception("Smoldyn generator requires initial conditions.")
        for cp, param in self.model.initial_conditions:
            quantity = param.value
            if(len(cp.monomer_patterns) > 1): raise Exception("Smoldyn does not support bound monomer_patterns")
            species = cp.monomer_patterns[0].monomer.name
            if(len(cp.monomer_patterns[0].site_conditions) > 1): raise Exception("Only one state supported in Smoldyn")
            
            if(len(cp.monomer_patterns[0].site_conditions.keys()) > 0):
                site = cp.monomer_patterns[0].site_conditions.keys()[0]
                state = self.smoldyn_state(cp.monomer_patterns[0].site_conditions[site])
            else:
                state = MolecState.SOLN
            c = cp.compartment
            if cp.compartment is None:
                for mp in cp.monomer_patterns:
                    if mp.compartment is not None: c = mp.compartment
                if c is None: raise Exception("All species must be in a compartment in Smoldyn")
            if(c.parent is None):
                smolAddSolutionMolecules(self.sim, species, int(quantity), 0, 0) #?
            elif isinstance(c.geometry, SphericalSurface):
                smolAddSurfaceMolecules(self.sim, species, state, int(quantity), c.name, PanelShape.ALL, "all", 0)
            else:
                smolAddCompartmentMolecules(self.sim, species, int(quantity), c.name)

    def smoldyn_state(self, state):
        if type(state) == str:
            if state.lower() == "soln":
                state_code = MolecState.SOLN
            elif state.lower() == "front":
                state_code = MolecState.FRONT
            elif state.lower() == "back":
                state_code = MolecState.BACK
            elif state.lower() == "up":
                state_code = MolecState.UP
            elif state.lower() == "down":
                state_code = MolecState.DOWN
            elif state.lower() == "bsoln":
                state_code = MolecState.BSOLN
            elif state.lower() == "all":
                state_code = MolecState.ALL
            elif state.lower() == "none":
                state_code = MolecState.NONE
            elif state.lower() == "some":
                state_code = MolecState.SOME
            else:
                raise Exception("Smoldyn generator has encountered an unknown string in a species pattern site condition.")
        elif state == None:
            state_code = MolecState.NONE
        else:
            raise Exception("Smoldyn generator has encountered an unknown element in a species pattern site condition.")
        return state_code

    def generate_reaction_rules(self):

        if not self.model.rules:
            return

        for r in self.model.rules:
            reactants = self.format_reaction_pattern(r.reactant_pattern)
            products  = self.format_reaction_pattern(r.product_pattern)

            l = len(products['name'])

            #import code; code.interact(local=locals())
            names  = (c_char_p * l)()
            states = (c_int * l)()
            notnull = 0
            for i in range(l):
                names[i]  = products['name'][i]
                states[i] = products['state'][i]
                if(names[i] != ''): notnull = notnull + 1
            smolAddReaction(self.sim, r.name, 
                            reactants['name'][0], reactants['state'][0],
                            reactants['name'][1], reactants['state'][1],
                            notnull, names, states, r.rate_forward.value)

            # I'm confused!!!
#            if r.is_reversible:
#               
#            else:
#                smolSetReactionProducts(self.sim, r.name, RevParam.IRREV, 0.0, 0, 0)

            # Which is it???
            if r.is_reversible:
                if(len(products['name']) < 1 or len(products['name']) > 2):
                    raise Exception("Rule %s unsupported number of products" % r.name)
                l = len(reactants['name'])
                names  = (c_char_p * l)()
                states = (c_int * l)()
                notnull = 0
                for i in range(l):
                    names[i]  = reactants['name'][i]
                    states[i] = reactants['state'][i]
                    if(names[i] != ''): notnull = notnull + 1
                #import code; code.interact(local=locals())
                smolAddReaction(self.sim, r.name+"reverse",
                                products['name'][0], products['state'][0], products['name'][1], products['state'][1],
                                notnull, names, states,
                                r.rate_reverse.value)
                smolSetReactionProducts(self.sim, r.name+"reverse", RevParam.PGEMMAX, 0.2, (c_char_p)(), 0)
            else:
                smolSetReactionProducts(self.sim, r.name, RevParam.IRREV, 0.0, (c_char_p)(), 0)
                                    
            if r.compartment is not None:
                if isinstance(r.compartment.geometry, SphericalSurface):
                    smolSetReactionRegion(self.sim, r.name, "", r.compartment.name)
                else:
                    smolSetReactionRegion(self.sim, r.name, r.compartment.name, "")
    
    def format_reaction_pattern(self, rp):
        names  = [self.format_complex_names(cp)  for cp in rp.complex_patterns]
        states = [self.format_complex_states(cp) for cp in rp.complex_patterns]
        while len(names)<2:
            names.append("")
            states.append(MolecState.ALL)
        return {'name':names, 'state':states}

    def format_complex_names(self, cp):
        if(len(cp.monomer_patterns) > 1):
            raise Exception("Complex Monomer patterns not supported in SmolDyn")
        return ("" if cp.monomer_patterns[0].monomer.name == None else cp.monomer_patterns[0].monomer.name)

    def format_complex_states(self, cp):
        monomer = cp.monomer_patterns[0]
        if(len(monomer.site_conditions) > 1):
            raise Exception("Complex Monomer patterns not supported in SmolDyn")
        if(len(monomer.site_conditions) < 1):
            if(monomer.monomer.name == None or monomer.monomer.name == ''):
                return MolecState.ALL
            else:
                return MolecState.SOLN
        else:
            site = monomer.site_conditions.keys()[0]
            return self.smoldyn_state(monomer.site_conditions[site])