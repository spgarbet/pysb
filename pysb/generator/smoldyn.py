import pysb

from pysb.geometry import *

class SmoldynGenerator(object):

    def __init__(self, model):
        self.model = model
        self.__content = None

    def get_content(self):
        if self.__content == None:
            self.generate_content()
        return self.__content

    def generate_content(self):
        self.__content = 'graphics opengl\n'
        self.generate_molecule_types()
        self.generate_parameters()
        self.generate_compartments()
        self.generate_species()
        self.generate_reaction_rules()
        self.__content += 'end_file'

    def generate_parameters(self):
        for p in self.model.parameters:
            if p.name in ['time_start', 'time_stop', 'time_step'] :
                self.__content += ("%s %f\n") % (p.name, p.value)
            if p.name in ['display_size'] :
                self.__content += ("%s all %d\n") % (p.name, p.value)
        self.__content += "\n"

    def generate_compartments(self):
        if not self.model.compartments:
            return
        for c in self.model.compartments:
            # Master frame (universe boundaries)
            # FIXME: For now defaults to reflective, must be square hypersurface
            # FIXME: Also thickness and colors are defaults. Ugh
            if c.parent is None:
                if not isinstance(c.geometry, SquareSpace):
                    raise Exception("SmolDyn specification must have a master SquareSpace to work within")
                self.__content += ("dim %d\n\n") % (c.dimension)
                for i in range(0, c.dimension):
                    l = c.geometry.location[i] - c.geometry.shape.side/2
                    r = c.geometry.location[i] + c.geometry.shape.side/2
                    self.__content += ("boundaries %d %d %d\n") % (i,l,r)
                self.__content += "frame_thickness 1\n\n"
            # Deal with surfaces
            elif isinstance(c.geometry, SphericalSurface):
                self.__content += ("start_surface %s\n") % (c.name)
                if(c.action is not None):
                    self.__content += "action "
                    self.__content += "\naction ".join(c.action)
                    self.__content += "\n"
                self.__content += "panel sph"
                for coord in c.geometry.location:
                    self.__content += (" %d") % (coord)
                self.__content += (" %d 20 20\n") % (c.geometry.shape.radius)
                # FIXME: How to deal with thickness / colors
                self.__content += "thickness 1\ncolor both grey 0.5\n"
                self.__content += "end_surface\n"
            elif isinstance(c.geometry, SphericalSpace):
                self.__content += ("start_compartment %s\n") % (c.name)
                self.__content += ("surface %s\n") % (c.parent.name)
                self.__content += "point "
                for coord in c.geometry.location:
                    self.__content += (" %d") % (coord)
                self.__content += "\nend_compartment\n"
            elif isinstance(c.geometry, SquareSpace):
                self.__content += ("start_compartment %s\n") % (c.name)
                self.__content += ("surface %s\n") % (c.parent.name)
                self.__content += "point"
                for coord in c.geometry.location:
                    self.__content += (" %d") % (coord)
                self.__content += "\nend_compartment\n"
            else:
                raise Exception(("Unimplemented Geometry for Compartment %s") % (c.name))
            self.__content += "\n"

    def generate_molecule_types(self):
        colors = ['maroon', 'red','orange','yellow','olive','green','purple','fuchsia','lime','teal','aqua','blue','navy','black','grey','silver']
        
        self.__content += "species "
        for m in self.model.monomers:
            self.__content += "  %s" % (m.name)
        self.__content += "\n\n"
        for m in self.model.monomers:
            self.__content += "difc %s(all) %lf\n" % (m.name, m.difc)
        self.__content += "\n\n"
        for i,m in enumerate(self.model.monomers):
            self.__content += "color %s(all) %s\n" % (m.name, colors[i % len(colors)])
        self.__content += "\n\n"

    def generate_species(self):
        if not self.model.initial_conditions:
            raise Exception("Smoldyn generator requires initial conditions.")
        for cp, param in self.model.initial_conditions:
            self.__content += self.format_speciespattern(cp, param.value)
            self.__content += "\n"
        self.__content += "\n"
    
    def format_speciespattern(self, cp, quantity):
        ret = '.'.join([self.format_speciessites(mp) for mp in cp.monomer_patterns])
        c = cp.compartment
        if cp.compartment is None:
            for mp in cp.monomer_patterns:
                if mp.compartment is not None: c = mp.compartment
            if c is None: raise Exception("All species must be in a compartment in Smoldyn")
        if isinstance(c.geometry, SphericalSurface):
            ret = 'surface_mol %d %s %s all all' % (quantity, ret, c.name)
        else:
            #import code
            #code.interact(local=locals())

            if c.parent is None:
                ret = ('mol %d %s ' % (quantity, ret)) + ' '.join(["%s" % el for el in c.geometry.location])
            else:
                ret = 'compartment_mol %d %s %s' % (quantity, ret, c.name)
        return ret

    def format_speciessites(self, mp):
        # sort sites in the same order given in the original Monomer
        site_conditions = sorted(mp.site_conditions.items(),
                                 key=lambda x: mp.monomer.sites.index(x[0]))
        site_pattern_code = ','.join([self.format_site_species(site, state) for (site, state) in site_conditions])
        ret = '%s(%s)' % (mp.monomer.name, site_pattern_code)
        return ret

    def format_site_species(self, site, state):
        if type(state) == str:
            state_code = state
        elif state == None:
            state_code = ''
        else:
            raise Exception("Smoldyn generator has encountered an unknown element in a species pattern site condition.")
        return '%s' % (state_code)

    def generate_reaction_rules(self):

        if not self.model.rules:
            return

        #if cp.compartment is None:
        #    for mp in cp.monomer_patterns:
        #        if mp.compartment is not None: c = mp.compartment
        #    if c is None: raise Exception("All species must be in a compartment in Smoldyn")
        #if isinstance(c.geometry, SphericalSurface):
        #    ret = 'surface_mol %d %s %s all all' % (quantity, ret, c.name)

        for r in self.model.rules:
            reaction = "reaction"

            if r.compartment is not None:
                if isinstance(r.compartment.geometry, SphericalSurface):
                    reaction = "reaction_surface %s" % r.compartment.name
                else:
                    reaction = "reaction_cmpt %s" % r.compartment.name

            label   = r.name
            react_p = r.reactant_pattern
            prod_p  = r.product_pattern
            reactants_code = self.format_reactionpattern(react_p)
            products_code  = self.format_reactionpattern(prod_p)
            if r.is_reversible:
                arrow = '<->'
            else:
                arrow = '->'

            if len(products_code)==0:
                products_code = '0'

            if len(reactants_code)==0:
                reactants_code = '0'

            self.__content += ("%s %s %s %s %s %s\n") % \
                (reaction, label, reactants_code, arrow, products_code, r.rate_forward.value)
        self.__content += "\n"

    def format_complexpatternrule(self, cp, default=''):
        return '.'.join([self.format_monomerpattern(mp,default) for mp in cp.monomer_patterns])

    def format_monomer_site(self, monomer, site):
        ret = site
        if monomer.site_states.has_key(site):
            for state in monomer.site_states[site]:
                ret += '~' + state
        return ret

    def format_reactionpattern(self, rp):
        #import code
        #code.interact(local=locals())
        site_count = sum(sum(map(lambda x: map(lambda y: len(y.site_conditions),x.monomer_patterns), rp.complex_patterns), []))
        if site_count > 0: 
            default = 'bsoln'
        else:
            default = ''
        return ' + '.join([self.format_complexpatternrule(cp,default) for cp in rp.complex_patterns])

    def format_monomerpattern(self, mp, default=''):
        # sort sites in the same order given in the original Monomer
        site_conditions = sorted(mp.site_conditions.items(),
                                 key=lambda x: mp.monomer.sites.index(x[0]))
        site_pattern_code = ','.join([self.format_site_condition(site, state) for (site, state) in site_conditions])
        ret = mp.monomer.name
        if len(site_pattern_code)!=0:
            ret += '(%s)' % site_pattern_code
        elif len(default) > 0:
            ret += '(%s)' % default
        #if mp.compartment is not None:
        #    ret = '%s@%s' % (ret, mp.compartment.name)
        return ret

    def format_site_condition(self, site, state):
        # empty
        if state == None:
            state_code = ''
        # single bond
        #elif type(state) == int:
        #    state_code = '!' + str(state)
        # multiple bonds
        #elif type(state) == list and all(isinstance(s, int) for s in state):
        #    state_code = ''.join('!%d' % s for s in state)
        ## state
        elif type(state) == str:
            state_code = state
        # state AND single bond
        #elif type(state) == tuple:
        #    # bond is wildcard (zero or more unspecified bonds)
        #    if state[1] == pysb.WILD:
        #        state = (state[0], '?')
        #    state_code = '~%s!%s' % state
        # one or more unspecified bonds
        #elif state == pysb.ANY:
        #    state_code = '!+'
        else:
            raise Exception("Smoldyn generator has encountered an unknown element in a rule pattern site condition.")
        #return '%s%s' % (site, state_code)
        return '%s' % (state_code)
    
