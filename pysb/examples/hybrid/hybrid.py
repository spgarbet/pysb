from pysb import *
from pysb.generator import *
from pysb.geometry import *
from pysb.integrate import odesolve
from pysb.smoldynlib import *
from pysb.generator.smoldynlib import *
from pysb import bng
import numpy
import itertools

# ---------------------------------------------------------------------------
# H E L P E R   F U N C T I O N S
# ---------------------------------------------------------------------------

# Returns True iff a ComplexPattern in contained within a listlike object (using cp.is_equivalent_to())
# sil == species_in_list
def sil(cp,l):
    for i in l:
        if cp.is_equivalent_to(i): return True
    return False

# Finds the index of a ComplexPattern within a listlike object.  If list does not contain cp, returns -1
def get_cp_index(cp,l):
    for i in range(len(l)):
        if cp.is_equivalent_to(l[i]): return i
    return -1

# Temporary function for error logging
def flogger(flog,title,l):
    flog.write(title+"\n")
    for i in l:
        flog.write(repr(i) + "\n")
    flog.write("\n")


# Takes a list of ComplexPatterns and formats them like BNG (e.g. 'S(b=1)%E(b=1)' -> 'E(b=1) % S(b=1)')
def alphabetize(list_of_cp_strings):
    ret = []
    for cp_str in list_of_cp_strings:
        cp_str = cp_str.replace(' ','')
        monomer_patterns = cp_str.split('%')
        monomer_patterns.sort()
        cp = monomer_patterns[0]
        for i in monomer_patterns[1:]:
            cp += ' %% %s' %(i)
        ret.append(cp)
    return ret

# ---------------------------------------------------------------------------
# R U N   S I M U L A T I O N
# ---------------------------------------------------------------------------

def run_hybrid(base_model, t, smoldyn_species):
    flog = open('/home/dwooten/logs/pysb_log' + str(numpy.random.randint(10000)),'w')
    
    bng.generate_equations(base_model)
    
    smol_species_str = alphabetize([i[0] for i in smoldyn_species])
    ode_species_str = [repr(i) for i in base_model.species if not repr(i) in smol_species_str]

    ode_species = [spec for spec in base_model.species if repr(spec) in ode_species_str]
    smol_species = [spec for spec in base_model.species if repr(spec) in smol_species_str]



    # ---------------------------------------
    # C R E A T E   N E W   O D E   M O D E L
    # ---------------------------------------
    
    Model(name="ode")
    # Copy ALL monomers to ODE
    # FIXME: In the future, just copy ODE monomers over
    for m in base_model.monomers:
        Monomer(m.name, sites=m.sites, site_states=m.site_states, difc=m.difc)
            
    # Copy reaction rates and rules to ODE (only reactions where ALL reactants are in ODEs)
    for rxn in base_model.reactions:
        reactants = [base_model.species[i] for i in rxn['reactants']]
        if not False in [i in ode_species for i in reactants]:  #If ALL reactants are listed in ode_species (if some are in Smoldyn, skip for now)
            
            # Get the reaction rule from the base model
            rule = base_model.rules[rxn['rule']]        # The associated rule
            
            if not rxn['reverse']:  # If it's a forward reaction
                rate_p = rule.rate_forward
                reactant_pattern = rule.reactant_pattern
                product_pattern = rule.product_pattern
                name = rule.name
            else:               # Reverse reaction
                rate_p = rule.rate_reverse
                reactant_pattern = rule.product_pattern
                product_pattern = rule.reactant_pattern
                name = rule.name + "reverse"
            rate_p = Parameter(rate_p.name,rate_p.value)
            Rule(name,reactant_pattern >> product_pattern,rate_forward=rate_p)
            
    # Determine which things "need" initial conditions
    # FIXME: This is definitely the wrong approach to take longterm.
    # The problem we are addressing is that BNG will not generate equations that can never happen,
    # such as C >> P + E, if C starts wit C0=0.  This is usually fine because you will also have
    # S + E >> C, and you will start with S and E, but in this case those may be in the ODE half,
    # and invisible here.
    need = []
    for rule in ode.rules:
        for cp in rule.reactant_pattern.complex_patterns:
            if not sil(cp,need): need.append(cp)
            
    # Copy ODE-SPECIFIC initial conditions to ODE
    for init in base_model.initial_conditions:
        if sil(init[0],ode_species):
            param=Parameter(init[1].name,init[1].value)
            Initial(init[0],param)
            
            got=need[get_cp_index(init[0],need)]
            need.remove(got)
    
    needcount=0
    for cp in need:
        if (ode.get_species_index(cp) == None):
            print "fmylife"
            print cp
        param=Parameter("autoinit"+str(needcount),0.00001)
        Initial(cp,param)
        needcount += 1
    
    
    bng.generate_equations(ode)
    
    
    # -----------------------------------------------
    # C R E A T E   N E W   S M O L D Y N   M O D E L
    # -----------------------------------------------
    
    Model(name="smoldyn")
    # Copy ALL monomers to Smoldyn
    # FIXME: In the future, just copy Smoldyn monomers (maybe?)
    for m in base_model.monomers:
        Monomer(m.name, sites=m.sites, site_states=m.site_states, difc=m.difc)
    
    main=Compartment('Main', None, geometry=SquareSpace(1, [0, 0, 0]))
    
    # Copy reaction rates and rules to Smoldyn
    for rxn in base_model.reactions:
        reactants = [base_model.species[i] for i in rxn['reactants']]
        if not False in [i in smol_species for i in reactants]:  #If ALL reactants are listed in smol_species
            r = repr(rxn['rate']).split("*")[0]
            rate_p = base_model.parameters.get(r)
            rate_p = Parameter(rate_p.name,rate_p.value)
            rule = base_model.rules[rxn['rule']]
            if not rxn['reverse']:
                reactant_pattern = rule.reactant_pattern
                product_pattern = rule.product_pattern
            else:
                reactant_pattern = rule.product_pattern
                product_pattern = rule.reactant_pattern
            Rule(rule.name,reactant_pattern >> product_pattern,rate_forward=rate_p)
    
    # Determine which things "need" initial conditions
    # FIXME: This is definitely the wrong approach to take longterm.
    # The problem we are addressing is that BNG will not generate equations that can never happen,
    # such as C >> P + E, if C starts wit C0=0.  This is usually fine because you will also have
    # S + E >> C, and you will start with S and E, but in this case those may be in the ODE half,
    # and invisible here.
    need = []
    for rule in smoldyn.rules:
        for cp in rule.reactant_pattern.complex_patterns:
            if not sil(cp,need): need.append(cp)
        
    # Copy initial conditions to Smoldyn
    for init in base_model.initial_conditions:
        if sil(init[0],smol_species):
            param=Parameter(init[1].name,init[1].value)
            Initial(init[0]**main,param)
            # We have added this parameter, so remove it from the need list
            got=need[get_cp_index(init[0],need)]
            need.remove(got)
    
    for cp in need:
        param=Parameter(repr(cp).split('(')[0]+"autoinit",0.01)
        Initial(cp**main,param)
    
    
    # Setup Smoldyn time parameters
    Parameter('time_start', t[0])
    Parameter('time_stop',  t[-1]+t[1]-t[0])
    Parameter('time_step',  t[1]-t[0])
    
    # Generate Smoldyn species, equations
    bng.generate_equations(smoldyn)
    
    # Setup Smoldyn compartment (set volume = 1, molecules reflect)
    membrane_reflect = [[repr(m), MolecState.ALL, SurfAction.REFLECT] for m in smoldyn.species]
    main.action = membrane_reflect
    
    # Generate Smoldyn sim
    g = SmoldynlibGenerator(smoldyn)
    g.generate_sim()
    
    
    lookup = {}
    for species in base_model.species:
        lookup[repr(species)]={'base':base_model.get_species_index(species)}
        
    for species in ode.species:
        lookup[repr(species)]['ode']=ode.get_species_index(species)
        
    for species in smoldyn.species:
        lookup[g.format_species_name(species)]['smoldyn']=smoldyn.get_species_index(species)
    
    # -------------------------------------------------------------------------
    # R U N   S I M U L A T I O N
    # -------------------------------------------------------------------------
    
    # Set up initials
    result = [list(odesolve(ode, numpy.linspace(0,0,1))[0])]    # This only gets the reactants and products relevant to ODE, add Smoldyn results as well
    result[0] += [0.0]*(len(base_model.species) - len(ode.species))   #Right now initialize all smoldyn results to 0 (because it works, for now)
    smolUpdateSim(g.sim)
    
    ode_species_count = len(ode.species)
    
    # Load initial vector for ODE solver
    ode_init = [result[0][ode.get_species_index(i)] for i in ode.species]
    
    prev_sec = 0
    for sec in t[1:]:
        
        # FIXME: This throws ALL older ODE solution information away, other than the previous timestep
        print "Running ODE t=%0.3f" %(sec)
        ode_step_z = odesolve(ode,[prev_sec,sec],y0=ode_init)
        ode_step = list(ode_step_z[1]) + [0]
        
        print "Running Smoldyn t=%0.3f" %(sec)
        smolRunSimUntil(g.sim,sec)
        
        
        # Move ODE products to Smoldyn
        print "Moving ODE products to Smoldyn"
        for species in smol_species_str:

            # If this is not produced by an ODE reaction, just continue
            if not lookup[species].has_key('ode'): continue
            
            # Find the index of this species in the ODE model
            index = lookup[species]['ode']
            
            # Add those molecules to Smoldyn
            # Here just add the integer part of #molecules to Smoldyn (can't add fractional molecules)
            # Then leave ode_step[i] as the fractional part.  That way, if it takes multiple timesteps for ODE to
            # generate 1 molecule, it will do so.
            smolAddSolutionMolecules(g.sim,species,int(ode_step[index]),None,None)
            
            # Remove them from ODE
            ode_step[index] -= int(ode_step[index])
            
        # Move Smoldyn products to ODE
        print "Moving Smoldyn products to ODE"
        for species in ode_species_str:
           
            #If this is not in BOTH Smoldyn AND ODE, just continue
            if not lookup[species].has_key('smoldyn'): continue
            if not lookup[species].has_key('ode'): continue
            
            # Get the Smoldyn species name (string representation of the ComplexPattern without the compartment)
            smol_index = lookup[species]['smoldyn']
            smoldyn_name = g.format_species_name(smoldyn.species[smol_index])

            # Find where to put it in ODE
            ode_index = lookup[species]['ode']
            
            # Add them to ODE
            count = smolGetMoleculeCount(g.sim,smoldyn_name,MolecState.SOLN)
            ode_step[ode_index] += count
            
            # Remove them from Smoldyn
            smolRemoveMolecules(g.sim,smoldyn_name,MolecState.SOLN,count)
    
        # Aggregate all results from this step into a single list    
        this_step = [0]*len(base_model.species)
        for species in smol_species_str:
            index = lookup[species]['base']
            this_step[index] = smolGetMoleculeCount(g.sim,species,MolecState.SOLN)
        for species in ode_species_str:
            # FIXME: This is tricky.  If we have a system where a species is produced in Smoldyn,
            # but never actually used as a reactant, we could leave it in Smoldyn (which would waste CPU
            # cycles), move to ODE (which would cost nothing, but be awkward), or set aside outside of
            # the models.  The frustrating part is that because you probably won't specify it in
            # "smoldyn_species" (the list of strings), my code will ASSUME it belongs in "ode_species",
            # even though it doesn't have an ODE index.
            #
            # Temporary fix - make sure that even though it is listed in ode_species, that it has
            # an ODE index.  Otherwise, grab it from Smoldyn here, not ODE
            # UPDATE: right now this never happens, but I'm not 100% I fixed the root problem
            if not lookup[species].has_key('ode'):
                index = lookup[species]['base']
                print "Tricky %s" %(species)
                this_step[index] = smolGetMoleculeCount(g.sim, species, MolecState.SOLN)
            else:
                base_index = lookup[species]['base']
                ode_index = lookup[species]['ode']
                this_step[base_index] = ode_step[ode_index]
        
        
        result.append(this_step)
        
        ode_init = [ode_step[ode.get_species_index(i)] for i in ode.species]
        
        prev_sec = sec
    
    for obs in base_model.observables:
        for res in result:
            concentrations = [res[i] for i in obs.species]
            res.append(sum([concentrations[i] * obs.coefficients[i] for i in range(len(concentrations))]))

    dtype = zip(['_s%d'%(i) for i in range(len(base_model.species))], itertools.repeat('<f8'))
    dtype += zip([obs.name for obs in base_model.observables], itertools.repeat('<f8'))
    
    ret_array = numpy.ndarray(len(t), dtype=dtype)
    for sec in range(len(t)):
        for i in range(len(dtype)):
            ret_array[sec][i] = result[sec][i]
    return ret_array
