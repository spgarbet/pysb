from pysb.integrate import odesolve
from pylab import *
from pysb.smoldynlib import *
from pysb import bng
import numpy as np
import itertools

from pysb.examples.mixed.smoldyn import model,g as smoldyn,g
from pysb.examples.mixed.ode     import model   as ode

def smoldyn_mixed_simulation(smoldyn, ode, t, ports):

    # Get species names for ODE model
    bng.generate_equations(ode)
    
    # Generate a dict of repr(species)->species_index (used in interpreting port calls)
    ode_species_count = len(ode.species)
    smoldyn_monomer_count = len(smoldyn.model.monomers)
    ode_observable_count = len(ode.observables)
    ode_species_index = {}
    for species in ode.species:
        ode_species_index[repr(species).replace(" ", "")] = ode.get_species_index(species)
    
    # Set up initial values
    ode_initial = odesolve(ode, np.linspace(0,0,1))
    smolRunSimUntil(g.sim, 0.0)
    
    # Set up results array and load with initial values
    dtype = ode_initial.dtype.descr + zip([m.name for m in smoldyn.model.monomers], itertools.repeat('<f8'))
    ret_array = np.ndarray(len(t), dtype=dtype)
    for ode_species in range(ode_species_count+ode_observable_count):
        ret_array[dtype[ode_species][0]][0] = ode_initial[dtype[ode_species][0]][0]
    for smoldyn_monomer in [m.name for m in smoldyn.model.monomers]:
        ret_array[smoldyn_monomer][0] = smolGetMoleculeCount(g.sim, smoldyn_monomer, MolecState.ALL)
        
    prev_sec = 0
    
    for step in range(1,len(t)):
        sec = t[step]
        dt = sec - prev_sec
        y0 = [ret_array[step-1][i] for i in range(ode_species_count)]
        
        y = odesolve(ode, [0, dt], y0=y0)
        smolRunSimUntil(g.sim, sec)
        
        extra = [0] * ode_species_count
        
        for port in ports:
            for transfer in ports[port]:
                if not ode_species_index.has_key(transfer[1]):
                    raise Exception("Port ODE species must be concrete: %s" % transfer[1])
                extra[ode_species_index[transfer[1].replace(" ","")]] += \
                      smolGetPortMolecules(g.sim, port, transfer[0], MolecState.ALL, 1)
                
        for i in range(ode_species_count):
            ret_array[step][i] = y[1][i] + extra[i]
            
        for i in range(ode_observable_count):
            ret_array[step][i + ode_species_count] = y[1][i + ode_species_count]
            for spec in range(len(ode.observables[i].species)):
                ret_array[step][i] += extra[ode.observables[i].species[spec]] * ode.observables[i].coefficients[spec]
        
        for i in range(smoldyn_monomer_count):
            ret_array[step][i + ode_species_count + ode_observable_count] = smolGetMoleculeCount(g.sim, smoldyn.model.monomers[i].name, MolecState.SOLN)

        print sec
        prev_sec = sec
        
    return ret_array

ports = {"toODEfront": [("OP", "IE()")], \
         "toODEback":  [("OP", "IE()")]}
t = linspace(0, 0.25, 200*0.25+1)

result = smoldyn_mixed_simulation(smoldyn, ode, t, ports)

#print ""
#for name in [m.name for m in smoldyn.model.monomers]:
#    print name + " " + str(smolGetMoleculeCount(g.sim,name, MolecState.SOLN))
#
#smolRemoveMolecules(g.sim, "OE", MolecState.SOLN, 10)
#
#print ""
#for name in [m.name for m in smoldyn.model.monomers]:
#    print name + " " + str(smolGetMoleculeCount(g.sim,name, MolecState.SOLN))

plot(t, result['OS'])
plot(t, result['OE'])
plot(t, result['OC'])
plot(t, result['OP'])
plot(t, result['ISobs'])
plot(t, result['IEobs'])
plot(t, result['ICobs'])
plot(t, result['IPobs'])

#plot(t, result['OS'] +\
#        result['OE'] +\
#        2*result['OC'] +\
#        result['OP'] +\
#        result['ISobs'] +\
#        result['IEobs'] +\
#        2*result['ICobs'] +\
#        result['IPobs'])

show()
