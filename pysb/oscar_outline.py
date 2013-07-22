# This is "stubbed" outline
# The return types in terms of arrays / dictionaries are well defined
# Don't change this file unless you've though *carefully* about the impact, and
# the information flow through the system. 
# Don't write over this file, keep it as a reference as you create the real functions.
# Using any of these functions, you can skip around in tasks orders, as they 
# return the expected results for the tyson model.
# Use the output of each to compare with the functions you create

from sympy.solvers import solve
from sympy import Symbol
from sympy import symbols
from sympy import symarray
from sympy.functions.elementary.complexes import Abs
from sympy import solve_poly_system
from sympy import log
from sympy.functions.special.delta_functions import Heaviside
from sympy import simplify
from sympy import Mul
from sympy import log
import pysb
import pysb.bng
import sympy
import re
import sys
import os
import pygraphviz
import networkx
import copy
from sympy.parsing.sympy_parser import parse_expr
from collections import Mapping






def find_slaves(model, t, ignore=15, epsilon=1e-6):
    return ['s0', 's1', 's5']
    slaves = []

    generate_equations(model)
    x = odesolve(model, t)
    x = x[ignore:] # Ignore first couple points
    t = t[ignore:]
    names = [n for n in filter(lambda n: n.startswith('__'), x.dtype.names)]
    x = x[names] # Only concrete species are considered
    names = [n.replace('__','') for n in names]
    x.dtype = [(n,'<f8') for n in names]

    for i, eq in enumerate(model.odes): # i is equation number
        eq   = eq.subs('s%d' % i, 's%dstar' % i)
        sol  = solve(eq, Symbol('s%dstar' % i)) # Find equation of imposed trace
        for j in range(len(sol)):  # j is solution j for equation i
            prueba = zeros(len(x))
            for p in model.parameters: sol[j] = sol[j].subs(p.name, p.value) # Substitute parameters
            # This is the loop that kills performance
            #prueba = Abs(sol[j].evalf(subs=x) - x)
            for l, tt in enumerate(t):
                prueba[l] = Abs(sol[j].evalf(subs={n:x[tt][n] for n in names}) - x[tt]['s%d'%i])
            if (prueba.max() <= epsilon): slaves.append("s%d" % i)
            print prueba.max()
    return slaves

# The output type may change, as needed for a graph package
# Large time (interacting with BNG)

# This is a function which builds the edges according to the nodes
def r_link(graph, s, r, **attrs):
    nodes = ('s%d' % s, 's%d' % r)
    if attrs.get('_flip'):
        del attrs['_flip']
        nodes = reversed(nodes)
    attrs.setdefault('arrowhead', 'normal')
    graph.add_edge(*nodes, **attrs)

def find_cycles(model):
    """
    Render the reactions produced by a model into the "dot" graph format.

    Parameters
    ----------
    model : pysb.core.Model
        The model to render.

    Returns
    -------
    sorted graph edges
    """

    pysb.bng.generate_equations(model)

    graph = networkx.DiGraph(rankdir="LR")
    ic_species = [cp for cp, parameter in model.initial_conditions]
    for i, cp in enumerate(model.species):
        species_node = 's%d' % i
        slabel = re.sub(r'% ', r'%\\l', str(cp))
        slabel += '\\l'
        color = "#ccffcc"
        # color species with an initial condition differently
        if len([s for s in ic_species if s.is_equivalent_to(cp)]):
            color = "#aaffff"
        graph.add_node(species_node,
                       label=species_node,
                       shape="Mrecord",
                       fillcolor=color, style="filled", color="transparent",
                       fontsize="12",
                       margin="0.06,0")
    for i, reaction in enumerate(model.reactions):       
        reactants = set(reaction['reactants'])
        products = set(reaction['products']) 
        attr_reversible = {}
        for s in reactants:
            for p in products:
                r_link(graph, s, p, **attr_reversible)
 #   networkx.draw(graph) 
 #   plt.show() 
    return networkx.simple_cycles(graph) #graph.edges() returns the edges

#This function finds conservation laws from the conserved cycles
def mass_conserved(model):
    c = find_cycles(model)
    h = []
    g = []
    for i, item in enumerate(c):
        b = 0
        u = 0
        del item[len(item)-1]
        for j, specie in enumerate(item):
            b += model.odes[int(re.findall(r'\d+', c[i][j])[0])]
        if b == 0:
            g.append(item)
            for l,k in enumerate(item):
                u += sympy.Symbol(c[i][l])    
            h.append(u-sympy.Symbol('C%d'%i))
            print 'cycle%d'%i, 'is conserved'
            
    return h, g

# Might need a "Prune" equation function

# Large time sink, tropicalization step is needed in here, i.e. maximum
def slave_equations(model, t, ignore=15, epsilon=1e-6):
    slaves = find_slaves(model, t, ignore=15, epsilon=1e-6)
    slave_conserved_eqs = {}
    for i, j in enumerate(slaves):
        slave_conserved_eqs[j] = model.odes[int(re.findall(r'\d+', slaves[i])[0])] 

#        slave_conserved_eqs.setdefault(j,[]).append(model.odes[int(re.findall(r'\d+', slaves[i])[0])])
        
    # Solve the slave equations here
    # Stubbed computation
#    slave_eq = {
#                 's0': (Symbol('C1')-Symbol('s5')-Symbol('s6'))*Symbol('k9')/(Symbol('k8')+Symbol('k9')),
#                 's1': Symbol('k1')*(Symbol('k8')+Symbol('k9'))/(Symbol('k3')*Symbol('k8')*(Symbol('C1')-Symbol('s5')-Symbol('s6'))),
#                 's4': (Symbol('C1')-Symbol('s5')-Symbol('s6'))*Symbol('k8')/(Symbol('k8')+Symbol('k9'))
#               } # Stub
    return slave_conserved_eqs

def pruned_equations(model, t, ignore=15, epsilon=1e-6, rho=1):

    k8, s5, k9, s0, k3, s1, k1, s2, k2, k3, k6, s6 = symbols('k8 s5 k9 s0 k3 s1 k1 s2 k2 k3 k6 s6')
    a = [k8*s5-k9*s0, -k8*s5+k9*s0, k1*s2-k2*s1-k3*s0*s1]
    return a

    generate_equations(model)
    x = odesolve(model, t)
    names = [n for n in filter(lambda n: n.startswith('__'), x.dtype.names)]
    x = x[names] # Only concrete species are considered
    names = [n.replace('__','') for n in names]
    x.dtype = [(n,'<f8') for n in names]
    conservation = mass_conserved(model)[0]

    pruned_eqs = slave_equations(model, t, ignore=15, epsilon=1e-6)
    eq = copy.deepcopy(pruned_eqs)
    for i, j in enumerate(eq):
        ble = re.findall(r'\b\S+\b', str(j)) #Creates a list of the monomials of each slave equation
        for l, m in enumerate(ble): #Compares the monomials to find the pruned system
            m_ready=parse_expr(m)
            m_elim=parse_expr(m)
            for p in model.parameters: m_ready = m_ready.subs(p.name, p.value) # Substitute parameters
            elim=zeros(len(x))
            for k in range(len(ble)):
                if (k+l+1) <= (len(ble)-1):
                    ble_ready = parse_expr(ble[k+l+1])
                    ble_elim = parse_expr(ble[k+l+1])
                    for p in model.parameters: ble_ready = ble_ready.subs(p.name, p.value) # Substitute parameters
                    for s, tt in enumerate(t):
                        elim[s] = m_ready.evalf(subs={n:x[tt][n] for n in names}) - ble_ready.evalf(subs={n:x[tt][n] for n in names})  
                    if elim.max() > 0 and abs(elim.max()) > rho:
                       pruned_eqs[i] = pruned_eqs[i].subs(ble_elim, 0)
                    else: pass
                    if elim.max() < 0 and abs(elim.max()) > rho:
                       pruned_eqs[i] = pruned_eqs[i].subs(m_elim, 0) 
                    else: pass
                else: pass 
    for i, l in enumerate(conservation): #Add the conservation laws to the pruned system
        pruned_eqs.append(l)
    return pruned_eqs
 

def diff_alg_system(model):
    sol_dict = {}
    index_slaves = []
    slaves = find_slaves(model, t, ignore=15, epsilon=1e-6)
    var_ready = []
    eqs_to_add = copy.deepcopy(model.odes)
    eqs_to_add_dict = {}
    var = find_slaves(model, t, ignore=15, epsilon=1e-6)
    eqs = pruned_equations(model, t, ignore=15, epsilon=1e-6, rho=1)
    w = mass_conserved(model)[1]
    cycle_eqs = mass_conserved(model)[0]
    for i in cycle_eqs:
        eqs.append(i)
    for i in w: #Adds the variable of s2 cycle, it is required because the solver doesnt know if s2 or C2 is theconstant 
        if len(i) == 1:
            var.append(i[0])   
    for j in var:
        var_ready.append(Symbol(j))
    sol = solve_poly_system(eqs, var_ready)
    for i, j in enumerate(var_ready):
        sol_dict[j] = sol[0][i] 

    for i, j in enumerate(eqs_to_add):
        eqs_to_add_dict[Symbol('s%d'%i)] = j
    for i in slaves:
        del eqs_to_add_dict[Symbol('%s'%i)]
        
    eqs_to_add_ready = copy.deepcopy(eqs_to_add_dict)    

    for i in eqs_to_add_dict: #Changes s2 to (d/dt)s2
        eqs_to_add_ready[Symbol('(d/dt)%s'%i)] = eqs_to_add_ready.pop(i)
    for l in eqs_to_add_ready.keys(): #Substitutes the values of the algebraic system
        for i, j in enumerate(sol_dict): eqs_to_add_ready[l]=eqs_to_add_ready[l].subs(sol_dict.keys()[i], sol_dict.values()[i])
    
    return eqs_to_add_ready


def tropicalization(model):

    eqs_for_tropicalization = diff_alg_system(model) 
    tropicalized = {}

    for i in eqs_for_tropicalization.keys():
        for par in model.parameters: eqs_for_tropicalization[i] = simplify(eqs_for_tropicalization[i].subs(par.name, par.value)) # Substitute parameters 

    for j in sorted(eqs_for_tropicalization.keys()):
        if type(eqs_for_tropicalization[j]) == Mul: print solve(log(j), dict = True) #If Mul=True there is only one monomial
        elif eqs_for_tropicalization[j] == 0: print 'there are not monomials'
        else:            
            ar = eqs_for_tropicalization[j].args #List of the terms of each equation  
            asd=0 
            for l, k in enumerate(ar):
                p = k
                for f, h in enumerate(ar):
                   if k != h:
                      p *= Heaviside(log(abs(k)) - log(abs(h)))
                asd +=p
            tropicalized[j] = asd
    return tropicalized    

