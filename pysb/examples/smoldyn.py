from pysb import *
from pysb.geometry import *

Model()

Parameter('display_size', 3)
Parameter('time_start', 0)
Parameter('time_stop',  250)
Parameter('time_step',  0.01)

# Concentrations in number per cell
Parameter('EGF_0',      1.2e6)
Parameter('EGFR_0',     1.8e5)
Parameter('Grb2_0',     1.5e5)
Parameter('Sos_0',      6.2e4)

# Biomolecular rate constants are in (# per cell)^-1 s^-1,
#  obtained by dividing constants in M^-1 s^-1 by Na*V,
#  where Na is Avogadro's number and V is the volume
#  of the relevant compartment (the cytoplasm for all cases here).
# Unimolecular rate constants are in s^-1
Parameter('kp1',      1.667e-06) # ligand-monomer binding
Parameter('km1',           0.06) # ligand-monomer dissociation

Parameter('kp2',      5.556e-06) # aggregation of bound monomers
Parameter('km2',            0.1) # dissociation of bound monomers

Parameter('kp3',            0.5) # dimer transphosphorylation   
Parameter('km3',          4.505) # dimer dephosphorylation        

Parameter('kp4',      8.333e-07) # binding of Grb2 to receptor
Parameter('km4',           0.05) # dissociation of Grb2 from receptor

Parameter('kp5',      5.556e-06) # binding of Grb2 to Sos
Parameter('km5',           0.06) # dissociation of Grb2 from Sos

Parameter('kdeg',          0.01)

Monomer('ACA', {'Orient':['up', 'down']})  # Up, Down state (Down is active)
Monomer('ATP')
Monomer('cAMP')
Monomer('cAR1',{ 'Y': ['U','P']}) 

Parameter('aca_0', 30.0)
Parameter('car1_0', 30.0)
Initial(ACA(st="down") ** Membrane1, aca_0 )
Initial(cAR1(st="up")  ** Membrane1, car1_0)


m=Compartment('Main',      None,     geometry=SquareSpace(100, [0, 0]))
Compartment('Membrane1', Main,     geometry=SphericalSurface(10, [0,0]) )
Compartment('Cell1',     Membrane1, geometry=SphericalSpace(10, [0, 0]))

reaction_surface cellwall1  r2a  ACA(down) + ATP(bsoln) -> cAMP(bsoln) + ACA(down)  2.3

Rule("r2a", ACA(st="down")** Membrane1 + ATP(st="bsoln")** Membrane1 >> cAMP(st="bsoln")** Membrane1 + ACA(st="down")** Membrane1,  km6)
Rule("r2a", (ACA(st="down") + ATP(st="bsoln") >> cAMP(st="bsoln") + ACA(st="down"))**m ,  km6)

#####################################
# EXAMPLES FROM BNG EGFR MODEL
# 
# Ligand-receptor binding      
#EGFR(l,r) + EGF(r) <-> EGFR(l!1,r).EGF(r!1) kp1, km1
#Rule('ligand_receptor_binding',
#     EGFR(l=None, r=None) + EGF(r=None) <>
#     EGFR(l=1, r=None)    % EGF(r=1),
#     kp1, km1)
#
## Receptor-aggregation 
##EGFR(l!+,r) + EGFR(l!+,r) <-> EGFR(l!+,r!1).EGFR(l!+,r!1) kp2,km2
#Rule('receptor_aggregation',
#     EGFR(l=ANY, r=None) + EGFR(l=ANY, r=None) <>
#     EGFR(l=ANY, r=1)    % EGFR(l=ANY, r=1),
#     kp2, km2)
#
## Transphosphorylation of EGFR by RTK
##EGFR(r!+,Y1068~U) -> EGFR(r!+,Y1068~P)  kp3
#Rule('transphos_egfr',
#     EGFR(r=ANY, Y1068='U') >>
#     EGFR(r=ANY, Y1068='P'),
#     kp3)
#
## Dephosphorylayion
##EGFR(Y1068~P) -> EGFR(Y1068~U)  km3
#Rule('dephos_egfr',
#     EGFR(Y1068='P') >>
#     EGFR(Y1068='U'),
#     km3)
#
## Grb2 binding to pY1068
##EGFR(Y1068~P) + Grb2(SH2)   <-> EGFR(Y1068~P!1).Grb2(SH2!1)   kp4,km4
#Rule('grb2_bind_egfr',
#     EGFR(Y1068='P')     + Grb2(SH2=None) <>
#     EGFR(Y1068=('P',1)) % Grb2(SH2=1),
#     kp4, km4)
#
## Grb2 binding to Sos
##Grb2(SH2,SH3) + Sos(PR) <-> Grb2(SH2,SH3!1).Sos(PR!1) kp5,km5
#Rule('grb2_bind_sos',
#     Grb2(SH2=None, SH3=None) + Sos(PR=None) <>
#     Grb2(SH2=None, SH3=1)    % Sos(PR=1),
#     kp5, km5)

