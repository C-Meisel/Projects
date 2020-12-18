#Initilize variables for the PCEC model
#as of 30Nov20 most of these variables are hard coded and need to be found
#The solution vector is initialized here
#Everything here is saved into a pointer class

#\/\/^\/^\/\/^\ Imports /\/\/^\/^\/\/^\#
import numpy as np
import math
from Charge_Transfer_Init import *

#\/\/^\/^\/\/^\ Parameters /\/\/^\/^\/\/^\#
"Importnat Variables For the SV to work"
i_ext = 40 #A
C_dl_neg = 6e5 # F/m2 this makes it so my function does not go to negative infinity
C_dl_pos = 2e2 # F/m2 (need to look up)
t_final = 10000 #seconds

"Physical Constants:"
F = 96485 #C/mol e-
R = 8.314 #J/mol*K

"Equation values"
T = 823 #K

#Node thicknesses:
dy1_neg = 980e-6 # m
dy2_neg = 20e-6 # m
dy_ely = 10e-6 # m
dy_pos = 20e-6 # m

"Initial Gas Concentrations/parameters"
#For now this only applies to the negatrode
#-----Negatrode: For now there is a non reactive node and a reaction node
#Initial Pressures
P_neg_gd = 101325 # Pa
P_neg_rxn = 101325 #81343 # Pa

# Initial mol fractions of negatrode gasses
X_k_gd = np.array([0.50, 0.50, 1e-12]) #H2, N2, Steam ([0.50, 0.49, 0.01])
X_k_rxn = np.array([0.50, 0.50, 1e-12]) #H2, N2, Steam ([0.40, 0.55, 0.05])

#Gas Concentrations
C_k_gd_neg0 = X_k_gd*((P_neg_gd)/(R*T)) #H2, N2, Steam
C_k_rxn_neg0 = X_k_rxn*((P_neg_rxn)/(R*T)) #H2, N2, Steam

#Gas properties:
mu = 2.08e-5 #kg/m-s #I am going to use your value for this (dynamic viscosity)
D_k_an = np.array([0.3e-5, 2.798e-5, 1.9e-5]) #m2/s, H2, N2, Steam (gas diffusion, your values)

"n_values:"
n_neg_p = -1
n_neg_o = -2
n_pos_p = 2
n_pos_o = 2

"Potentials (I will place in more accurate numbers later) The anode is the reference"
phi_neg_0 = 0 #this will by my reference electrode
phi_elyte_0 = 0.5 # I will calculate more accurately later
phi_pos_0 = 1.05 # I will calculate more accurately later

dphi_int_neg_0 = phi_elyte_0-phi_neg_0 #Sets the applied potential on the cell\
dphi_int_pos_0 = phi_pos_0-phi_elyte_0

"Activity concentrations"
#Negatrode Activity Concentrations: 
C_H_Ni = X_H_Ni*C_Ni_s #(mol/m^2)
C_H2O_Ni = X_H2O_Ni*C_Ni_s #(mol/m^2)
C_vac_Ni = X_vac_Ni*C_Ni_s #(mol/m^2)
C_Hx_elyte = X_Hx_elyte*lat_site_BCZYYb #(mol/m^3)
C_Ox_elyte = X_Ox_elyte*lat_site_BCZYYb #(mol/m^3)
C_vac_elyte = X_vac_elyte*lat_site_BCZYYb #(mol/m^3)
#Positrode Activity Concentrations: (mol/m^2)
C_Hx_BF_s = X_Hx_BF*C_BCFZY_s
C_H2O_BF_s = X_H2O_BF*C_BCFZY_s
C_vac_BF_s = X_vac_BF*C_BCFZY_s
C_O_BF_s = X_O_BF*C_BCFZY_s
C_Ox_BF_s = X_Ox_elyte*C_BCFZY_s

"Material Parameters"
#Lattice Site parameters:
lat_site_BCZYYb = 19484 #mol/m^3 7111
#Zhu, H., Ricote, S., Duan, C., O’Hayre, R. P. & Kee, R. J.  Defect Chemistry and Transport within Dense 
#BaCe0.7Zr0.1Y0.1Yb0.1O3−δ(BCZYYb) Proton-Conducting Membranes . J. Electrochem. Soc. 165, F845–F853 (2018).
C_Ni_s = 2.6e-05 #Surface site Concentrations mol/m^2 (again this is just from hw4)
C_BCFZY_s = 2e-5 #Made this up, though I do have the bulk value

"geometric parameters"
n_brugg = -0.5 #bruggman factor assuming alpha is -1.5
#anode
eps_Ni = 0.159 #see calculations
eps_elyte_neg = 0.191 #See Calculations
eps_gas_neg = 1-eps_Ni-eps_elyte_neg
d_Ni_neg = 1*10**-5 #(m)rough estimate from SEM images (average diameter of Ni in negatrode)
d_elyte_neg = 5*10**-6 #(m) rough estimate from SEM images (average diameter of BCZYYb in negatrode)
d_part_avg = (d_Ni_neg+d_elyte_neg)/2 #just taking a linear average of the two particle sizes
r_int = 2*10**-6 #(m) rough estimate from SEM images, interface region between particles, on SEM images it looks almost like the radius
#Cathode
d_BCFZY = 500*10**-9 #(m) rough estimate from SEM images
eps_BCFZY = 0.45 #just assuming 55% porosity need to look up this value could measure more accurately
eps_gas_pos = 1-eps_BCFZY

"Thermodynamic values (first 5 taken from homework 4, last one I had to make up)"
g_H_Ni_o = -7.109209e+07      # standard-state gibbs energy for H adsorbed on Ni surface (J/kmol)
g_H2O_Ni_o = -3.97403035e+08  # standard-state gibbs energy for H2O adsorbed on Ni surface (J/kmol)
g_Vac_Ni_o = 0.0              # standard-state gibbs energy for Ni surface vacancy (J/kmol)
g_Vac_elyte_o = 0.0           # standard-state gibbs energy for electrolyte oxide vacancy (J/kmol)
g_Ox_elyte_o = -2.1392135e+08 # standard-state gibbs energy for electrolyte oxide O2- (J/kmol)
g_Hx_elyte_o = -2.1392135e+07 # standard-state gibbs energy for electrolyte protons H+ (J/kmol)

#/\/\/\/\/\ Initializing Solution Vector /\/\/\/\/\
#SV_0 = np.hstack([dphi_int_neg_0, C_k_gd_neg0, C_k_rxn_neg0, dphi_int_pos_0 ])
#SV_0 = np.hstack([dphi_int_neg_0, C_k_rxn_neg0, dphi_int_pos_0 ])
#SV_s = np.array(['dphi_dl_neg','C_H2_gd','C_N2_gd','C_steam_gd','C_H2_rxn','C_N2_rxn','C_steam_rxn','dphi_dl_pos'])
#SV_0 = np.array([dphi_int_neg_0,dphi_int_pos_0])
SV_0 = np.hstack([dphi_int_neg_0, C_H_Ni, C_H2O_Ni, C_vac_Ni, C_Hx_elyte, 
    C_Ox_elyte, C_vac_elyte, C_k_rxn_neg0, dphi_int_pos_0 ])
# SV including all negatrode reactions:

print(SV_0)

#/\/\/\/\/\ Making the parameter class /\/\/\/\/\
class pars:
    C_k_gd_neg0 = C_k_gd_neg0 #H2, N2, Steam
    #important parameters
    i_ext = i_ext
    C_dl_neg = C_dl_neg
    C_dl_pos = C_dl_pos
    time_span = np.array([0,t_final])
    T = T

    "Gas diffusion parameters"
    #Calculations
    tau_fac_neg = eps_gas_neg**n_brugg
    Kg_neg = (eps_gas_neg**3*d_part_avg**2)/(72*tau_fac_neg*(1-eps_gas_neg)**2)
    #Node thicknesses:
    dy_neg1 = dy1_neg
    dy_neg2 = dy2_neg
    eps_gas_neg = eps_gas_neg
    eps_g_dy_Inv_rxn = 1/(dy2_neg*eps_gas_neg)
    #Gas properties:
    dyn_vis_gas = mu 
    D_k_gas_neg = D_k_an

    #beta values
    beta_o = 0.5 
    beta_p = 0.5

    #Interface Potentials
    dphi_int_neg_0 = phi_elyte_0-phi_neg_0 #Sets the applied potential on the cell\
    dphi_int_pos_0 = phi_pos_0-phi_elyte_0

    "Chemical parameters: (For these I just used yours, not sure where/how to find them) (I also kept hte positrode and negatrode values the same for now)"
    #Negatrode OER
    k_fwd_star_neg_o = k_fwd_star_neg_o
    k_rev_star_neg_o = k_rev_star_neg_o
    #Negatrode HER 
    k_fwd_star_neg_p =  k_fwd_star_neg_p
    k_rev_star_neg_p = k_rev_star_neg_p
    #Positrode ORR
    k_fwd_star_pos_o = k_fwd_star_pos_o
    k_rev_star_pos_o = k_rev_star_pos_o
    #Positrode HRR 
    k_fwd_star_pos_p = k_fwd_star_pos_p
    k_rev_star_pos_p = k_rev_star_pos_p

    "n_values for the negatrode charge transfer reactions:"
    n_neg_p = n_neg_p
    n_neg_o = n_neg_o
    n_pos_p = n_pos_p
    n_pos_o = n_pos_o

    "Geometric parameters"
    #Anode Geometric Parameters
    L_TPB = 2*math.pi*r_int
    A_surf_Ni_neg = 4*math.pi*(d_Ni_neg/2)**2
    LTPB_invANi = L_TPB/A_surf_Ni_neg
    A_fac_Ni = 0.5*eps_Ni*3.*dy2_neg/d_Ni_neg #Converted from your A_fac reaction A of Ni per A negatrode
    A_surf_elyte_neg = 4*math.pi*(d_elyte_neg/2)**2
    tau_fac_neg = eps_gas_neg**n_brugg #tortuosity factor
    Kg_neg = (eps_gas_neg**3*d_part_avg**2)/(72*tau_fac_neg*(1-eps_gas_neg)**2) #gas permeability, see calculations for more details
    #Cathode Geometric Parameters
    A_surf_BCFZY = 4*math.pi*(d_BCFZY/2)**2
    tau_fac_pos = eps_gas_pos**n_brugg
    Kg_pos = (eps_gas_pos**3*d_part_avg**2)/(72*tau_fac_pos*(1-eps_gas_neg)**2)

    #Positrode Surface Activity Concentrations: (mol/m^2)
    C_Hx_BF_s = C_Hx_BF_s
    C_H2O_BF_s = C_H2O_BF_s
    C_vac_BF_s = C_vac_BF_s
    C_O_BF_s = C_O_BF_s
    C_Ox_BF_s = C_Ox_BF_s

    "Positrode Product calculations" #Calculates the product terms in the mass action equations
    prod_fwd_pos_o = C_O_BF_s**-nu_O_BF_pos_o #- signs are needed to cancel out the sign convention of the stoichiometric coefficients
    prod_fwd_pos_p = C_Hx_BF_s**-nu_Hx_BF_pos_p * C_O_BF_s**-nu_O_BF_pos_p
    prod_rev_pos_o = C_Ox_BF_s**nu_Ox_BF_pos_o
    prod_rev_pos_p = C_H2O_BF_s**nu_H2O_BF_pos_p * C_vac_BF_s**nu_vac_BF_pos_p

    "Stoichiometric values For the charge transfer reactions:"
    #negatrode proton reaction:
    nu_H_Ni_neg_p = nu_H_Ni_neg_p
    nu_vac_ely_neg_p = nu_vac_ely_neg_p
    nu_Hx_ely_neg_p = nu_Hx_ely_neg_p
    nu_vac_Ni_neg_p = nu_vac_Ni_neg_p
    #negatrode oxide reaction:
    nu_H_Ni_neg_o = nu_H_Ni_neg_o
    nu_H2O_Ni_neg_o = nu_H2O_Ni_neg_o
    nu_vac_Ni_neg_o = nu_vac_Ni_neg_o
    nu_vac_ely_neg_o = nu_vac_elyte_neg_o 
    nu_Ox_ely_neg_o = nu_Ox_elyte_neg_o
    #postirode proton reaction:
    nu_Hx_BF_pos_p = nu_Hx_BF_pos_p
    nu_O_BF_pos_p = nu_O_BF_pos_p
    nu_H2O_BF_pos_p = nu_H2O_BF_pos_p
    nu_vac_BF_pos_p = nu_vac_BF_pos_p 
    #positrode oxide reaction:
    nu_O_BF_pos_o = nu_O_BF_pos_o
    nu_Ox_BF_pos_o = nu_Ox_BF_pos_o
    nu_vac_BF_pos_o = nu_vac_BF_pos_o
    
#/\/\/\/\/\ Making the pointer class/\/\/\/\/\
#specifies where in SV certain variables are stored
class ptr:
    dphi_int_neg = 0
    
    #Negatrode surface concentrations:
    C_H_Ni = dphi_int_neg+1
    C_H2O_Ni = C_H_Ni+1
    C_vac_Ni = C_H2O_Ni+1
    C_Hx_elyte = C_vac_Ni+1
    C_Ox_elyte = C_Hx_elyte+1
    C_vac_elyte = C_Ox_elyte+1

    ##C_k in negatrode GDL: starts just after dphi_int_neg, is same size as C_k_gd_neg0:
    #C_k_gd_neg = np.arange(dphi_int_neg+1, dphi_int_neg+1+C_k_gd_neg0.shape[0])
    
    #C_k_rxn_neg = np.arange(dphi_int_neg+1, dphi_int_neg+1+C_k_gd_neg0.shape[0])
    
    #C_k_rxn_neg = np.arange(C_k_gd_neg[-1]+1, C_k_gd_neg[-1]+1+C_k_gd_neg0.shape[0])

    C_k_rxn_neg = np.arange(C_vac_elyte+1, C_vac_elyte+1+C_k_gd_neg0.shape[0])
    
    dphi_int_pos = C_k_rxn_neg[-1]+1


SV_0 = np.hstack([dphi_int_neg_0, C_H_Ni, C_H2O_Ni, C_vac_Ni, C_Hx_elyte, 
    C_Ox_elyte, C_vac_elyte, C_k_rxn_neg0, dphi_int_pos_0 ])