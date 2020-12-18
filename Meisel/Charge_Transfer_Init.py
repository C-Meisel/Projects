"This is where I will store all of the initial variables for the Charge Transfer Reactions"
"The params file was just getting too crowded"
import numpy as np

"beta values (need to also look up)"
#beta_o = 0.5 
#beta_p = 0.5

"Chemical parameters: (For these I just used yours, not sure where/how to find them) (I also kept hte positrode and negatrode values the same for now)"
#Negatrode OER
k_fwd_star_neg_o = 4.16307062e+1 # Chemical forward rate constant, m^4/mol^2/s
k_rev_star_neg_o = 4.0650045e-1 #Chemical reverse rate constant, m^4/mol^2/s
#Negatrode HER also neeed to look these up, but im assuming they are much faster than the oxide ones
k_fwd_star_neg_p = 4.16307062e+3 # Chemical forward rate constant, m^4/mol^2/s
k_rev_star_neg_p = 4.0650045e+1 #Chemical reverse rate constant, m^4/mol^2/s
#Positrode ORR
k_fwd_star_pos_o = 4.16307062e+0 # Chemical forward rate constant, m^4/mol^2/s
k_rev_star_pos_o = 4.0650045e-2 #Chemical reverse rate constant, m^4/mol^2/s
#Positrode HRR also neeed to look these up, but im assuming they are much faster than the oxide ones
k_fwd_star_pos_p = 4.16307062e+2 # Chemical forward rate constant, m^4/mol^2/s
k_rev_star_pos_p = 4.0650045e+0 #Chemical reverse rate constant, m^4/mol^2/s

"Concentrations/activities: I need to look these up so I used yours from HW4."
#-----Negatrode:
#Mol fractions (no units)
X_H_Ni = 0.6 #HW4
X_H2O_Ni = 0.2 #HW4
X_vac_Ni = 0.2 #HW4
X_Ox_elyte = 0.8 #I know this is 0.8
X_Hx_elyte = 0.1 #I am unsure of the split between Hx and oxygen vacancies
X_vac_elyte = 0.1 
#-----Positrode:
#Mol fractions (no units) #I made these up, all I know is that 80% of the oxygen lattice sites contain Oxygen:
X_Hx_BF = 0.05
X_H2O_BF = 0.05
X_vac_BF = 0.05
X_O_BF = 0.05
X_Ox_BF = 0.8

"Stoichiometric values For the charge transfer reactions:"
#negatrode proton reaction:
nu_H_Ni_neg_p = -1
nu_vac_ely_neg_p = -1
nu_Hx_ely_neg_p = 1
nu_vac_Ni_neg_p = 1
#negatrode oxide reaction:
nu_H_Ni_neg_o = -2
nu_H2O_Ni_neg_o = 1
nu_vac_Ni_neg_o = 1
nu_vac_elyte_neg_o = 1
nu_Ox_elyte_neg_o = -1
#postirode proton reaction:
nu_Hx_BF_pos_p = -2
nu_O_BF_pos_p = -1
nu_H2O_BF_pos_p = 1
nu_vac_BF_pos_p = 1
#positrode oxide reaction:
nu_O_BF_pos_o = -1
nu_Ox_BF_pos_o = 1
nu_vac_BF_pos_o = 1

"Material Parameters"
#Lattice Site parameters:
lat_site_BCZYYb = 19484 #mol/m^3 7111
#Zhu, H., Ricote, S., Duan, C., O’Hayre, R. P. & Kee, R. J.  Defect Chemistry and Transport within Dense 
#BaCe0.7Zr0.1Y0.1Yb0.1O3−δ(BCZYYb) Proton-Conducting Membranes . J. Electrochem. Soc. 165, F845–F853 (2018).
C_Ni_s = 2.6e-05 #Surface site Concentrations mol/m^2 (again this is just from hw4)
C_BCFZY_s = 2e-5 #Made this up, though I do have the bulk value
