#This is where I will store my functions to be used in my model


#/\/\/\/\/\ Imports: /\/\/\/\/\
import numpy as np
import math

#Constants:
F = 96485 #C/mol
R = 8.313 #J/mol*k


#/\/\/\/\/\/\ Gas transport functions /\/\/\/\/\/\/\
#So far only contains one function (likely will stay that way)
def electrode_gas_transport(s1,s2,gasProps): #calculates gas diffusion betwee two nodes
    N_k  = np.zeros_like(s1['C_k'])
    #Setting the volume fractions of each layer of the negatrode:
    f1 = s1['dy']/(s1['dy'] + s2['dy'])
    f2 = 1-f1
    C_int = f1*s1['C_k'] + f2*s2['C_k'] #mol/m^3
    #re-finding the mol fractions of the gas constituents
    X_k_1 = s1['C_k']/np.sum(s1['C_k'])
    X_k_2 = s2['C_k']/np.sum(s2['C_k'])
    X_k_int = f1*X_k_1 + f2*X_k_2
    #Calculating the pressure values
    P_1 = np.sum(s1['C_k'])*R*gasProps['T']
    P_2 = np.sum(s2['C_k'])*R*gasProps['T']
    #Calculating V_k_diff
    D_k_eff = gasProps['eps_g']*gasProps['D_k']/gasProps['t_fac'] #(m^2/s)eps_g_neg and tau_fac_neg will be solved for before the function
    dY = 0.5*(s1['dy'] + s2['dy']) #getting the average thickness for each layer
    V_k_diff = -D_k_eff*(X_k_2 - X_k_1)/(dY*X_k_int)
    #Calculating Vconv and V_k)tot
    V_conv = -gasProps['Kg']*(P_2 - P_1)/dY/gasProps['mu'] #Kg_neg will be solved for before the function
    V_k_diff = -D_k_eff*(X_k_2 - X_k_1)/dY/X_k_int #m/s

    V_k  = V_conv + V_k_diff #m/s
        
    N_k = C_int*X_k_int*V_k #mol/(m^2*s)
    return N_k

#/\/\/\/\/\/\ Ion/electron diffusion modeling /\/\/\/\
#def mixed_conduction(s1,s2,matProps): #calculates the ionic and electronic diffusion of the three species

#/\/\/\/\/\/\ Main Modeling function /\/\/\/\/\/\/\
def residual(t, SV, pars, ptr):
    dSV_dt = np.empty_like(SV) #initializing residual (Zeroes_like gave me errors)
    
    #----- Negatrode ------------------------------------------------
    "Charge Transfer"
    dphi_neg = SV[ptr.dphi_int_neg]
    
    #"Activity concentrations"
    #Negatrode Activity Concentrations: 
    C_H_Ni = SV[ptr.C_H_Ni] #(mol/m^2)
    C_H2O_Ni = SV[ptr.C_H2O_Ni] #(mol/m^2)
    C_vac_Ni = SV[ptr.C_vac_Ni] #(mol/m^2)
    C_Hx_elyte = SV[ptr.C_Hx_elyte] #(mol/m^3)
    C_Ox_elyte = SV[ptr.C_Ox_elyte] #(mol/m^3)
    C_vac_elyte = SV[ptr.C_vac_elyte] #(mol/m^3)

    # Product Reactions:
    "Negatrode Product calculations" #Calculates the product terms in the mass action equations
    prod_fwd_neg_o = C_Ox_elyte**-pars.nu_Ox_ely_neg_o * C_H_Ni**-pars.nu_H_Ni_neg_o  #- signs are needed to cancel out the sign convention of the stoichiometric coefficients
    prod_fwd_neg_p = C_H_Ni**-pars.nu_H_Ni_neg_p * C_vac_Ni**-pars.nu_vac_ely_neg_p
    prod_rev_neg_o = C_vac_elyte**pars.nu_vac_ely_neg_o * C_H2O_Ni**pars.nu_H2O_Ni_neg_o * C_vac_Ni**pars.nu_vac_Ni_neg_o
    prod_rev_neg_p = C_Hx_elyte**pars.nu_Hx_ely_neg_p * C_vac_Ni**pars.nu_vac_Ni_neg_p

    #Mass Action Equations For Charge Transfer:
    i_Far_neg_o = pars.n_neg_o*F*(pars.k_fwd_star_neg_o*math.exp((-pars.beta_o*pars.n_neg_o*F*dphi_neg)/(R*pars.T))*prod_fwd_neg_o
        -pars.k_rev_star_neg_o*math.exp(((1-pars.beta_o)*pars.n_neg_o*F*dphi_neg)/(R*pars.T))*prod_rev_neg_o)
    
    i_Far_neg_p = pars.n_neg_p*F*(pars.k_fwd_star_neg_p*math.exp((-pars.beta_p*pars.n_neg_p*F*dphi_neg)/(R*pars.T))*prod_fwd_neg_p
        -pars.k_rev_star_neg_p*math.exp(((1-pars.beta_p)*pars.n_neg_p*F*dphi_neg)/(R*pars.T))*prod_rev_neg_p)
    
    i_Far_neg = i_Far_neg_o + i_Far_neg_p
    
    #Final calculations for the change in potential difference on the negatrode
    i_dl_neg = pars.i_ext - i_Far_neg
    ddphi_int_neg = -i_dl_neg/pars.C_dl_neg
    #print(ddphi_int_neg)
    dSV_dt[ptr.dphi_int_neg] = ddphi_int_neg
    #print(dSV_dt[ptr.dphi_int_neg])

    "Surface charge transfer reactions (finish equations lower down)"
    #sdot calculations (mol/s): (need to make sure this isnt mol/m^2-s)
    sdot_H_Ni = ((pars.nu_H_Ni_neg_p*i_Far_neg_p)/(pars.n_neg_p*F) +
        (pars.nu_H_Ni_neg_o*i_Far_neg_o)/(pars.n_neg_o*F))
    sdot_Hx_elyte = (pars.nu_Hx_ely_neg_p*i_Far_neg_p)/(pars.n_neg_p*F)#same as VOH
    sdot_Ox_elyte = (pars.nu_Ox_ely_neg_o*i_Far_neg_o)/(pars.n_neg_o*F) #same as VO
    sdot_H2O_Ni = (pars.nu_H2O_Ni_neg_o*i_Far_neg_o)/(pars.n_neg_o*F)
    sdot_vac_Ni = ((pars.nu_vac_Ni_neg_p*i_Far_neg_p)/(pars.n_neg_p*F) +
        (pars.nu_vac_Ni_neg_o*i_Far_neg_o)/(pars.n_neg_o*F)) #surface reaction site for a proton
    sdot_vac_elyte = ((pars.nu_vac_ely_neg_p*i_Far_neg_p)/(pars.n_neg_p*F) +
        (pars.nu_vac_ely_neg_o*i_Far_neg_o)/(pars.n_neg_o*F)) #site for a proton
    
    #Hydrogen adsorption (gas phase reaction):
    sdot_H2_gas = 0.5*(((pars.nu_H_Ni_neg_p*i_Far_neg_p)/(pars.n_neg_p*F)) +
     ((i_Far_neg_o*pars.nu_H_Ni_neg_o)/(pars.n_neg_o*F)))#(mol/s)/pars.A_surf_Ni_neg #(mol/(m2*s))int not an array
    #Water desorption(gas phase reaction):
    sdot_H20_gas_neg = (i_Far_neg_o*pars.nu_H2O_Ni_neg_o)/(pars.n_neg_o*F)#/pars.A_surf_Ni_neg #int not an array

    #DCk/dt (mol/m^2-s)
    dSV_dt[ptr.C_H_Ni] = 2*sdot_H2_gas-sdot_H_Ni*pars.LTPB_invANi #DCk_dt_H
    dSV_dt[ptr.C_Hx_elyte] = sdot_Hx_elyte*pars.LTPB_invANi #DCk_dt_Hx 
    #J_Ox = 0.000454 #Temporary from HW4
    dSV_dt[ptr.C_Ox_elyte]=  sdot_Ox_elyte-sdot_H2O_Ni*pars.LTPB_invANi #DCk_dt_Ox  maybe its this J_Ox-sdot_H2O_Ni*pars.LTPB_invANi
    dSV_dt[ptr.C_H2O_Ni]= sdot_H2O_Ni*pars.LTPB_invANi-sdot_H20_gas_neg #DCk_dt_H2O
    dSV_dt[ptr.C_vac_Ni] = -(2*sdot_H2_gas + sdot_H20_gas_neg) + sdot_vac_Ni*pars.LTPB_invANi #DCk_dt_vac_Ni
    dSV_dt[ptr.C_vac_elyte]= sdot_vac_elyte*pars.LTPB_invANi# DCk_dt_vac_elyte

    "Negatrode Gas Transport"
    #Getting parameters from the SV
    #C_k_gd_neg = SV[ptr.C_k_gd_neg]
    C_k_gd_neg = pars.C_k_gd_neg0 #SV[ptr.C_k_gd_neg]
    C_k_rxn_neg = SV[ptr.C_k_rxn_neg]
    #Making dictionaries for the gas diffusion equation:
    s1 = {'C_k':C_k_gd_neg,'dy':pars.dy_neg1}
    s2 = {'C_k':C_k_rxn_neg,'dy':pars.dy_neg2}
    gasProps_neg = {'Kg':pars.Kg_neg,'t_fac':pars.tau_fac_neg,'D_k':pars.D_k_gas_neg, 
                        'mu':pars.dyn_vis_gas,'T':pars.T,'eps_g':pars.eps_gas_neg}
    #Running gas transport function
    N_k_i = electrode_gas_transport(s1,s2,gasProps_neg) #mol/m^3
    "Negatrode Gas Phase Reactions"

    #final gas phase equation
    sdot_k = np.array([sdot_H2_gas,0,sdot_H20_gas_neg]) #(mol/(m2*s))hydrogen, nitrogen, water
    #print(sdot_k)
    #print(N_k_i)
    dCk_dt = (N_k_i + sdot_k*pars.A_fac_Ni)*pars.eps_g_dy_Inv_rxn
    #print(dCk_dt)
    #print(SV[ptr.C_k_gd_neg])
    dSV_dt[ptr.C_k_rxn_neg] = dCk_dt
    #q
    #----- Positrode ------------------------------------------------
    dphi_pos = SV[ptr.dphi_int_pos]
    #print(dphi_pos)
    #MA Exp equations:
    exp_fwd_pos_o = math.exp((-pars.beta_o*pars.n_pos_o*F*dphi_pos)/(R*pars.T))
    exp_rev_pos_o = math.exp(((1-pars.beta_o)*pars.n_pos_o*F*dphi_pos)/(R*pars.T))
    #print(exp_fwd_pos_o,exp_rev_pos_o)
    #Mass-Action equations
    i_Far_pos_o = pars.n_pos_o*F*(pars.k_fwd_star_pos_o*exp_fwd_pos_o*pars.prod_fwd_pos_o
        -pars.k_rev_star_pos_o*exp_rev_pos_o*pars.prod_rev_pos_o)
    i_Far_pos_p = pars.n_pos_p*F*(pars.k_fwd_star_pos_p*math.exp((-pars.beta_p*pars.n_pos_p*F*dphi_pos)/(R*pars.T))*pars.prod_fwd_pos_p
        -pars.k_rev_star_pos_p*math.exp(((1-pars.beta_p)*pars.n_pos_p*F*dphi_pos)/(R*pars.T))*pars.prod_rev_pos_p)
    i_Far_pos = i_Far_pos_o + i_Far_pos_p
    
    i_dl_pos = pars.i_ext - i_Far_pos
    ddphi_int_pos = -i_dl_pos/pars.C_dl_pos
    dSV_dt[ptr.dphi_int_pos] = ddphi_int_pos
    return dSV_dt