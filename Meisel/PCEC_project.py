"""
    This file runs and executes a PCEC model, calculating the faradaic efficiency and ohmic resistance as a function of the electrolyte thickness. 

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Initializes the model by calling pcec_params.py
            a - pcec_params.py reads in the user inputs from pemfc_inputs.py, 
            b - pcec_params.py then creates an initial solution vector 'SV_0'
            c - pcec_params.py returns SV_0 and a class 'pars' holding all simulation parameters.
        
        2 - Calls the residual function and integrates over the user-defined    
            time span. (will change this to only looks at steady state in the future)

        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""
# Import necessary modules:
from scipy.integrate import solve_ivp #integration function for ODE system.
from pcec_functions import residual # point the model to the residual function
from pcec_params import pars, SV_0, ptr #imports parameters, the initial SV and the pointer class
from matplotlib import pyplot as plt


def pcec_model():
    # The use of the 'lambda' function is required here so that we can pass the 
    #   class variablels 'pars' and 'ptr.'  Otherwise, we can only pass the 
    #   time span and our initial solution SV_0:
    solution = solve_ivp(lambda t, y: residual(t, y, pars, ptr),
        pars.time_span, SV_0, rtol=1e-8, atol=1e-6, method='BDF')
    return solution

solution = pcec_model()
print(solution)
#for var in solution.y:
    #plt.plot(solution.t,var)
    #plt.show()

#plt.plot(solution.t,solution.y[0])
#plt.xlabel("time (s)")
#plt.ylabel('Dphi_dl_neg')

plt.plot(solution.t,solution.y[4])
plt.xlabel("time (s)")
plt.ylabel('C_H+_neg')
plt.show()
