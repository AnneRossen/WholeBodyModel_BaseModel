#!/usr/bin/env python
# %%
from LoadParameters import *
from differential_equations import *
from plot_save_results import *
from digitalisation import *
from scipy.integrate import solve_ivp
import sys 

######################################################
####### Corresponds to Peter's Driver_Model.m #######
######################################################

# Set basal values for glucose, insulin and glucagon
G0 = 5 #[mmol/L]
I0 = 15.1765 #[mU/L]
Gamma0 = 100 #[ng/L]

# Time and Food eaten at time t=0 [GLC,AA,TGL]
t_span = [0, 60*72] #[min]
Food = np.array([0.6,0.24,0.16])*100 #[g] (total 100g) 
Food *= np.array([1/180,1/89.1,1/860])*1000 # Converting [g] to [mmol]

######################################################

# Load parameters
p = Parameters()
x0 = p.get_initial_values(G0, I0, Gamma0)
x0[-15:-12] = Food
#vars(p)

# Solve differential equation system
scale_factor =  np.max(np.abs(x0))
p.efficiency_convertion()
solution = solve_ivp(Model, t_span, x0, args=(p,), method='BDF'
                     #t_eval=t_eval, rtol=1e-3, atol=1e-6
                     )

# Extract solutions
t = solution.t.flatten()
sol = solution.y
print(solution.message)

######################################################

# Display Results: "Figure 8"
plot_results(t, sol, p,
            organ_s = ["Liver"],
            metabolite_s = ["GLC", "AA", "TGL", "GLY", "FFA", "GLR"],
            digitalisation_points = [Fig4_1L, Fig4_1R, Fig4_2L, Fig4_2R, Fig4_3L, Fig4_3R])


# Display GI sub-model
plot_results(t[t<25*60], sol[:,t<25*60], p,
            organ_s = ["Heart"],
            metabolite_s = ["GLC", "TGL", "INS"],
            digitalisation_points = [McQuiad_GLC, McQuiad_TG, McQuiad_INS])


# Display Results: ALL
plot_results(t, sol, p,
            organ_s = p.organs,
            metabolite_s = p.metabolites)            
