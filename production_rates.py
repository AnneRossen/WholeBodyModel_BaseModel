#!/usr/bin/env python
import numpy as np
import sys, math


######################################################
## Corresponds to Peter's ProductionRates.m and r.m ##
######################################################

def get_r(p, organ):
    """
    Calulates reaction rate r(C) using Michaelis Menten kinetics
    
    Args:
        p (class): parameters
        organ (str): organ to calculate r for
    """
    C = p.C[organ]
    Vm = p.Vm[organ]
    Km = p.Km[organ]
    rates = np.zeros(p.S.shape[0])
   
    for i in np.arange(0, len(rates), 1):
        if Vm[i] != 0: # Reactions that occour
            idx_m = np.where(p.S[i,:] < 0)[0]
            
            # Standard form Michaelis-Menten (V*[A] / K1+[A])
            if len(idx_m) == 1: 
                numerator = Vm[i] * C[idx_m[0]]
                denominator = Km[i][0] + C[idx_m[0]]
                rates[i] = numerator/denominator
    
            # Two substrate Michaelis-Menten (V*[A]*[B] / K1+K2*[A]+K3*[B]+[A]*[B])
            elif len(idx_m) > 1: 
                numerator = Vm[i] * C[idx_m[0]] * C[idx_m[1]]
                denominator = Km[i][0] + Km[i][1]*C[idx_m[0]] + Km[i][2]*C[idx_m[1]] + C[idx_m[0]]*C[idx_m[1]]
                rates[i] = numerator/denominator
    return rates


def ProductionRates(p, J, L, II, QQ):
    """
    Calculates concentration dependent parameters and rates
    
    Args:
        p (class): parameters
        J (list): Jejenum [GLC, AA, TGL] concentration
        L (list): Ileum [GLC, AA, TGL] concentration
        II (float): ?
        QQ (float): ?
    """
    #print("Caluclating production rates")
    for organ in p.organs:
        # Michaelis-Menten reaction rates
        p.r[organ] = get_r(p, organ)
        
        # Insulin and Glucagon ratios.
        p.I[organ] = p.C[organ][-2]
        p.Gamma[organ] = p.C[organ][-1]
    
    # Insulin submodel
    p.I["X"] = (18 * p.C["Heart"][0])**p.I["beta_(PIR1)"] / (p.I["beta_(PIR2)"]**p.I["beta_(PIR1)"] + p.I["beta_(PIR3)"] * (18*p.C["Heart"][0])**p.I["beta_(PIR4)"])
    p.I["Y"] = p.I["X"]**p.I["beta_(PIR5)"]
    p.I["P_inf"] = p.I["Y"]
    p.I["S"] = (p.I["M1"]*p.I["Y"] + p.I["M2"]*min(0, p.I["X"]-II)) * QQ
    p.I["r_PIR"] = p.I["S"] / p.I["S_B"] * p.I["r_B_PIR"]
    
    # Insulin clearance
    p.I["r_(LIC)"] = p.I["F_(LIC)"] * (p.Q["Q_A"]*p.I["Heart"] + p.Q["Gut"]*p.I["Gut"] + p.I["r_PIR"])
    p.I["r_(KIC)"] = p.I["F_(KIC)"] * p.Q["Kidney"] * p.I["Heart"]
    p.I["r_(MIC)"] = p.I["F_(PIC)"] * p.Q["Muscle"] * p.I["Heart"]
    p.I["r_(AIC)"] = p.I["F_(PIC)"] * p.Q["Adipose"] * p.I["Heart"]

    # Glucagon submodel
    p.gamma["r_(PGammaC)"] = p.gamma["r_(MGammaC)"] * p.Gamma["Heart"]
    p.gamma["M_G_(PGammaR)"] = 2.93 - 2.10 * math.tanh(4.18 * (p.C["Heart"][0] / p.gamma["G_B_H"] - 0.61))
    p.gamma["M_I_(PGammaR)"] = 1.31 - 0.61 * math.tanh(1.06 * (p.I["Heart"] / p.gamma["I_B_H"] - 0.47))
    p.gamma["r_(PgammaR)"] = p.gamma["M_G_(PGammaR)"] * p.gamma["M_I_(PGammaR)"] * p.gamma["r_(B_PGammaR)"]
    
    # SIMO-model
    [p.GI["GLC"], p.GI["AA"], p.GI["TGL"]] = p.GI["k_(gj)"]*J + p.GI["k_(gl)"]*L

    ##################################################
    
    # Liver Hormonal influence 
    p.r["Liver"][0]  *= (p.I["Liver"]/p.I["L_SS"]) **p.mu["Liver"][0]
    p.r["Liver"][2]  *= (p.I["Liver"]/p.I["L_SS"] * p.gamma["L_SS"]/p.Gamma["Liver"]) **p.mu["Liver"][2]
    p.r["Liver"][3]  *= (p.I["L_SS"]/p.I["Liver"] * p.Gamma["Liver"]/p.gamma["L_SS"]) **p.mu["Liver"][3] # NaN ISSUE
    p.r["Liver"][4]  *= (p.I["Liver"]/p.I["L_SS"] * p.gamma["L_SS"]/p.Gamma["Liver"]) **p.mu["Liver"][4]
    p.r["Liver"][5]  *= (p.I["L_SS"]/p.I["Liver"] * p.Gamma["Liver"]/p.gamma["L_SS"]) **p.mu["Liver"][5] # NaN ISSUE
    p.r["Liver"][12] *= (p.I["Liver"]/p.I["L_SS"]) **p.mu["Liver"][12]
    p.r["Liver"][22] *= (p.I["Liver"]/p.I["L_SS"]) **p.mu["Liver"][22]
    p.r["Liver"][29]  = (p.I["r_PIR"] - p.I["r_(LIC)"]) /p.V["Liver"]
    p.r["Liver"][30]  = (p.gamma["r_(PgammaR)"] - p.gamma["r_(PGammaC)"]) /p.V["Liver"]
    
    # Kidney Hormonal influence
    p.r["Kidney"][29] = -p.I["r_(KIC)"]/p.V["Kidney"]
    
    # Muscle Hormonal influence
    p.r["Muscle"][0]  *= (p.I["Muscle"]/p.I["MP_SS"]) **p.mu["Muscle"][0]
    p.r["Muscle"][2]  *= (p.I["Muscle"]/p.I["MP_SS"]) **p.mu["Muscle"][2]
    p.r["Muscle"][4]  *= (p.I["Muscle"]/p.I["MP_SS"]) **p.mu["Muscle"][4]
    p.r["Muscle"][5]  *= (p.I["MP_SS"]/p.I["Muscle"]) **p.mu["Muscle"][5]
    p.r["Muscle"][12] *= (p.I["Muscle"]/p.I["MP_SS"]) **p.mu["Muscle"][12]
    p.r["Muscle"][25] *= (p.I["Muscle"]/p.I["MP_SS"]) **p.mu["Muscle"][25]
    p.r["Muscle"][26] *= (p.I["MP_SS"]/p.I["Muscle"]) **p.mu["Muscle"][26]
    p.r["Muscle"][29]  = -p.I["r_(MIC)"]/p.V["Muscle"]
    
    # Adipose Hormonal influence
    p.r["Adipose"][0]  *= (p.I["Adipose"]/p.I["AP_SS"]) **p.mu["Adipose"][0]
    p.r["Adipose"][20] *= (p.I["Adipose"]/p.I["AP_SS"]) **p.mu["Adipose"][20]
    p.r["Adipose"][27] *= (p.I["Adipose"]/p.I["AP_SS"]) **p.mu["Adipose"][27]
    p.r["Adipose"][28] *= (p.I["AP_SS"]/p.I["Adipose"] * p.Gamma["Adipose"]/p.gamma["AP_SS"]) **p.mu["Adipose"][28]
    p.r["Adipose"][29]  = -p.I["r_(AIC)"]/p.V["Adipose"]
    
    # Production rate vectors (S*r)
    for organ in p.organs:
        p.R[organ] = p.S.T @ p.r[organ]

        # Error handling
        if np.any(np.isnan(p.r[organ])) or np.any(np.isinf(p.r[organ])):
            print("Warning! NaN or Inf detected in r for organ: ", organ)
            print(p.r[organ])
            sys.exit()        
    return p

