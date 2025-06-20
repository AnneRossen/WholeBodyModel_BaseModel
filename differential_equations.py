#!/usr/bin/env python
from production_rates import ProductionRates
import numpy as np

######################################################
####### Corresponds to Peter's Model.m #######
######################################################

def Model(t, X, p):    
    """
    Whole-body model differential equation system
    
    Args:
        p (class): parameters
    """
    # Equations for circulating metabolites
    idx_X = 0
    for organ in p.organs:
        p.C[organ] = X[idx_X : (idx_X+18)]
        idx_X += 18
        
    # Equations for SIMO model
    S = X[-15:-12] # Stomach 
    J = X[-12:-9]  # Jejenum
    R = X[-9:-6]   # Delay
    L = X[-6:-3]   # Ileum
    
    # Equations for insulin release
    P = X[-3]
    II = X[-2]
    QQ = X[-1]
    
    ##################################################
    
    # Calculate production rates vectors and glucagon/insulin submodel
    p = ProductionRates(p, J, L, II, QQ)
    
    # SIMO-vectors for macronutrient uptake
    SIMO = dict()
    SIMO["GI"] = np.zeros(18)
    SIMO["GI"][0] = p.GI["GLC"]
    SIMO["GI"][9] = p.GI["AA"]
    SIMO["MP"] = np.zeros(18)
    SIMO["MP"][11] = p.GI["TGL"]*0.5 # fat diposited directly into muscle and asipose tissue
    SIMO["AP"] = np.zeros(18)
    SIMO["AP"][11] = p.GI["TGL"]*0.5
    
    ##################################################
    
    # Main equations
    dC_Bdt =  (p.m * p.Q["Brain"] @ (p.C["Heart"]-p.C["Brain"])) / p.V["Brain"] + p.R["Brain"]
    dC_Hdt =  (p.m @ (p.Q["Brain"]*p.C["Brain"] + p.Q["Liver"]*p.C["Liver"] + p.Q["Kidney"]*p.C["Kidney"] + p.Q["Muscle"]*p.C["Muscle"] + p.Q["Adipose"]*p.C["Adipose"] - p.Q["Heart"]*p.C["Heart"])) / p.V["Heart"] + p.R["Heart"]
    dC_Gdt =  (p.m * p.Q["Gut"] @ (p.C["Heart"]-p.C["Gut"]) + SIMO["GI"]) / p.V["Gut"] + p.R["Gut"]
    dC_Ldt =  (p.m @ (p.Q["Q_A"]*p.C["Heart"] + p.Q["Gut"]*p.C["Gut"] - p.Q["Liver"]*p.C["Liver"])) / p.V["Liver"] + p.R["Liver"]
    dC_Kdt =  (p.m * p.Q["Kidney"] @ (p.C["Heart"]-p.C["Kidney"])) / p.V["Kidney"] + p.R["Kidney"]
    dC_MPdt = (p.m * p.Q["Muscle"] @ (p.C["Heart"]-p.C["Muscle"]) + SIMO["MP"]) / p.V["Muscle"] + p.R["Muscle"]
    dC_APdt = (p.m * p.Q["Adipose"] @ (p.C["Heart"]-p.C["Adipose"]) + SIMO["AP"]) / p.V["Adipose"] + p.R["Adipose"]
    dCdt = np.concatenate([dC_Bdt, dC_Hdt, dC_Gdt, dC_Ldt, dC_Kdt, dC_MPdt, dC_APdt]).flatten()

    # SIMO model
    dSdt = -p.GI["k_(js)"]*S
    dJdt =  p.GI["k_(js)"]*S - p.GI["k_(gj)"]*J - p.GI["k_(rj)"]*J
    dRdt =  p.GI["k_(rj)"]*J - p.GI["k_(lr)"]*R
    dLdt =  p.GI["k_(lr)"]*R- p.GI["k_(gl)"]*L
    dSIMOdt = np.concatenate([dSdt, dJdt, dRdt, dLdt]).flatten()
    
    # Insulin submodel
    dPdt =  p.I["alpha"] * (p.I["P_inf"]-P)
    dIIdt = p.I["beta"] * (p.I["X"]-II)
    dQQdt = p.I["K"] * (p.I["Q0"]-QQ) + p.I["y"]*P - p.I["S"]
    dIdt = np.concatenate([np.array([dPdt]), np.array([dIIdt]), np.array([dQQdt])]).flatten()

    dzdt = np.concatenate([dCdt, dSIMOdt, dIdt])
    
    # Error handling
    print(f"t={t}")
    if np.any(np.isnan(dzdt)) or np.any(np.isinf(dzdt)):
        print("Warning! NaN or Inf detected in dzdt at t =", t)
        print("dzdt shape:", dzdt.shape, " Min/Max:", np.min(dzdt), np.max(dzdt))
        sys.exit()
   
    return dzdt
