#!/usr/bin/env python
import pandas as pd
import numpy as np
import math, json

###################################################################
## Corresponds to Peter's LoadParmModel.m and ModelBasalValues.m ##
###################################################################

def np_convert(df):
    """
    Convert pandas df into a dictionary like numpy
    """
    structured_df = np.array(
        [tuple(row) for row in df.to_numpy()],
        dtype=[(col, "f8") for col in df.columns]  # Assign column names
        ) 
    return structured_df    

class Parameters:
    """
    Class to contain all parameters
    """
    def __init__(self):
        self.organs = None
        self.metabolites = None
        self.S = None
        self.m = None
        self.T = None
        self.Q = None
        self.V = None
        self.Vm = None
        self.Km = None
        self.GI = None
        self.I = None
        self.gamma = None
        self.Gamma = dict()
        self.mu = None
        self.C = None
        self.r = None
        self.R = None
        self.load_param()
    
    
    def load_param(self):
        """
        Loads values from excel data sheets into pandas dfs.
        Before utilizing efficiency_convertion() you will have a nice overview of the data
        """
        # Stoichiometric_matrix
        S = pd.read_excel('data/reaction_kinetics_OG.xlsx', sheet_name="stoichiometric_matrix", header=0, index_col=0)
        S = S.drop(S.columns[0], axis=1)

        # Circulating metabolites
        m = pd.read_excel('data/reaction_kinetics_OG.xlsx', sheet_name="circulating_metabolites", header=0, index_col=0)
        m = pd.DataFrame(np.diag(m.values.flatten()), index=m.columns, columns=m.columns)
                
        # T-vector: Index of reaction
        T = pd.read_excel('data/reaction_kinetics_OG.xlsx', sheet_name="tissue_reactions", header=0, index_col=0)
        T = T.drop(T.columns[0], axis=1)

        # Q Flowrates and V Volumes of compartments in main model
        Q = pd.read_excel('data/parameters_OG.xlsx', sheet_name="FlowRate", header=0, index_col=0, usecols=[0, 1])["Value"]
        V = pd.read_excel('data/parameters_OG.xlsx', sheet_name="Volume", index_col=0, usecols=[0, 1])["Value"]

        # Model Kinetics (maximum velocities)
        Vm = pd.read_excel('data/reaction_kinetics_OG.xlsx', sheet_name="Vm", header=0, index_col=0)
        Vm = Vm.drop(Vm.columns[0], axis=1)
        organs = list(Vm.columns)
        for organ in organs:
            Vm[organ] = Vm[organ] / V[organ]
        
        # Model Kinetics (limiting velocities)
        Km = pd.read_excel('data/reaction_kinetics_OG.xlsx', sheet_name="Km", header=0, index_col=0)
        Km = Km.drop(Km.columns[0], axis=1)
        Km2 = dict()
        for organ in organs:
            Km2[organ] = np.zeros((len(Km), 3))
            for i in range(0, len(Km)):
                if isinstance(Km[organ][i], str):
                    Km2[organ][i] = json.loads(Km[organ][i])
                else:
                    Km2[organ][i][0] = Km[organ][i]
        
    
        # SIMO-Model parameters
        GI = pd.read_excel('data/parameters_OG.xlsx', sheet_name="SIMOmodel", index_col=0, usecols=[0, 1])["Value"]

        # Insulin sub model parameters
        I = pd.read_excel('data/parameters_OG.xlsx', sheet_name="Insulin", index_col=0, usecols=[0, 1])["Value"]

        # Glucagon sub model parameters
        gamma = pd.read_excel('data/parameters_OG.xlsx', sheet_name="Glucagon", index_col=0, usecols=[0, 1])["Value"]

        # Hormonal Control parameters
        mu = pd.read_excel('data/reaction_kinetics_OG.xlsx', sheet_name="mu", header=0, index_col=0)
        mu = mu.drop(mu.columns[0], axis=1)
        
        # Save
        self.organs = np.array(organs)
        self.metabolites = np.array(S.columns)
        self.S = S
        self.m = m
        self.T = T
        self.Q = Q
        self.V = V
        self.Vm = Vm
        self.Km = Km2
        self.GI = GI
        self.I = I
        self.gamma = gamma
        self.mu = mu
        

    def get_initial_values(self, G0, I0, Gamma0):
        """
        Calculate initial values for model initialization
        """
        gamma = self.gamma
        I = self.I
        Q = self.Q 

        # Hearth glucose and insulin at t=0
        gamma["G_B_H"] = G0
        gamma["I_B_H"] = I0    
        
        # Circulating insulin initial values
        I["Heart"] = I0
        I["Brain"] = I0
        I["Gut"] = I0
        I["Kidney"] = I0 * (1-I["F_(KIC)"])
        I["Muscle"] = I0 * (1-I["F_(PIC)"])
        I["Adipose"] = I0 * (1-I["F_(PIC)"])
        I["Liver"] = 1 / Q["Liver"] * (Q["Heart"]*I["Heart"] - Q["Brain"]*I["Brain"] - Q["Kidney"]*I["Kidney"] - Q["Muscle"]*I["Muscle"] - Q["Adipose"]*I["Adipose"])
    
        # Glucagon submodel
        gamma["r_(B_PGammaR)"] = gamma["r_(MGammaC)"] * Gamma0
        gamma["r_(PGammaC)"] = gamma["r_(MGammaC)"] * Gamma0
        gamma["M_G_(PGammaR)"] = 2.93 - 2.10 * math.tanh(4.18 * (G0 / gamma["G_B_H"] - 0.61))
        gamma["M_I_(PGammaR)"] = 1.31 - 0.61 * math.tanh(1.06 * (I["Heart"] / gamma["I_B_H"] - 0.47))
        gamma["r_PgammaR"] = gamma["M_G_(PGammaR)"] * gamma["M_I_(PGammaR)"] * gamma["r_(B_PGammaR)"]
        
        # Circulating glucagon
        gamma["Heart"] = Gamma0
        gamma["Brain"] = Gamma0
        gamma["Gut"] = Gamma0
        gamma["Kidney"] = Gamma0
        gamma["Muscle"] = Gamma0
        gamma["Adipose"] = Gamma0
        gamma["Liver"] = (gamma["r_PgammaR"] - gamma["r_(PGammaC)"] + Q["Q_A"] * gamma["Heart"] + Q["Gut"] * gamma["Gut"]) / Q["Liver"]
        
        # Hormonal control
        gamma["L_SS"] = gamma["Liver"]
        gamma["AP_SS"] = gamma["Adipose"]
        I["L_SS"] = I["Liver"]
        I["MP_SS"] = I["Muscle"]
        I["AP_SS"] = I["Adipose"]
        
        ###########################################################
        
        # 141 differential equations
        x0 = np.ones(len(self.organs)*18+15)

        # Metabolites and intialize C
        init_metabolites = pd.read_excel('data/inital_values_OG.xlsx', header=0, index_col=0)
        init_metabolites = init_metabolites.drop(init_metabolites.columns[0], axis=1)
        n = init_metabolites.shape[0]
        C = pd.DataFrame(0, index=self.S.columns, columns=init_metabolites.columns)
        idx_x = 0
        for organ in self.organs:
            x0[idx_x : (idx_x+n)] = np.append(init_metabolites[organ][:-2], [I[organ], gamma[organ]])
            C[organ] = list(init_metabolites[organ])
            idx_x += n
            
        # SIMO model
        x0[-15:-3] = 0 #Set to zero as no meal is anticipated to be in system initially

        # Insulin submodel
        X_B = (18*G0)**I["beta_(PIR1)"] / (I["beta_(PIR2)"]**I["beta_(PIR1)"] + I["beta_(PIR3)"] * (18*G0)**I["beta_(PIR4)"])
        P_inf = X_B**I["beta_(PIR5)"]
        Y_B = P_inf
        H = I["K"]
        # Differential equations affecting insulin release
        x0[-3] = P_inf                                               #P
        x0[-2] = X_B                                                 #II
        x0[-1] = (H * I["Q0"] + I["y"] * P_inf)/(H + I["M1"] * Y_B)  #QQ
    
        ###########################################################
        
        # Basal values kept for later use
        I["S_B"] = I["M1"] * Y_B * x0[-1]
        I["r_B_PIR"] = Q["Liver"] / (1 - I["F_(LIC)"]) * I["Liver"] - Q["Gut"] * I["Gut"] - Q["Q_A"] * I["Heart"]
        
        # Save
        self.gamma = gamma
        self.I = I
        self.Q = Q 
        self.C = C
        self.r = pd.DataFrame(0, index=self.S.index, columns=self.T.columns)
        self.R = pd.DataFrame(0, index=self.S.columns, columns=self.T.columns)
        return x0
    
    def efficiency_convertion(self):
        """
        Pandas df present data in a nice way, but np is much more effiecient 
        """
        #TODO: dummy proof this function
        try: 
            self.S = self.S.to_numpy()
            self.m = self.m.to_numpy()
            self.T = np_convert(self.T)
            self.Q = self.Q.to_dict()
            self.V = self.V.to_dict()
            self.Vm = np_convert(self.Vm)
            #self.Km = self.Km.to_numpy()
            self.GI = self.GI.to_dict()
            self.I = self.I.to_dict()
            self.gamma = self.gamma.to_dict()
            self.mu = np_convert(self.mu)
            self.C = np_convert(self.C)
            self.r = np_convert(self.r)
            self.R = np_convert(self.R)
        except:
            print("Error while converting parameters")
