# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:32:50 2023

@author: barto
"""
import pandas as pd
import numpy as np
from update_gen_comb import update_gen_comb  
timespan = pd.read_excel('DUMMY_beta_01_eq.xlsx', sheet_name='Meteo', usecols='A', skiprows=2, header=None, nrows=800).values.flatten()
tau = 15 
N = len(timespan)
def carica_gen_comb(excel_filename, sheetname_gen_comb):
    dati_gen_comb = pd.read_excel(excel_filename, sheet_name='gen_comb', usecols='B', skiprows=5, nrows=6, header=None).values.flatten()
    gen_comb = {}
    gen_comb['comb'] = pd.read_excel(excel_filename, sheet_name=sheetname_gen_comb, usecols='B', skiprows=4, header=None).iloc[0, 0]
    gen_comb['sheetname'] = sheetname_gen_comb
    gen_comb['PCI'] = dati_gen_comb[0]  # PCI(k): potere calorifico inferiore del combustibile [kJ/kg]
    gen_comb['cF'] = dati_gen_comb[1]#*(tau/60);
    gen_comb['PN'] = dati_gen_comb[2]  # PN(k) : potenza nominale [kWe]
    gen_comb['N_acc_max'] = dati_gen_comb[3]  # Numero massimo di accensioni [-]
    gen_comb['c_acc'] = dati_gen_comb[4]  # Coefficiente di penalizzazione accensione [€/N_acc]
    gen_comb['P_min'] = dati_gen_comb[5]  # P_{\min}(k) : potenza minima tecnica [kWe]

    eta_values_raw = pd.read_excel(excel_filename, sheet_name=sheetname_gen_comb, usecols="B:C", skiprows=20, nrows=10, header=None)
    gen_comb['eta_values'] = [eta_values_raw.iloc[:, 0].values / gen_comb['PN'], eta_values_raw.iloc[:, 1].values]

    gen_comb['alpha'] = pd.read_excel(excel_filename, sheet_name=sheetname_gen_comb, usecols="M", skiprows=4, header=None, nrows=800).values.flatten()
    gen_comb['n_x'] = N
    gen_comb = update_gen_comb(gen_comb, gen_comb['alpha'])
    
    return gen_comb
   