# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:10:58 2023

@author: barto
"""
import pandas as pd
import numpy as np
from update_gen_foto import update_gen_foto

timespan = pd.read_excel('DUMMY_beta_01_eq.xlsx', sheet_name='Meteo', usecols='A', skiprows=2, nrows=800, header=None).values.flatten()
N = len(timespan)

def carica_gen_foto(excel_filename, sheetname_gen_foto):

    dati_gen_foto = pd.read_excel(excel_filename, sheet_name=sheetname_gen_foto, usecols="B", skiprows=4, nrows=8,header=None).values.flatten()
    
    gen_foto = {}
    gen_foto['sheetname'] = sheetname_gen_foto
    gen_foto['type'] = 'gen_foto'
    gen_foto['DNIr'] = dati_gen_foto[0]  # DNI_{r}: irraggiamento di riferimento [1000 W/m^2]
    gen_foto['gamma'] = dati_gen_foto[1]  # \gamma(k): coefficiente di correzione della potenza [°C^-1]
    gen_foto['Tr'] = dati_gen_foto[2]  # T_r: temperatura di riferimento [°C]
    gen_foto['c'] = dati_gen_foto[3]  # c(k): coefficiente di tilt per il pannello [-]
    gen_foto['r'] = dati_gen_foto[4]  # r(k): coefficiente correttivo di Ross [°C m^2/W]
    gen_foto['eta'] = dati_gen_foto[5]  # \eta(k): rendimento complessivo di sistema [-]
    gen_foto['Delta_Tr'] = dati_gen_foto[6]  # \Delta T_{r}(k): differenza di temperatura standard [°C]
    gen_foto['PN'] = dati_gen_foto[7]  # PN(k): potenza nominale [kWe]
    
    gen_foto['G'] = np.zeros(N)
  
    
    # Assegna un alpha di default
    alpha = pd.read_excel(excel_filename, sheet_name=sheetname_gen_foto, usecols="M", skiprows=4, nrows=800,header=None)
    gen_foto['alpha'] = alpha.values.flatten()  # questo dovrebbe essere un feasible initial guess
    
    # Chiamata alla funzione di aggiornamento
    gen_foto = update_gen_foto(gen_foto, gen_foto['alpha'])
    
    # Definizione del numero di variabili di ottimizzazione
    gen_foto['n_x'] = N
    
    return gen_foto
