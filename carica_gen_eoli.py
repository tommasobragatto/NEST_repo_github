# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:20:37 2023

@author: barto
"""
import pandas as pd
import numpy as np
from update_gen_eoli import update_gen_eoli

N = 96 #len(timespan)
def carica_gen_eoli(excel_filename, sheetname_gen_eoli):   
    dati_gen_eoli = pd.read_excel(excel_filename, sheet_name='gen_eoli', usecols="B", skiprows=4, nrows=4, header=None).values.flatten()
# Dati di targa della macchina della K-esima macchina
    gen_eoli = {}
    gen_eoli['sheetname'] = sheetname_gen_eoli
    gen_eoli['v_min'] = dati_gen_eoli[0]  # Accedi al valore nella prima riga della prima colonna
    gen_eoli['rp'] = dati_gen_eoli[1]     # Accedi al valore nella seconda riga della prima colonna
    gen_eoli['v_max'] = dati_gen_eoli[2]  # Accedi al valore nella terza riga della prima colonna
    gen_eoli['PN'] = dati_gen_eoli[3]     # Accedi al valore nella quarta riga della prima colonna


#  LA FUNZIONE CARATTERISTICA
    fk_values_raw = pd.read_excel(excel_filename, sheet_name=sheetname_gen_eoli, usecols="B:C", skiprows=11, nrows=10, header=None)
    gen_eoli['fk_values'] = np.vstack((
    [gen_eoli['v_min'], 0],
    fk_values_raw,
    [gen_eoli['rp'], gen_eoli['PN']],
    [gen_eoli['v_max'], gen_eoli['PN']]
))
   # gen_eoli['fk_values'] = fk_values

# Assegna un alpha di default
    alpha = pd.read_excel(excel_filename, sheet_name=sheetname_gen_eoli, usecols="M", skiprows=4, nrows=800, header=None)
    gen_eoli['alpha'] = alpha.values.flatten()  # questo dovrebbe essere un feasible initial guess

# Chiamata alla funzione di aggiornamento
    gen_eoli = update_gen_eoli(gen_eoli, gen_eoli['alpha'])

# Definizione del numero di variabili di ottimizzazione
    gen_eoli['n_x'] = N
    
    return gen_eoli
    
   

