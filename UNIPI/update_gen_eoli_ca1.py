# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:17:25 2023

@author: barto
"""
import numpy as np
from casadi import*
import pandas as pd
import casadi as ca

N = 96
# Funzione che restituisce un nuovo valore basato sull'argomento
def update_gen_eoli_ca1(gen_eoli, alpha):
   
    v = pd.read_excel('DUMMY_beta_01_eq.xlsx', sheet_name='Meteo', usecols='D', skiprows=2, header=None, nrows=800).values.flatten()
    
    print("Valore di velocità", v)

    # Carico i parametri di ingresso
    v_min = gen_eoli['v_min']  # v_{\min}: velocità minima del vento [m/s]
    v_max = gen_eoli['v_max']  # v_{\max}: velocità massima del vento [m/s]
      
    PN = gen_eoli['PN']
    Pmin = 0
    fk_values = gen_eoli['fk_values']
        # v(i): velocità del vento [m/s]

    # Inizializza array fk
    fk = np.zeros(N)
    #fk = SX.zeros(N)
    for i in range(N):
        # Controlla le velocità massime e minime
        fk[i] = ca.if_else(ca.logic_or(v[i] <= v_min, v[i] >= v_max),
                            0,
                            ca.if_else(ca.logic_and(v[i] > v_min, v[i] < fk_values[0][0]),
                                       Pmin + ((fk_values[0][1] - Pmin) / (fk_values[0][0] - v_min)) * (v[i] - v_min),
                                       ca.if_else(ca.logic_and(v[i] > fk_values[-1][0], v[i] < v_max),
                                                  (fk_values[-1][1] + ((PN - fk_values[-1][1]) / (v_max - fk_values[-1][0])) * (v[i] - fk_values[-1][0])),
                                                  0)))

        # Controlla gli intervalli nella curva caratteristica
        for j in range(len(fk_values) - 1):
            fk[i] = ca.if_else(ca.logic_and(v[i] >= fk_values[j][0], v[i] <= fk_values[j + 1][0]),
                                fk_values[j][1] + ((fk_values[j + 1][1] - fk_values[j][1]) / (fk_values[j + 1][0] - fk_values[j][0])) * (v[i] - fk_values[j][0]),
                                fk[i])

  
   
    
    phi = fk
    gen_eoli['alpha'] = alpha
    gen_eoli['G'] = phi * alpha

    gen_eoli['W'] = gen_eoli['G']
    return gen_eoli


   