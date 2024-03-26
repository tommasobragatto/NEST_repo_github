# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:00:23 2023

@author: barto
"""
import numpy as np
from casadi import*
import pandas as pd
import casadi as ca

N = 96

def update_gen_foto(gen_foto, alpha):

    # Carico i parametri di ingresso
    DNIr = gen_foto['DNIr']  # DNI_r: irraggiamento di riferimento [1000 W/m^2]
    gamma = gen_foto['gamma']  # gamma(k): coefficiente di correzione della potenza [°C^-1]
    Tr = gen_foto['Tr']  # T_r: temperatura di riferimento [°C]
    c = gen_foto['c']  # c(k): coefficiente di tilt per il pannello [-]
    r = gen_foto['r']  # r(k): coefficiente correttivo di Ross [°C m^2/W]
    eta = gen_foto['eta']  # eta(k): rendimento complessivo di sistema [-]
    Delta_Tr = gen_foto['Delta_Tr']  # Delta T_r(k): differenza di temperatura standard [°C]
    PN = gen_foto['PN']  # PN(k): potenza nominale [kWe]
    Ta = pd.read_excel('DUMMY_beta_01_eq.xlsx', sheet_name='Meteo', usecols='C', skiprows=2, header=None, nrows=800).values.flatten() #T_a(i): temperatura ambiente [°C]
    DNI = pd.read_excel('DUMMY_beta_01_eq.xlsx', sheet_name='Meteo', usecols='B', skiprows=2, header=None, nrows=800).values.flatten() #Radiazione solare diretta [W/m2]
    #print("Valore della temperatura Ta", Ta)
   
    DNIa = c * DNI
    Tc = Ta + r * DNIa
    Tb = Tc - (DNIa / DNIr) * Delta_Tr
    phi = (DNIa / DNIr) * PN * (1 + gamma * (Tb - Tr)) * eta

    gen_foto['type'] = 'gen_foto'
    gen_foto['G'] = phi * alpha
    #gen_foto['Nabla_G'] = np.diag(phi)
    gen_foto['alpha'] = alpha

    # Uscita W e Nabla_W
    gen_foto['W'] = gen_foto['G']
    #gen_foto['Nabla_W'] = gen_foto['Nabla_G']

    return gen_foto

  