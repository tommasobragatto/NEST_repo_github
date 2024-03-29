# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:43:28 2023

@author: barto
"""
import pandas as pd
import numpy as np
import casadi as ca
from piecewise_linear_ca1 import piecewise_linear_ca1
#from carica_gen_foto_SANDIA import carica_gen_foto_SANDIA
def update_gen_comb_ca(gen_comb, alpha):

    
    gen_comb['alpha'] = alpha
    
    # Carico i parametri di ingresso
    PN = gen_comb['PN']  # PN(k) : potenza nominale [kWe]
    P_min = gen_comb['P_min']  # P_{min}(k) : potenza minima tecnica [kWe]
    N_acc_max = gen_comb['N_acc_max']  # M_{max}(k): numero massimo di accensioni [-]
    
    # Calcolo G e F
    alpha_min = P_min / PN
    gen_comb['F'], gen_comb['G'], F_sig, G_sig, gen_comb['theta'] = calcola_FG(gen_comb)
    
    # Calcolo Psi 
    Psi, Psi_sig, N_acc, N_acc_sig = calcola_Psi(alpha, alpha_min, N_acc_max)
    gen_comb['Psi'] = Psi
    gen_comb['N_acc'] = N_acc

    # Uscita W
    gen_comb['W'] = gen_comb['G']
    
    return gen_comb

   

def sigmoid(alpha, alpha_min):
    theta = (alpha >= alpha_min)  # theta esatto (discontinuo)
    a = 20
    theta_sig = 1 / (1 + np.exp(-a * (alpha - alpha_min)))  # theta approssimato (continuo)
    return theta, theta_sig


def calcola_Psi(alpha, alpha_min, N_acc_max):
    #  funzione sigmoid 
    theta, theta_sig = sigmoid(alpha, alpha_min)

    # Converti gli oggetti CasADi in array NumPy
    theta_np = ca.vertcat(theta)
    
    # Usa la funzione di CasADi per calcolare la differenza
    Delta_theta = ca.diff(theta_np)
    
    # Calcola il valore di Psi sommando i valori assoluti delle differenze in theta
    Psi = (ca.sum1(ca.fabs(Delta_theta))) / (2 * N_acc_max) - 1

    # Converti anche theta_sig in array NumPy se necessario
    theta_sig_np = ca.vertcat(theta_sig)

    # Usa la funzione di CasADi per calcolare la differenza in theta_sig
    Delta_theta_sig = ca.diff(theta_sig_np)

    # Calcola il valore di Psi_sig sommando i valori assoluti delle differenze in theta_sig
    Psi_sig = (ca.sum1(ca.fabs(Delta_theta_sig))) / (2 * N_acc_max) - 1

    # Conta il numero di valori positivi in Delta_theta, rappresenta il numero di incrementi
    N_acc = ca.sum1(Delta_theta > 0)

    # Conta il numero di valori positivi in Delta_theta_sig, rappresenta il numero di incrementi in theta_sig
    N_acc_sig = ca.sum1(Delta_theta_sig > 0)

    # Restituisce i valori calcolati
    return Psi, Psi_sig, N_acc, N_acc_sig



def calcola_FG(gen_comb):
    alpha = gen_comb['alpha']
    PN = gen_comb['PN']  # PN(k) : potenza nominale [kWe];
    P_min = gen_comb['P_min']  # P_{min}(k) : potenza minima tecnica [kWe];
    PCI = gen_comb['PCI']  # PCI(k): potere calorifico inferiore del combustibile [kJ/kg]
    eta_values = gen_comb['eta_values']  # Curva caratteristica rendimento elettrico della macchina \eta_e(\alpha) [-];
    
    alpha_min = P_min / PN
    theta, theta_sig = sigmoid(alpha, alpha_min)
    
    phi1 = PN
    G = phi1 * theta * alpha
    G_sig = phi1 * theta_sig * alpha
    
    phi2 = PN / (1000 * PCI) * 3600
    
   
  

    # Convertiamo la lista eta_values in un array NumPy
    eta_values_array = np.array(eta_values)

    # Ora è possibile utilizzare la sintassi di slicing
    eta = piecewise_linear_ca1(eta_values_array[:, 2], eta_values_array[:, 1], alpha)

    # eta = piecewise_linear(eta_values[:, 2], eta_values[:, 1], alpha)
   
    eta = eta / 100  # In tabella e' in percentuale ma qui serve in [0,1].
    
    F = phi2 * theta * alpha / eta  # Consumo combustibile in [kg/h]
    F_sig = phi2 * theta_sig * alpha / eta
    
    return F, G, F_sig, G_sig, theta









