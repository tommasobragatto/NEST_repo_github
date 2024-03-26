# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 15:06:06 2023

@author: barto
"""
from casadi import *
import casadi as ca
import pandas as pd
from carica_gen_foto import carica_gen_foto
from update_gen_foto import update_gen_foto
from carica_gen_comb import carica_gen_comb 
from update_gen_comb import update_gen_comb 
from update_gen_comb_ca import update_gen_comb_ca
from calcola_MPP_Ac_Power import calcola_MPP_Ac_Power
from carica_gen_eoli import carica_gen_eoli 
from update_gen_eoli import update_gen_eoli 
from update_gen_eoli_ca1 import update_gen_eoli_ca1 
from carica_carico_L1 import carica_carico_L1
import matplotlib.pyplot as plt
from openpyxl import load_workbook
from openpyxl.utils.dataframe import dataframe_to_rows

# Specificare il nome del file Excel e il nome del foglio all'interno del file
excel_filename = 'DUMMY_beta_01_eq.xlsx'
sheet_name1 = 'gen_foto'
sheet_name3 = 'carico_L1'
sheet_name4 ='Tariffa'
sheet_name6 = 'gen_eoli'
timespan = pd.read_excel(excel_filename, sheet_name='Meteo', usecols='A', skiprows=2, header=None, nrows=800, names=['time'])
N = len(timespan)    #ORIZZONTE DI SIMULAZIONE E NUMERO VARIABILI DI OTTIMIZZAZIONE PER CIASCUN DISPOSITIVO, N=96

#------------------------ INIZIALIZZAZIONE -----------------------------------------------------------------------------------------------------------------------------------

# Carica le condizioni iniziali per ogni dispositivo
x_0 = carica_gen_foto(excel_filename, sheet_name1)
x_0num = carica_gen_foto(excel_filename, sheet_name1) 
z_0 = carica_carico_L1(excel_filename, sheet_name3)
e_0 = carica_gen_eoli(excel_filename, sheet_name6)
e_0num = carica_gen_eoli(excel_filename, sheet_name6)
z_C = z_0['C']
#Stampare i valori iniziali dei set-point di ogni dispositivo: 
print("Voce di x_0:", x_0['alpha'])      # generatore fotoltaico
print("Voce di e_0:", e_0['alpha'])      # generatore eolico
print("valori di z_C:", z_C) 


#------------------------ CREAZIONE VARIABILI SIMBOLICHE -----------------------------------------------------------------------------------------------------------------------------------

#Creazione delle variabili simboliche, occorre specificare le dimensioni delle variabili ovvero N.
alpha_x_sym = SX.sym("alpha_x_sym",N)      #var.sim. per generatore fotovoltaico
#alpha_y_sym = SX.sym("alpha_y_sym",N)      #var.sim. per generatore a combustibile
alpha_e_sym = SX.sym("alpha_e_sym",N)      #var.sim. per generatore eolico
#beta_b_sym = SX.sym("beta_b_sym",N)        #var.sim. per accumulatore elettrico_BMS
alpha_all =vertcat(alpha_x_sym, alpha_e_sym) 
#------------------------ AGGIORNAMENTO -----------------------------------------------------------------------------------------------------------------------------------

#Aggiornamento degli oggetti che diventano simbolici
x_updated = update_gen_foto(x_0, alpha_x_sym)        #aggiornamento per generatore fotovoltaico_Sandia
#y_updated = update_gen_comb_ca(y_0, alpha_y_sym)            #aggiornamento per generatore a combustibile
e_updated = update_gen_eoli_ca1(e_0, alpha_e_sym)            #aggiornamento per generatore eolico
#b_updated = update_accumulatore_BMS_ca(b_0, beta_b_sym)     #aggiornamento per accumulatore elettrico_BMS


#------------------------ IMPOSTAZIONE PROBLEMA -----------------------------------------------------------------------------------------------------------------------------------

#Potenze dei generatori, dell'accumulatore e del carico che diventano simbolici
x_G_sym = x_updated['W']    #voce dizionario potenza generatore fotovoltaico_SANDIA->simbolica
e_G_sym = e_updated['W']    #voce dizionario potenza generatore eolico->simbolica
Z_C_sym=  z_0['W']          #voce dizionario potenza carico L1->simbolica


#CALCOLO DELLA POTENZA EROGATA complessivamente dai dispositivi
W_sym=x_G_sym+e_G_sym-z_C

#LETTURA E ASSEGNAZIONE DEI VALORI PREZZO DURANTE LA GIORNATA ALLA VARIABILE df_tariffa
df_tariffa = pd.read_excel(excel_filename, sheet_name4, skiprows=3, header=None, nrows=800, usecols=[11, 12])
print(len(df_tariffa))

# Comando che utilizza casadi.if_else per selezionare i valori di C_wf in base al segno di W
# dove cWf è uguale a “-prezzo di vendita” per gli istanti in cui W>0 e a “-prezzo di acquisto” per gli istanti in cui W<0
C_wf = ca.if_else(W_sym > 0, -df_tariffa.iloc[:, 0].values, ca.if_else(W_sym == 0, 0, -df_tariffa.iloc[:, 1].values))
C1=1e-5


# Minimizzazione della funzione costo f:
f= sum1(W_sym * C_wf) + C1*sum1(alpha_all)

#------------------------ VINCOLI -----------------------------------------------------------------------------------------------------------------------------------

#VINCOLI VARIABILI OTTIMIZZAZIONE

#vincoli inferiori variabili di ottimizzazione
lbxgf=DM.zeros(N)         #vincolo inferiore generatore fotovoltaico
lbxeo=DM.zeros(N)          #vincolo inferiore generatore eolico

#vincoli inferiori variabili di ottimizzazione TOTALI
lbxtot=vertcat(lbxgf,lbxeo)
print("vincoli inferiori variabili:", lbxtot.shape) #STAMPA LE DIMENSIONI DEI VINCOLI TOTALI DELLE VARIABILI INFERIORI

#vincoli superiori variabili di ottimizzazione
ubxgf=DM.ones(N)          #vincolo superiore generatore fotovoltaico
ubxeo=DM.ones(N)           #vincolo superiore generatore eolico
 

#vincoli superiori variabili di ottimizzazione TOTALI
ubxtot=vertcat(ubxgf,ubxeo)
print("vincoli superiori variabili:", ubxtot.shape) #STAMPA LE DIMENSIONI DEI VINCOLI TOTALI DELLE VARIABILI SUPERIORI

#------------------------ RISOLUTORE -----------------------------------------------------------------------------------------------------------------------------------

nlp = {'x': vertcat(alpha_x_sym, alpha_e_sym), 'f': f}

# Ricava i valori numerici per le variabili iniziali
x0_values = vertcat(x_0num['alpha'], e_0num['alpha'])

# Creazione del solver
solver = nlpsol('solver', 'ipopt', nlp)

# Risoluzione del problema
res = solver(x0=x0_values, lbx=lbxtot, ubx=ubxtot)
# Estrazione dei risultati

# Estrazione dei risultati, e cioè delle variabili di ottimizzazione per ciascun dispositivo
alpha_x_result = res['x'][:N].full().flatten()  # Estrae i risultati di alpha_x_sym, variabile ottimizzazione generatore fotovoltaico_Sandia 
alpha_e_result = res['x'][N:2*N].full().flatten()  # Estrai i risultati di alpha_y_sym, variabile ottimizzazione generatore combustibile

#Gli oggetti vengono aggiornati, le variabili di ottimizzazione determinate vengono inserite come argomento.
x_updated = update_gen_foto(x_0num,alpha_x_result)
e_updated = update_gen_eoli(e_0num, alpha_e_result)
timespan['rounded_time'] = timespan['time'].apply(lambda x: (x.hour + x.minute / 60))
#------------------------ GRAFICI -----------------------------------------------------------------------------------------------------------------------------------

# GRAFICO PROFILO DEI SET-POINT DEI DISPOSITIVI
plt.figure(1)
plt.plot(timespan['rounded_time'], alpha_x_result, marker='o', linestyle='-', label='alpha_x_sym')
plt.plot(timespan['rounded_time'], alpha_e_result, marker='h', linestyle='-', label='alpha_e_sym')
plt.legend()

# Aggiunta di etichette agli assi del primo grafico
plt.xlabel('Time(hr)')
plt.ylabel('Alpha Values')
plt.title('Device set-point profile as a function of time')

# GRAFICO PROFILO DELLE POTENZE EROGATE
plt.figure(2)
plt.plot(timespan['rounded_time'],x_updated['W'], marker='o', linestyle='-', label='Potenza fotovoltaico')
plt.plot(timespan['rounded_time'],e_updated['W'], marker='h', linestyle='-', label='Potenza eolico')
plt.plot(timespan['rounded_time'],z_C, marker='s', linestyle='-', label='Potenza carico L1')
plt.legend()

# Aggiunta di etichette agli assi del secondo grafico
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')
plt.title('Power profile as a function of time')


plt.figure(3)
plt.plot(timespan['rounded_time'],x_updated['W'], marker='o', linestyle='-', label='Potenza fotovoltaico')
plt.plot(timespan['rounded_time'],z_C, marker='s', linestyle='-', label='Potenza carico L1')
plt.legend()

# Aggiunta di etichette agli assi del secondo grafico
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')
plt.title('PV__GEN Power profile as a function of time')



plt.figure(4)
plt.plot(timespan['rounded_time'],e_updated['W'], marker='h', linestyle='-', label='Potenza eolico')
plt.plot(timespan['rounded_time'],z_C, marker='s', linestyle='-', label='Potenza carico L1')
plt.legend()

# Aggiunta di etichette agli assi del secondo grafico
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')
plt.title('Wind_GEN Power profile as a function of time')


plt.figure(5)
plt.plot(timespan['rounded_time'],x_updated['W']+e_updated['W'], marker='s', linestyle='-', label='potenza totale generatori')
plt.plot(timespan['rounded_time'],z_C, marker='s', linestyle='-', label='Potenza carico L1')
plt.legend()
# Aggiunta di etichette agli assi
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')








