# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 15:06:06 2023

@author: barto
"""
from casadi import *
import pandas as pd
from carica_gen_foto import carica_gen_foto
from carica_gen_comb import carica_gen_comb 
from update_gen_foto import update_gen_foto
from update_gen_comb import update_gen_comb 
from calcola_MPP_Ac_Power import calcola_MPP_Ac_Power
from carica_carico_L1 import carica_carico_L1
from update_gen_comb_ca import update_gen_comb_ca
import matplotlib.pyplot as plt

# Specificare il nome del file Excel e il nome del foglio all'interno del file
excel_filename = 'DUMMY_beta_01_eq.xlsx'
sheet_name1 = 'gen_foto'
sheet_name2 = 'gen_comb'
sheet_name3 = 'carico_L1'
timespan = pd.read_excel(excel_filename, sheet_name='Meteo', usecols='A', skiprows=2, header=None, nrows=800, names=['time'])
N = len(timespan)
# Carica le condizioni iniziali 
x_0 = carica_gen_foto(excel_filename, sheet_name1)
y_0 = carica_gen_comb(excel_filename, sheet_name2)
x_0num = carica_gen_foto(excel_filename, sheet_name1)
y_0num = carica_gen_comb(excel_filename, sheet_name2)
z_0 = carica_carico_L1(excel_filename, sheet_name3)
print("Voce di x_0:", x_0['alpha'])
print("Voce di y_0:", y_0['alpha'])
print("Voce di z_0:", z_0['C'])
print("Dimensioni di x_0 alpha:", x_0['alpha'].shape)
print("Dimensioni di y_0 alpha:", y_0['alpha'].shape)
z_C = z_0['C']
print("valori di z_C:", z_C) # Stampa i valori del carico

#Creazione delle variabili simboliche, occorre specificare le dimensioni delle variabili.
alpha_x_sym = SX.sym("alpha_x_sym",N)
alpha_y_sym = SX.sym("alpha_y_sym",N)

x_updated = update_gen_foto(x_0, alpha_x_sym)
y_updated = update_gen_comb_ca(y_0, alpha_y_sym)
x_G_sym = x_updated['G']
y_G_sym = y_updated['G']

alpha_all =vertcat(alpha_x_sym, alpha_y_sym) #mette in colonna le variabili di ottimizzazione (i setpoint dei dispositivi)

C1=1e-5

f=sum1(y_updated['F'])+C1*sum1(alpha_all) #funzione obiettivo

#aggiungere una voce di dizionario pern tenere conto del contributo alla funzione obiettivo.
g1 = x_G_sym + y_G_sym - z_C #vincolo
g2 = y_updated['Psi'] 
g = vertcat(g1,g2)
#aggiungere una voce di dizionario pern tenere conto del contributo alla funzione obiettivo.
# nlp = {'x': vertcat(alpha_x, alpha_y), 'f': f, 'g': g}
nlp = {'x': vertcat(alpha_x_sym, alpha_y_sym), 'f': f, 'g': g}

# Ricava i valori numerici per le variabili iniziali
x0_values = vertcat(x_0num['alpha'], y_0num['alpha'])


# Creazione del solver
solver = nlpsol('solver', 'ipopt', nlp)
# vincoli variabili ottimizzazione
lbw=DM.zeros(2*N)
ubw=DM.ones(2*N)

# vincoli processo - bilancio energetico
ubg1=DM.zeros(N)
lbg1=DM.zeros(N)
# vincoli processo - accensioni generatore elettrico
ubg2=0
lbg2=-1000

ubgtot = vertcat(ubg1,ubg2)
lbgtot = vertcat(lbg1,lbg2)
# Risolutore del problema
res = solver(x0=x0_values, lbx=lbw, ubx=ubw, lbg=lbgtot, ubg=ubgtot)
# Estrazione dei risultati
alpha_x_result = res['x'][:N].full().flatten()  # Estrai i risultati di alpha_x_sym
alpha_y_result = res['x'][N:].full().flatten()  # Estrai i risultati di alpha_y_sym
#res['x']: Accede al vettore delle variabili ottimizzate dalla soluzione.
#[:N]: Estrae i primi N elementi associati ad alpha_x_sym.
#.full(): è una funzione di CasADi che serve ad ottenere l'equivalente "completo" dell'oggetto simbolico.
# in questo caso converte l'oggetto Casadi in un array NumPy).
#.flatten(): Appiattisce l'array, convertendo una possibile struttura multidimensionale in un array monodimensionale.
x_updated = update_gen_foto(x_0num,alpha_x_result)
y_updated = update_gen_comb(y_0num, alpha_y_result)

timespan['rounded_time'] = timespan['time'].apply(lambda x: (x.hour + x.minute / 60))

# Creazione dei grafici 
plt.figure(1)
plt.plot(timespan['rounded_time'], alpha_x_result, marker='o', linestyle='-', label='alpha_x_sym fotovoltaico')
plt.plot(timespan['rounded_time'], alpha_y_result, marker='s', linestyle='-', label='alpha_y_sym generatore combustibile')
plt.legend()
plt.xlabel('Time(hr)')
plt.ylabel('Alpha Values')

plt.figure(2)
plt.plot(timespan['rounded_time'],x_updated['W'], marker='o', linestyle='-', label='Potenza fotovoltaico')
plt.plot(timespan['rounded_time'],y_updated['W'], marker='s', linestyle='-', label='Potenza generatore combustibile')
plt.plot(timespan['rounded_time'],z_C, marker='s', linestyle='-', label='Potenza carico L1')
plt.legend()
# Aggiunta di etichette agli assi
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')


plt.figure(3)
plt.plot(timespan['rounded_time'],x_updated['W'], marker='o', linestyle='-', label='Potenza fotovoltaico')
plt.plot(timespan['rounded_time'],y_updated['W'], marker='s', linestyle='-', label='Potenza generatore combustibile')
#plt.plot(timespan['rounded_time'],z_C, marker='s', linestyle='-', label='Potenza carico L1')
plt.legend()
# Aggiunta di etichette agli assi
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')


plt.figure(4)
plt.plot(timespan['rounded_time'],x_updated['W'], marker='o', linestyle='-', label='Potenza fotovoltaico')

plt.legend()
# Aggiunta di etichette agli assi
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')


plt.figure(6)
plt.plot(timespan['rounded_time'],y_updated['W'], marker='s', linestyle='-', label='Potenza generatore combustibile')
plt.legend()
# Aggiunta di etichette agli assi
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')

plt.figure(7)
plt.plot(timespan['rounded_time'],x_updated['W']+y_updated['W'], marker='s', linestyle='-', label='Potenza totale generatori')
plt.plot(timespan['rounded_time'],z_C, marker='s', linestyle='-', label='Potenza carico L1')
plt.legend()
# Aggiunta di etichette agli assi
plt.xlabel('Time(hr)')
plt.ylabel('Power[kWe]')


# Aggiunta di una legenda
plt.legend()

# Visualizzazione del grafico
plt.show()

# # Calcolo dei profili x_G_sym e y_G_sym
x_G_result = res['g'][:N].full().flatten()
y_G_result = res['g'][N:].full().flatten()
