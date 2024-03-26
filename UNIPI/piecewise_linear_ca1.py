# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 11:17:03 2023

@author: barto
"""
from casadi import*
import numpy as np
import casadi as ca
def piecewise_linear_ca1(tab_y, tab_x, vec_x):   
                                    
    vec_y = SX.zeros(vec_x.size())    #Viene inizializzato un vettore simbolico vec_y della stessa dimensione di vec_x 
                                      #(quindi con lo stesso numero di elementi) e tutti i suoi elementi sono inizializzati a zero.

    x = MX.sym('x')  # Viene creata una variabile simbolica x utilizzando la libreria CasADi. 
                     #Questa variabile simbolica rappresenta la variabile indipendente della funzione.
    for j in range(vec_x.numel()):               #for j in range(vec_x.size()[0]):
        x_val = vec_x[j]    # Estrae il valore corrente da vec_x all'indice j.
         
       
        # Questa parte di codice gestisce i casi in cui x è al di fuori dei limiti di tab_x
        y = if_else(x <= tab_x[0],
                       (tab_y[1] - tab_y[0]) / (tab_x[1] - tab_x[0]) * (x - tab_x[0]) + tab_y[0],
                       if_else(x >= tab_x[-1],     #L'utilizzo di -1 nell'indicizzazione tab_x[-1] si riferisce all'ultimo elemento nella lista tab_x. 
                                  (tab_y[-1] - tab_y[-2]) / (tab_x[-1] - tab_x[-2]) * (x - tab_x[-1]) + tab_y[-1],
                                  0  # Valore fittizio, verrà sovrascritto nel ciclo sottostante
                                  ))

        # Questo ciclo annidato gestisce i segmenti tra i punti di tabulazione. 
        # Utilizza la funzione if_else per definire la funzione a tratti tra i punti di tabulazione.
        for i in range(len(tab_x) - 1):  
            idx_i = ca.logic_and ((tab_x[i] <= x), (x <= tab_x[i+1]))
            y = if_else((idx_i),
                           (tab_y[i+1] - tab_y[i]) / (tab_x[i+1] - tab_x[i]) * (x - tab_x[i]) + tab_y[i],
                           y
                           )
        # La seguente riga di codice crea una Funzione CasADi per valutare l'espressione. 
        piecewise_func = Function('piecewise_func', [x], [y]) #x è la variabile simbolica che rappresenta l'input della funzione a tratti.
                                #y è l'espressione che definisce la funzione a tratti in base a come è stata calcolata nel ciclo precedente
                               
        vec_y[j] = piecewise_func(x_val)

    return vec_y

