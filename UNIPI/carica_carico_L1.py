# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 02:29:48 2023

@author: barto
"""

import pandas as pd

def carica_carico_L1(excel_filename, sheetname_carico_L1):

    carico_L1 = {}
    carico_L1['sheetname'] = sheetname_carico_L1

    # Carica i dati dalla colonna B di righe 4 a 803
    carico_L1['C'] = pd.read_excel(excel_filename, sheet_name=sheetname_carico_L1, usecols='B', skiprows=3, nrows=800, header=None).values.flatten()
    # Calcola la colonna W come negativo di C
    carico_L1['W'] = -carico_L1['C']

    # Stampa il dizionario carico_L1
    print(carico_L1)

    # Restituisci il dizionario (opzionale)
    return carico_L1
