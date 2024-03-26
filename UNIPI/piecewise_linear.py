# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 16:32:18 2023

@author: barto
"""

import numpy as np

def piecewise_linear(tab_y, tab_x, vec_x):
    vec_y = np.empty_like(vec_x, dtype=float)
    
    for j in range(len(vec_x)):
        x = vec_x[j]
        
        if x <= tab_x[0]:
            y = (tab_y[1] - tab_y[0]) / (tab_x[1] - tab_x[0]) * (x - tab_x[0]) + tab_y[0]
        elif x >= tab_x[-1]:
            y = (tab_y[-1] - tab_y[-2]) / (tab_x[-1] - tab_x[-2]) * (x - tab_x[-1]) + tab_y[-1]
        else:
            for i in range(len(tab_x) - 1):
                if tab_x[i] <= x <= tab_x[i+1]:
                    y = (tab_y[i+1] - tab_y[i]) / (tab_x[i+1] - tab_x[i]) * (x - tab_x[i]) + tab_y[i]
                    break

        vec_y[j] = y

    return vec_y
