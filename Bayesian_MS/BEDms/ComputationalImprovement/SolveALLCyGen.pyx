# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 18:10:35 2020

@author: s1778490
"""

import numpy as np
import Model1_Cython as M1
import Model2_Cython as M2
import Model3_Cython as M3

def solveALLCy2 (ts, pD, sp, inputs, ivss, pre, model, realt):
    
    AllSol = np.empty((len(ts),len(pD[:,1])*2))
    AllSolTest = np.empty((len(realt),len(pD[:,1])*2))
    
    for drawInd in range(0,len(pD[:,1])):
        
        p = pD[drawInd,:]
        
        if model == 'M1':
            temp = M1.solve_coupled_ode1(ts, p, sp, inputs, ivss, pre)
        elif model == 'M2':
            temp = M2.solve_coupled_ode2(ts, p, sp, inputs, ivss, pre)
        elif model == 'M3':
            temp = M3.solve_coupled_ode3(ts, p, sp, inputs, ivss, pre)
        else:
            print('Please, select a correct model as M1, M2 or M3 for the desired one')

        AllSolTest[:,drawInd] = temp[:,2][realt]
        AllSolTest[:,drawInd+(len(pD[:,1]))] = temp[:,3][realt]

    return(AllSolTest)
    