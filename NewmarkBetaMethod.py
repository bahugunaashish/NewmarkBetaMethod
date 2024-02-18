# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 17:43:03 2024

@author: Ashish Bahuguna
         Clemson Univerity
         abahugu@clemson.edu

Newmark-Beta Method 
This code is python version of MATLAB code developed by Dr. Pang, 
                                                        Clemson Univeristy 
                                                        Clemson, SC,
                                                        USA

"""
import numpy as np
from ProcessNGA import *
import matplotlib.pyplot as plt
#%%
def NewmarkBetaMethod(mass, stiffness, damping, force, dt, method, flag =0):
    '''
    

    Parameters
    ----------
    mass : TYPE Float
        DESCRIPTION. mass of the structure
        
    stiffness : TYPE float
        DESCRIPTION. stiffness of structure
        
    damping : TYPE float 
        DESCRIPTION. damping of structure calculated as: 
        c = 2 * zeta * (k * m) ** 0.5; where k is stffness, m is mass, zeta is damping ratio
        
    force : TYPE array
        DESCRIPTION. conatins values of force as
        force = mass * accel
        
    dt : TYPE float
        DESCRIPTION. time interval
        
    method : TYPE string
        DESCRIPTION. method of solving as:
            
        Avergae Mehtod  
        gamma = 1/2
        beta = 1/4 
        
        Linear Method
        gamma = 1/2
        beta = 1/6
        


    Returns
    -------
    U : TYPE float array
        DESCRIPTION. Displacment 
        
    V : TYPE float array
        DESCRIPTION. Velocity
        
    A : TYPE float array
        DESCRIPTION. Acceleration
        
    dynStiffness : TYPE float
        DESCRIPTION. dynamic stiffness
        
    a : TYPE float 
        DESCRIPTION. constant value as:
        a = (1 / (beta * dt)) * mass + (gamma / beta) * damping  # table 5.4.2 Eq 1.4
       
    b : TYPE float
        DESCRIPTION. constant as
        b = (1 / (2*beta)) * mass + dt * (gamma / (2 * beta) - 1) * damping
        

    '''
    try:
        nsteps = len(force)
        U = np.array([])  # Displacement
        V = np.array([])  # velocity
        A = np.array([])  # Acceleration
        dforce = np.diff(force)  # difference between consecutive values
    
        if method == 'Average':
            gamma = 1/2
            beta = 1/4
        elif method == 'Linear':
            gamma = 1/2
            beta = 1/6
        else:
            print("Choose from given Options\n 1 - Average \n 2- Linear")
    
        u0 = 0.
        v0 = 0.
        A0 = (force[0] - damping * v0 - stiffness * u0) / mass  # Table 5.4.2, Eq 1.1
    
        U = np.append(U, u0)
        V = np.append(V, v0)
        A = np.append(A, A0)
        # Eq. stiffness
    
        dynStiffness = stiffness + (gamma / (beta * dt)) * damping + mass / (beta * dt ** 2)  # Table 5.4.2, Eq 1.3
    
        # Constant a and b
        # print('=======================================================')
        a = (1 / (beta * dt)) * mass + (gamma / beta) * damping  # table 5.4.2 Eq 1.4
        b = (1 / (2*beta)) * mass + dt * (gamma / (2 * beta) - 1) * damping
        if flag ==0:
            print(f'Constant \na = {a:.3f} \nb = {b:.3f}')
            print('=======================================================')
        # Calculation for each steps
        
            print('============================================================================================================')
            print(f'{"ti":<3}       {"pi":<10}        {"Ai":<10} {"dpi":<8}     {"dp_dyn":<6}   {"dUi":<8}   {"dVi":<8} {"dAi":<8} {"Vi":<10} {"Ui":<10}')
            print('============================================================================================================')
        ti = np.array([])
        for i in range(nsteps-1):
            dP_dyn = dforce[i] + a * V[i] + b * A[i]  # Table 5.4.2 Eq 2.1
    
            dU = dP_dyn/dynStiffness  # Eq 2.2
    
            dV = (gamma / (beta * dt)) * dU - (gamma / beta) * V[i] + dt * (1 - gamma / (2 * beta)) * A[i]
    
            dA = (1. / (beta * dt ** 2)) * dU - (1. / (beta * dt)) * V[i] - (1 / (2 * beta)) * A[i]
            
            
            U = np.append(U, (U[i] + dU))
            V = np.append(V, (V[i] + dV))
            A = np.append(A, (A[i] + dA))
            # print(U[i], U[i+1])
            
            if flag ==0:
                print(f'{i/10:.1f}   {force[i]:>10.4f}   {A[i]:>15.4f}   {dforce[i]:>8.4f}   {dP_dyn:>10.4f} {dU:>8.4f} {dV:>8.4f} {dA:>8.4f} {V[i]:>10.4f} {U[i]:>10.4f}')
            ti = np.append(ti, i*dt)
        return U, V, A, dynStiffness, a, b, ti
    except IOError:
        print("process FAILED!: Provide the method string is case-sensitive")
    
     


