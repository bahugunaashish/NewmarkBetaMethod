# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 19:37:25 2024

@author: abahugu
"""

import numpy as np
from ProcessNGA import *
import matplotlib.pyplot as plt
from NewmarkBetaMethod import *
from numpy.linalg import inv
from numpy import linalg as LA
from matplotlib.pyplot import figure
np.set_printoptions(precision=3, suppress=True, linewidth=100)
#%% Input 
k = 10.
print(f'Stiffness : {k} kip/in')
m = 0.2533
print(f'mass : {m} kip-sec2/in')
zeta = 0.05
print(f'zeta : {zeta}')

Tn = 2 * np.pi / (k / m) ** 0.5
print(f'Period : {Tn:.2f} sec')
wn = (k / m) ** 0.5

c = 2 * zeta * (k * m) ** 0.5
print(f'damping : {c:.2f} kip-s/in')

dt = 0.1
print(f'dt : {dt} sec')
t = np.linspace(0., 1., num=11)
p = 10 * np.sin(np.pi * t / 0.6)

p[t >= 0.6] = 0

plt.plot (t, p, color = 'blue', label = 'P(t)')
plt.xlabel('Time (s)')
plt.ylabel('P(t) kip')
plt.grid(color = 'grey', linestyle = 'dashed', linewidth = 0.5)
plt.xlim([0., abs(t[-1])])
plt.ylim([0., max(p)])
plt.legend()
plt.show()
# print(t)
# print(p)
method = 'Average'
# print(t, p)
print('<<<<<<<<<<<<<<<<<Verification of python code of NewMark Beta Method >>>>>>>>>>>>>>>>>>>>> ')

U, V, A, dynStiffness, a, b, ti = NewmarkBetaMethod(m, k, c, p, dt, method, flag=0)

print('=============================================================================================================')
print('****** all Values are matching with AKC Table E5.3******\n\n\n')
