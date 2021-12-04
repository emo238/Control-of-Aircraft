# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 09:00:14 2021

@author: Emeric Villette
"""
#%% INITIALISATION

from __future__ import unicode_literals
from matplotlib.pyplot import *
import control
from control.matlab import *
from math import *
from scipy.interpolate import interp1d
from pylab import *
from matplotlib.widgets import Slider
import numpy as np
import scipy.interpolate
from sisopy31 import*

#%% AIRCRAFT CHARACTERISTIC

l_ref = 5.24        #Reference length
l_t = (3/2)*l_ref   #Total length
m = 8400            #Mass in kg
c = 52/100          #Aircraft centering (center of gravity position)(as % of total length)
S = 34              #Reference surface (Wings) in m^2
r_g = 2.65          #Radius of gyration in m
g = 9.81            #Gravity

M = 1.19            #Mach
h = 2700            #Altitude (ft)

T = 282.801         #Temperature in Kelvin
P = 91821.4         #Pressure Pa
rho = 1.13110         #Density (kg/m^3)
c_sound = 337.121   #Speed of sound (m/s)


#%% GRAPH VALUE FOR M = 1.19

C_x0 = 0.0375       #Drag coeficient for null incidence
C_z_alpha = 2.85    #Lift gradient coeficient wrt alpha
C_z_deltam = 0.9    #Lift gradient coeficient wrt delta m
delta_m0 = 0         #Equilibrium fin deflection for null lift
alpha_0 = 0.006     #Incidence for Null lift and  Null fin deflection
f = 0.608           #Aerodynamic center of body and wings
f_delta = 0.9       #Aerodynamic center of fins
k = 0.34            #Polar coeficient
C_mq = -0.42        #Damping coeficient 

#%% STUDY OF THE UNCONTROLLED AIRCRAFT

#Determine the equilibrium points conditions around the chosen operating point

X_f = f*l_t
X_fdeltam = f_delta*l_t
X_G = c*l_t

X = X_f - X_G
Y = X_fdeltam - X_G

V_eq = M*c_sound
Q = 0.5*rho*V_eq**2

alpha_eq_list = [0]
F_p_eq = [0]

eps = 10**-10
    
i=0

prout = True

while prout :
    C_z_eq = (1/(Q*S))*(m*g-F_p_eq[i]*np.sin(alpha_eq_list[i]))
    C_x_eq = C_x0 + k*(C_z_eq**2)
    C_x_deltam = 2*k*C_z_eq*C_z_deltam
    num_1 = C_x_eq*np.sin(alpha_eq_list[i]) + C_z_eq*np.cos(alpha_eq_list[i])
    den_1 = C_x_deltam*np.sin(alpha_eq_list[i])+C_z_deltam*np.cos(alpha_eq_list[i])
    delta_m_eq = delta_m0 - (num_1/den_1)*(X/(Y-X))
    
    a_l_p_h_a = alpha_0 + C_z_eq/C_z_alpha - C_z_deltam*delta_m_eq/C_z_alpha
    F_equilibre = (Q*S*C_x_eq)/(np.cos(a_l_p_h_a))
    
    alpha_eq_list.append(a_l_p_h_a)
    F_p_eq.append(F_equilibre)
    
    if (abs(a_l_p_h_a-alpha_eq_list[i])<eps):
        prout = False
        
    i+=1

print("Nombre d'itération : ",i)
print("alpha_equilibre = ",alpha_eq_list[len(alpha_eq_list)-1])
alpha_eq=alpha_eq_list[len(alpha_eq_list)-1]

#Build a small signals model (Matrix construction)

F_eq = (0.5*rho*V_eq**2*S*C_x_eq)/(np.cos(alpha_eq))
C_x_alpha = 2*k*C_z_eq*C_z_alpha
gamma_eq=0
F_tho = 0
C_m_alpha = (X/l_ref)*(C_x_alpha*np.sin(alpha_eq)+C_z_alpha*np.cos(alpha_eq))
I_yy=m*r_g**2
C_m_delta_m = (Y/l_ref)*(C_x_deltam*np.sin(alpha_eq)+C_z_deltam*np.cos(alpha_eq))

X_v = (2*Q*S*C_x_eq)/(m*V_eq)
X_alpha = F_eq/(m*V_eq)*np.sin(alpha_eq)+(Q*S*C_x_alpha)/(m*V_eq)
X_gamma = (g*np.cos(gamma_eq))/V_eq
X_deltam = (Q*S*C_x_deltam)/(m*V_eq)
X_tho = (-F_tho*np.cos(alpha_eq))/(m*V_eq)

m_v = 0
m_alpha = (Q*S*l_ref*C_m_alpha)/I_yy
m_q = (Q*S*l_ref**2*C_mq)/(V_eq*I_yy)
m_delta_m = (Q*S*l_ref*C_m_delta_m)
m_tho = 0

Z_v = (2*Q*S*C_z_eq)/(m*V_eq)
Z_alpha = (F_eq*np.cos(alpha_eq))/(m*V_eq)+(Q*S*C_z_alpha)/(m*V_eq)
Z_gamma = (g*np.sin(gamma_eq))/V_eq
Z_delta_m = (Q*S*C_z_deltam)/(m*V_eq)
Z_tho = (F_tho*np.sin(alpha_eq))/(m*V_eq)

A=np.matrix([[-X_v,-X_gamma,-X_alpha,0,0,0],
             [Z_v,0,Z_alpha,0,0,0],
             [-Z_v,0,Z_alpha,0,0,0],
             [0,0,m_alpha,m_q,0,0],
             [0,0,0,1,0,0],
             [0,V_eq,0,0,0,0]])
B=np.matrix([[0,-X_tho],
             [Z_delta_m,Z_tho],
             [-Z_delta_m,-Z_tho],
             [m_delta_m,m_tho],
             [0,0],
             [0,0]])

#print("A = ",A)
#print("B = ",B)

# submatrix ( matr i x s l i c i n g )
# note t h a t ind i ce s begin at 0 and
# the i n d i c e of the end i s not inc luded
Ar=A[2:4,2:4]
Br=B[2:4,0:1]
Cr=np.matrix([0.0,1.0])
Dr=0.0
#sss space state system
sys_q=ss(Ar,Br,Cr,Dr)
#transfer function
sys_q_tf=tf(sys_q)
#damping ratio
control.matlab.damp(sys_q)
#dcgain
dcgain(sys_q)
#step_response
figure(1)
Yq,Tq=control.matlab.step(sys_q)
plot(Tq,Yq,'b',lw=2)
grid(True)
title('Réponse indicielle q/Dm')
show()
#feedback
Kq=-0.115
Tqbo=feedback(Kq*sys_q,1)
sisotool(sys_q)
