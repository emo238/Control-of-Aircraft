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

np.set_printoptions(precision=4,suppress=True)

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
rho = 1.13110       #Density (kg/m^3)
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

#%% mdamp function

def mdamp(A):
    roots=np.linalg.eigvals(A)
    ri=[]
    a=[]
    b=[]
    w=[]
    xi=[]
    st=[]
    for i in range(0,roots.size):
        ri.append(roots[i])
        a.append(roots[i].real)
        b.append(roots[i].imag)
        w.append(math.sqrt(a[i]**2+b[i]**2))
        xi.append(-a[i]/w[i])
        if b[i]>0:
            signb='+'
        else:
            signb='-'
        st.append('%.5f'%(a[i])+signb+'j'+'%.5f'%(math.fabs(b[i]))+'xi='+'%.5f'%(xi[i])+'w='+'%.5f'%(w[i])+'rad/s')
    print("\n",st)
#%% STUDY OF THE UNCONTROLLED AIRCRAFT

#Determine the equilibrium points conditions around the chosen operating point

X_f = -f*l_t
X_fdeltam = -f_delta*l_t
X_G = -c*l_t

X = X_f - X_G
Y = X_fdeltam - X_G

V_eq = M*c_sound
Q = 0.5*rho*V_eq**2

alpha_eq_list = [0]
F_p_eq = [0]

eps = 10**-10
    
i=0

cond = True

while cond :
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
        cond = False
        
    i+=1

print("\nNombre d'itération : ",i)
print("\nalpha_equilibre = ",alpha_eq_list[len(alpha_eq_list)-1])
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
             [-Z_v,0,-Z_alpha,1,0,0],
             [0,0,m_alpha,m_q,0,0],
             [0,0,0,1,0,0],
             [0,V_eq,0,0,0,0]])
Btau=np.matrix([[0,-X_tho],
             [Z_delta_m,Z_tho],
             [-Z_delta_m,-Z_tho],
             [m_delta_m,m_tho],
             [0,0],
             [0,0]])
B=Btau[:,0:1]

C=np.eye(6)
D=np.zeros((6,1))


sys=control.ss(A,B,C,D)
control.damp(sys)

# print("A = ",A)
# print("B = ",B)

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

#%% Mode study

Ar=A[0:4,0:4]
Br=B[0:4,0:1]

# print("\nAr = ",Ar)
# print("\nBr = ",Br)

# Cr=np.matrix([0.0,1.0])
# Dr=0.0

mdamp(Ar)

# print("\n",matrix_rank(control.ctrb(Ar,Br)))

eigenValues, eigenVectors=np.linalg.eig(A)
# print("\nEigenvalues of A :")
# print(eigenValues)
# print("\nEigenvectors of A :")
# print(eigenVectors)

#%% Short period Mode
Ai=Ar[2:4, 2:4]
Bi=Br[2:4, 0]

mdamp(Ai)

Cia=np.matrix([[1, 0]])
Ciq=np.matrix([[0, 1]])
Di=np.matrix([[0]])

TaDm_ss=control.ss(Ai, Bi, Cia, Di)
print("\nTransfer function alpha/delta m : ")
TaDm_tf= control.tf(TaDm_ss)
print(TaDm_tf)
print("\nStatic gain of alpha/delta m =%f "%(control.dcgain(TaDm_tf)))
TqDm_ss= control.ss(Ai, Bi, Ciq, Di)
print("\nTransfer function q/delta_m : " )
TqDm_tf= control.ss2tf(TqDm_ss)
print(TqDm_tf)
print("\nStatic gain of q/delta m =%f "%(dcgain(TqDm_tf)))

figure(2)
Ya,Ta = control.matlab.step(TaDm_tf,arange(0,10,0.01))
Yq,Tq = control.matlab.step(TqDm_tf,arange(0,10,0.01))

plot(Ta,Ya,'b',Tq,Yq,'r',lw=2)
plot([0,Ta[-1]],[Ya[-1],Ya[-1]],'k--',lw=1)
plot([0,Ta[-1]],[0.95*Ya[-1],0.95*Ya[-1]],'k--',lw=1)
plot([0,Tq[-1]],[Yq[-1],Yq[-1]],'k--',lw=1)
plot([0,Tq[-1]],[1.05*Yq[-1],1.05*Yq[-1]],'k--',lw=1)
plot([0,Tq[-1]],[0.95*Yq[-1],0.95*Yq[-1]],'k--',lw=1)

minorticks_on()
grid(b=True, which='both')
#grid(True)
title(r'Step response $\alpha/\delta_m$ et $q/\delta_m$')
legend((r'$\alpha/\delta_m$',r'$q/\delta_m$'))

xlabel('Time(s)')
ylabel(r'$\alpha$ (rad) & $q$ (rad/s)')

Osa, Tra, Tsa = step_info(Ta,Ya)
Osq, Trq, Tsq = step_info(Tq,Yq)
yya = interp1d(Ta, Ya)
plot(Tsa,yya(Tsa),'bs')
text(Tsa,yya(Tsa)-0.2,Tsa)
yyq=interp1d(Tq, Yq)
plot(Tsq,yyq(Tsq),'rs')
text(Tsq,yyq(Tsq)-0.2,Tsq)
print('\nAlpha settling time 5%% =%f s'%Tsa)
print('q Settling time 5%% = %f s'%Tsq)

#%% Phugoid M

Ap = Ar[0:2,0:2]
Bp = Br[0:2,0:1]

mdamp(Ap)

Cpv=np.matrix([[1,0]])
Cpg=np.matrix([[0,1]])
Dp=np.matrix([[0]])

TvDm_ss=control.ss(Ap,Bp,Cpv,Dp)
print("\nTransfer function V/delta_m :")
TvDm_tf=control.tf(TvDm_ss)
print(TvDm_tf)

print("\nStatic gain of alpha/delta_m =%f"%(control.dcgain(TvDm_tf)))

TgDm_ss=control.ss(Ap,Bp,Cpg,Dp)
print("\nTransfer function gamma/delta_m :")
TgDm_tf=control.tf(TgDm_ss)
print(TgDm_tf)

print("\nStatic gain of gamma/delta_m =%f"%(control.dcgain(TgDm_tf)))

figure(3)
Yv,Tv=control.matlab.step(TvDm_tf,arange(0,700,0.1))
Yg,Tg=control.matlab.step(TgDm_tf,arange(0,700,0.1))

plot(Tv,Yv,'b',Tg,Yg,'r',lw=2)
plot([0,Tv[-1]],[Yv[-1],Yv[-1]],'k--',lw=1)
plot([0,Tv[-1]],[1.05*Yv[-1],1.05*Yv[-1]],'k--',lw=1)
plot([0,Tv[-1]],[0.95*Yv[-1],0.95*Yv[-1]],'k--',lw=1)
plot([0,Tg[-1]],[Yv[-1],Yg[-1]],'k--',lw=1)
plot([0,Tg[-1]],[1.05*Yg[-1],1.05*Yg[-1]],'k--',lw=1)
plot([0,Tg[-1]],[0.95*Yg[-1],0.95*Yg[-1]],'k--',lw=1)

minorticks_on()
grid(b=True, which='both')
#grid(True)
title(r'Step response $V/\delta_m$ et $\gamma/\delta_m$')
legend((r'$V/\delta_m$',r'$\gamma/\delta_m$'))

xlabel('Time(s)')
ylabel(r'$V$ (rad) & $\gamma$ (rad/s)')

Osv, Trv, Tsv = step_info(Tv,Yv)
Osg, Trg, Tsg = step_info(Tg,Yg)
yyv = interp1d(Tv, Yv)
plot(Tsv,yyv(Tsv),'bs')
text(Tsv,yyv(Tsv)-0.2,Tsv)
yyg=interp1d(Tg, Yg)
plot(Tsg,yyg(Tsg),'rs')
text(Tsg,yyg(Tsg)-0.2,Tsg)
print('\nAlpha settling time 5%% =%f s'%Tsv)
print('q Settling time 5%% = %f s'%Tsg)













