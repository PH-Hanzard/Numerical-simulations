# -*- coding: utf-8 -*-
from __future__ import division #Division retourne floating point number
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import scipy.constants as cte


def GenerPorte(centre,largeur_porte, span, pts):
    xi = centre - (span/2)
    xf = centre + (span/2)     
    lambd = np.linspace(xi,xf,pts)
    I=[]
    for i in range(lambd.size):
        if i < (((centre-(largeur_porte/2))-xi)/(xf-xi))*lambd.size:
            I.append(0.)
        elif i > (((centre+(largeur_porte/2))-xi)/(xf-xi))*lambd.size:
            I.append(0.)
        else:
            I.append(1)
    return lambd,np.asarray(I)
    
def Plot(x,y,couleur,titre,xlabel,ylabel):
    sns.set_context("talk")
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(x,y, label='',marker='',color=couleur)
    plt.title(titre)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    
def LambdaToOmega(Lambd):
    return 2 * np.pi * cte.c / Lambd
    
def OmegaToLambda(Omega):
    return 2 * np.pi * cte.c / Omega
    
def OmegaTild(Omega,Omega_centrale):
    return Omega - Omega_centrale
    
def LambdToK(Lambd):
    return 2 * np.pi / Lambd
    
def TransferLens(f_c,x,f,W):
    return np.exp(-2 * (f_c**2) * (x**2) / ((f**2) * (W**2)))
    
def TransferVipa(r1,r2,k,t,n_r,theta_in,theta_i,x,f):
    Delta_1 = 2 * t * n_r * np.cos(theta_in)
    Delta_2 = 2 * t * np.tan(theta_in) * np.cos(theta_i) * x / f
    Delta_3 = t * np.cos(theta_in) * (x**2) / (n_r * (f**2))
    Delta = Delta_1 - Delta_2 - Delta_3
    
    Res_1 = (1-(r1*r2))**2
    Res_2 = 4 * r1 * r2 * (np.sin(k * Delta / 2))**2
    
    Res = 1 / (Res_1 + Res_2)
    return Res
    
def TransferGrating(y,lambd,f,d,theta_d,omega_tild,theta_ig,W):
    alpha = ((lambd**2)*f) / (2 * np.pi * cte.c * d * np.cos(theta_d))
    w_0 = (np.cos(theta_ig) * f * lambd) / (np.cos(theta_d) * np.pi * W)
    
    Res_1 = (y - (alpha * omega_tild))**2
    
    Res = np.exp(-1*Res_1/(w_0**2))
    return Res


def theta_inn(theta_i,n_r):
    return theta_i/n_r
    
    
def Normalize(fct):
    return fct/max(fct)
    
def ToLog(x):
    return 10*np.log10(x)
#Attention : k , lambd, tehta_inn
#==============================================================================
# VARIABLES
#==============================================================================
#Generation porte
Lambda_centrale = 1550e-9
Delta_lambda = 30e-9
Span_lambda = 50e-9
NbPts = 5e5

#Fct transfert lentille
f_c = 100e-3
x = -1.5e-3
f = 180e-3
W = 1.2e-3

#Fct transfert VIPA
r1 = 1.
r2 = 0.97
n_r = 2.
t = 1.3e-3
theta_i = math.radians(4.2)
theta_in = theta_inn(theta_i,n_r)

#FCT transfert grating
y = 0
d = 1100*1e3 
theta_d = math.radians(50)
theta_ig = math.radians(70)



#==============================================================================
# EXECUTION
#==============================================================================


#cree porte en longueur d onde
Lambd,intensity = GenerPorte(Lambda_centrale,Delta_lambda, Span_lambda, NbPts)

#On passe en omega tild
OmegaTil = OmegaTild(LambdaToOmega(Lambd), LambdaToOmega(Lambda_centrale))



I_in = intensity
TransfertLens = TransferLens(f_c,x,f,W)
TransfertVipa = TransferVipa(r1,r2,LambdToK(Lambd),t,n_r,theta_in,theta_i,x,f)
TransfertGrating = TransferGrating(y,Lambd,f,d,theta_d,OmegaTil,theta_ig,W)

I_out = I_in * TransfertLens * TransfertVipa * TransfertGrating
plt.plot(Lambd*1e9,(ToLog(I_out)), label='',marker='',color='k')
    
print('x = ',x)
for i in range (30):
    x += 0.1e-3
    y += 0
    print('x = ',x)
    I_in = intensity
    TransfertLens = TransferLens(f_c,x,f,W)
    TransfertVipa = TransferVipa(r1,r2,LambdToK(Lambd),t,n_r,theta_in,theta_i,x,f)
    TransfertGrating = TransferGrating(y,Lambd,f,d,theta_d,OmegaTil,theta_ig,W)

    I_out = I_in * TransfertLens * TransfertVipa * TransfertGrating
    plt.plot(Lambd*1e9,ToLog(I_out), label='',marker='')
    


sns.set_context("talk")
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
#plt.plot(Lambd,I_out, label='',marker='',color='r')
plt.title(u'Transmission du Système Réseau + VIPA')
plt.xlabel(u'Longueur d onde (nm)')
plt.ylabel(u'Intensite (échelle log)')
plt.show()


