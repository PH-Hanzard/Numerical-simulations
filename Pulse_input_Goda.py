# -*-coding:utf-8 -*
from __future__ import division #Division retourne floating point number
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.constants
import scipy.fftpack
import seaborn as sns



#==============================================================================
# Recupere spectre longueur onde - convertit frequence - applique phase quadratique - ifft
# Peut afficher la forme temporelle a un point ou generer plusieurs images
#==============================================================================
data = np.genfromtxt('Goda_freq.dat', delimiter=';', skip_header=0, skip_footer=0, names=['x', 'y'])
#data['y']=data['y']/(max(data['y']))


#==============================================================================
# Pour spectre du spectro
#==============================================================================

#Calcul FWHM SPECTRE
#for k in xrange(0,len(data['y'])):
#    if data['y'][k]>0.4 and data['y'][k]<0.6  :
#        fwhm1_lambda=k
#        break
#for k in xrange(len(data['y'])-1,0,-1):
#    if data['y'][k]>0.4 and data['y'][k]<0.6  :
#        fwhm2_lambda=k
#        break
#FWHM_lambda = data['x'][fwhm2_lambda]-data['x'][fwhm1_lambda]
#print('FWHM lambda = ',FWHM_lambda,' nm')


frequency=scipy.constants.c/(data['x']*1e-9)
Intensite_frequence=data['y']


#116ps2/kmé violette - 153 PRA
beta_2=153e-24
beta_3=0.24e-36
w=2*np.pi*frequency
# Rezki : 192.853e12 -- PRA : 193.557

for i in range(3):
    w_0=2*np.pi*193.557e12
    z=i*2000
    beta_2_z=beta_2*z*1e-3
    beta_3_z=beta_3*z*1e-3
    #+ regler longueur d onde centrale
    
    #Calcul FWHM Therorique
    #Lambda_centre=0.5*(data['x'][fwhm2_lambda]+data['x'][fwhm1_lambda])*1e-9
    #Dispersion=(2*np.pi*scipy.constants.c/(Lambda_centre**2))*beta_2*1e-3
    #FWHM_theorique=Dispersion*z*FWHM_lambda*1e-9
    #print('FWHM Temporelle theorique : ',FWHM_theorique*1e9,' ns')
    
    def Qua(beta_2_z,beta_3_z,w):
            return np.exp((-1j*beta_2_z*((w-w_0)**2))*0.5)*np.exp((-1j*beta_3_z*((w-w_0)**3))/6.)  
            
    #    return np.exp((-1j*beta_2_z*(w-w_0)**2)/2) 
    I_Qua=np.sqrt(Intensite_frequence)*Qua(beta_2_z,beta_3_z,w)
    T_Qua = (scipy.ifft(I_Qua))
    TPS_Qua = scipy.fftpack.fftfreq(T_Qua.size,frequency[1]-frequency[0])
    
    T_Qua_final=((abs(T_Qua))**2)/(max(abs(T_Qua))**2)
#    T_Qua_final = TPS_Qua
    
    plt.figure(1)
    #plt.subplot(223)
    plt.ylabel('Intensite')
    plt.xlabel('temps (ns)')
    xlim=-1
    x2lim=1
    plt.axis([xlim, x2lim, 0, max(T_Qua_final)+0.1])
    
    plt.plot((TPS_Qua)*1e9,T_Qua_final,marker='',label='%i m'%z)






#plt.legend(bbox_to_anchor=(-0.55, 1), loc=2, borderaxespad=0.)
plt.show()
#plt.subplot(221)
#plt.plot(data['x'], data['y'], label='', color='r')
##plt.plot(xvals,yinterp, label='', color='b',marker='x')
#
#plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
#plt.title('Input spectrum')
#plt.xlabel("Longueur d'onde (nm)")
#plt.ylabel("Intensite (u.a.)")
#plt.axis([1548.4, 1549.6, 0, max(data['y']+(0.1*data['y']) )])
#
#plt.subplot(222)
#plt.plot(frequency*1e-12, (abs(Intensite_frequence))**1, label='',color='g')
#plt.title('Frequency spectrum')
#plt.xlabel("Frequency (THz)")
#plt.ylabel("Intensite (u.a.)")
#plt.axis([193.45, 193.62, 0, max(data['y']+(0.1*data['y']) )])
#
#plt.subplot(224)
#plt.ylabel('Intensite')
#plt.xlabel('Frequence (THz)')
#plt.axis([193.45, 193.62, 0, max(data['y']+(0.1*data['y']) )])
#plt.plot(frequency*1e-12,(abs(I_Qua))**2,color='g')
#mng=plt.get_current_fig_manager()
#mng.window.showMaximized()
##plt.style.use(‘ggplot’)
#plt.show()