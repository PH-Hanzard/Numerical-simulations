from __future__ import division #Division retourne floating point number
import numpy as np
#import math
import matplotlib.pyplot as plt
import scipy
import scipy.constants

#==============================================================================
# Cree spectre frequence - convertit en longueur onde pour fiberdesk - cree fichier
# On appliquera ensuite FFT. Donc remember : span temporel inverse a span frequence. Pour une bonne resolution de 
# l ipulsion de sortie le span en frequence doit etre assez important
# Deux caracteristiques a modifier :
#     Marge_frequence qui augmente ou diminue le span frequenciel donc l inverse en temporel
#     pas : nb de points
# Delta_t x Delta_f = N
#==============================================================================

Marge_frequence=5e12
Nb_Pts=3e6
fi=(193.45e12)-Marge_frequence
ff=(193.65e12)+Marge_frequence
pas=(ff-fi)/Nb_Pts
#pas=0.00000008e12
freq=np.arange(fi,ff,pas)
f_pulse=0.03e12
f_pulse2=f_pulse/2

E_0=1
E_02=np.sqrt(2)
#Triple impulsion gaussienne
def E(E_0,freq,f_pulse):
    return  E_0*np.exp(-2*np.log(2)*(((freq-(193.52e12))/f_pulse)**2))+(E_02*np.exp(-2*np.log(2)*(((freq-(193.58e12))/f_pulse2)**2)))
E=E(E_0,freq,f_pulse)


Lambda=scipy.constants.c/(freq)
Intensite_Lambda=E

Lambda=Lambda.tolist()
Lambda.reverse()
Lambda=np.array(Lambda)
Intensite_Lambda=Intensite_Lambda.tolist()
Intensite_Lambda.reverse()
Intensite_Lambda=np.array(Intensite_Lambda)

plt.figure(1)
plt.figure(1).suptitle('Pulse',fontsize=16,color='b')
plt.subplot(211)
plt.ylabel('Intensite')
plt.xlabel('frequence THz')
plt.axis([fi*1e-12, ff*1e-12, 0, 2.2])
plt.plot(freq*1e-12,abs(E)**2,color='b')
plt.subplot(212)
plt.ylabel('Intensite')
plt.xlabel('Longueur d onde (nm)')
#plt.axis([fi*1e-12, ff*1e-12, 0, 2.2])
plt.plot(Lambda*1e9,abs(Intensite_Lambda)**2,color='b')


Lambda=Lambda*1e9
Data=  np.column_stack((Lambda, abs(Intensite_Lambda)**2))
np.savetxt('Goda_freq.dat', Data,delimiter="; \t",fmt='%1.9e') 
plt.show()
