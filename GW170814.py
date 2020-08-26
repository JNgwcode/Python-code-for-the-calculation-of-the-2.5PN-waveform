# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 17:30:32 2020

@author: laugo
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi
from scipy.constants import G
from scipy.constants import c

def Mchirp(m1,m2):
    return (m1*m2)**(3/5)/(m1+m2)**(1/5)  
    
def f(tau,M_c):
    return 1/pi*(5/(256*tau))**(3/8)*(G*M_c/c**3)**(-5/8)

def Ncyc(M_c,f_min,f_max):
    return 1/(32*pi**(8/3))*(G*M_c/c**3)**(-5/3)*(f_min**(-5/3)-f_max**(-5/3))

def Phi(tau,M_c,Phi0):
    return -2*(5*G*M_c/c**3)**(-5/8)*tau**(5/8)+Phi0

def hplus(tau,M_c,r,iota,Phi0):
    return 1/r*(G*M_c/c**2)**(5/4)*(5/(c*tau))**(1/4)*np.cos(Phi(tau,M_c,Phi0))*(1+(np.cos(iota))**2)/2

def hcross(tau,M_c,r,iota,Phi0):
    return 1/r*(G*M_c/c**2)**(5/4)*(5/(c*tau))**(1/4)*np.sin(Phi(tau,M_c,Phi0))*np.cos(iota)

t_c=1
NP=2048
t=np.linspace(0,t_c-0.001,NP)

f_min=35
f_max=450

m_sun=1.989*10**(30)
m1=30.5*m_sun
m2=25.3*m_sun
M_c=Mchirp(m1,m2)

fGW=f(t_c-t,M_c)

print(M_c)
print(fGW)
print(Ncyc(M_c,f_min,f_max),Ncyc(M_c,fGW[0],fGW[2046]))

plt.figure('Frquency signal')
plt.plot(t-t_c, fGW, label="GW frequency", color="g")

plt.xlabel("Time (s)")
plt.ylabel("f (Hz)")
plt.title("Chirp of the GW170814")
plt.legend()
plt.show()


i=0
Phi0=0#pi/2+1.43
r=540*3.0875*10**16
hgw=(hplus(t_c-t,M_c,r,i,Phi0)+hcross(t_c-t,M_c,r,i,Phi0))#*10**(16)
print(hgw)

plt.figure('Amplitude')
plt.plot(t-t_c, hgw, label="GW amplitude", color="r")

plt.xlabel("Time (s)")
plt.ylabel("Strain (10^(-16))")
plt.title("Chirp of the GW170814")
plt.legend()
plt.show()


