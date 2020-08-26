# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 17:57:49 2020

@author: laugo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi
from scipy.constants import G
from scipy.constants import c

def redmass(m1,m2):
    return (m1*m2)/(m1+m2)

def nu(m1,m2):
    return redmass(m1,m2)/(m1+m2)

def Theta(m1,m2,tau):
    return nu(m1,m2)*c**3/(5*G*(m1+m2))*tau
    #checked
def x(m1,m2,tau):
    return 1/4*Theta(m1,m2,tau)**(-1/4)*(1+(743/4032+11/48*nu(m1,m2))*Theta(m1,m2,tau)**(-1/4)-1/5*pi*Theta(m1,m2,tau)**(-3/8))#+(19583/254016+24401/193536*nu(m1,m2)+31/288*(nu(m1,m2))**2)*(Theta(m1,m2,tau))**(-1/2))
    #checked
def gamma(m1,m2,tau,r):
    #return G*(m1+m2)/(r*c**2) lowest order
    return x(m1,m2,tau)*(1+(1-nu(m1,m2)/3)*x(m1,m2,tau))#+(1-65/12*nu(m1,m2))*(x(m1,m2,tau))**2)
    #checked
def phi(m1,m2,tau):
    return -1/nu(m1,m2)*(Theta(m1,m2,tau))**(5/8)*(1+(3715/8064+55/96*nu(m1,m2))*Theta(m1,m2,tau)**(-1/4)-3/4*pi*(Theta(m1,m2,tau))**(-3/8))#+(9275495/14450688+284875/258048*nu(m1,m2)+1855/2048*(nu(m1,m2))**2)*(Theta(m1,m2,tau))**(-1/2))
    #checked
def omega_s(m1,m2,tau,r):
    return (G*(m1+m2)/r**3*(1+(-3+nu(m1,m2))*gamma(m1,m2,tau,r)))**(1/2)#+(6+41/4*nu(m1,m2)+(nu(m1,m2))**2)*(gamma(m1,m2,tau,r))**2))**(1/2)
    #checked
def psi(m1,m2,r,tau,omega0):
    return phi(m1,m2,tau)-2*G*omega_s(m1,m2,tau,r)/c**3*np.log(omega_s(m1,m2,tau,r)/omega0) #Maggiore
    #return phi(m1,m2,tau)-2*(x(m1,m2,tau))**(3/2)*(1-nu(m1,m2)/2*x(m1,m2,tau))*np.log(omega_s(m1,m2,tau,r)/omega0) #Blanchet=Maggiore+correction: nu/2*x
    #checked
def Hplus0(m1,m2,r,tau,iota,omega0):
    return -(1+(np.cos(iota))**2)*np.cos(2*psi(m1,m2,r,tau,omega0))-1/96*(np.sin(iota))**2*(17+(np.cos(iota))**2)
    #checked
def Hcross0(m1,m2,r,tau,iota,omega0):
    return -2*np.cos(iota)*np.sin(2*psi(m1,m2,r,tau,omega0))
    #checked
def Hplus12(m1,m2,r,tau,iota,omega0):
    return -np.sin(iota)/8*(m1-m2)/(m1+m2)*((5+(np.cos(iota))**2)*np.cos(psi(m1,m2,r,tau,omega0))-9*(1+(np.cos(iota))**2)*np.cos(3*psi(m1,m2,r,tau,omega0)))
    #checked
def Hcross12(m1,m2,r,tau,iota,omega0):
    return -3/4*np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*(np.sin(psi(m1,m2,r,tau,omega0))-3*np.sin(3*psi(m1,m2,r,tau,omega0)))
    #checked
def Hplus1(m1,m2,r,tau,iota,omega0):
    return np.cos(2*psi(m1,m2,r,tau,omega0))*(19/6+3/2*(np.cos(iota))**2-1/3*(np.cos(iota))**4+nu(m1,m2)*(-19/6+11/6*(np.cos(iota))**2+(np.cos(iota))**4))-np.cos(4*psi(m1,m2,r,tau,omega0))*(4/3*(np.sin(iota))**2*(1+(np.cos(iota))**2)*(1-3*nu(m1,m2)))
    #checked
def Hcross1(m1,m2,r,tau,iota,omega0):
    return np.cos(iota)*np.sin(2*psi(m1,m2,r,tau,omega0))*(17/3-4/3*(np.cos(iota))**2+nu(m1,m2)*(-13/3+4*(np.cos(iota))**2))+np.cos(iota)*(np.sin(iota))**2*np.sin(4*psi(m1,m2,r,tau,omega0))*(-8/3*(1-3*nu(m1,m2)))
    #checked
def Hplus32(m1,m2,r,tau,iota,omega0):
    return np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(psi(m1,m2,r,tau,omega0))*(19/64+5/16*(np.cos(iota))**2-1/192*(np.cos(iota))**4+nu(m1,m2)*(-49/96+1/8*(np.cos(iota))**2+1/96*(np.cos(iota))**4))+np.cos(2*psi(m1,m2,r,tau,omega0))*(-2*pi*(1+(np.cos(iota))**2))+np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(3*psi(m1,m2,r,tau,omega0))*(657/128-45/16*(np.cos(iota))**2+81/128*(np.cos(iota))**4+nu(m1,m2)*(225/64-9/8*(np.cos(iota))**2-81/64*(np.cos(iota))**4))+np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(5*psi(m1,m2,r,tau,omega0))*(625/384*(np.sin(iota))**2*(1+(np.cos(iota))**2)*(1-2*nu(m1,m2)))
    #checked
def Hcross32(m1,m2,r,tau,iota,omega0):
    return np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(psi(m1,m2,r,tau,omega0))*(21/32+5/96*(np.cos(iota))**2+nu(m1,m2)*(-23/48+5/48*(np.cos(iota))**2))-4*pi*np.cos(iota)*np.sin(2*psi(m1,m2,r,tau,omega0))+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(3*psi(m1,m2,r,tau,omega0))*(-603/64+135/64*(np.cos(iota))**2+nu(m1,m2)*(171/32-135/32*(np.cos(iota))**2))+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(5*psi(m1,m2,r,tau,omega0))*(625/192*(1-2*nu(m1,m2))*(np.sin(iota))**2)
    #checked
def hplus(m1,m2,r,t,iota,omega0):
    return 2*G*redmass(m1,m2)*x(m1,m2,t)/(c**2*r)*(Hplus0(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**1/2*Hplus12(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))*Hplus1(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**3/2*Hplus32(m1,m2,r,t,iota,omega0))
    #checked
def hcross(m1,m2,r,t,iota,omega0):
    return 2*G*redmass(m1,m2)*x(m1,m2,t)/(c**2*r)*(Hcross0(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**1/2*Hcross12(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))*Hcross1(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**3/2*Hcross32(m1,m2,r,t,iota,omega0))
    #checked
t_c=0.8#0.2#
NP=3277+1#820+1#
t=np.linspace(0,t_c-0.001,NP)

f_min=35
f_max=450

m_sun=1.989*10**(30)
m1=14.2*m_sun# <--values from the tutorial 19.64*m_sun/6.71*m_sun 
m2=7.5*m_sun# values from the paper

i=0
r=440*3.0875*10**19
o0=2*pi*f_min

hgw=(hplus(m1,m2,r,t_c-t,i,o0)+hcross(m1,m2,r,t_c-t,i,o0))*10**(19)
data =[t-t_c, hgw/2]
datcsv = np.array(data).transpose()
# make a csv filename, header, and format
fncsv = 'GW151226_owndata15.csv'
#headcsv = eventname+' time-'+str(tevent)+     ' (s),H1_data_whitened,L1_data_whitened,H1_template_whitened,L1_template_whitened'
fmtcsv = ",".join(["%10.6f"]*2)
np.savetxt(fncsv, datcsv, fmt=fmtcsv)#, header=headcsv)



plt.figure('Amplitude')
plt.plot(t-t_c, hgw, label="GW amplitude 1.5", color="g")

plt.xlabel("Time (s)")
plt.ylabel("Strain (10^(-16))")
plt.title("Chirp of the GW151226")
plt.legend()
plt.show()