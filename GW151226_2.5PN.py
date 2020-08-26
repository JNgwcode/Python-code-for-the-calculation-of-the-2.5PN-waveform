# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:18:30 2020

@author: laugo
"""

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
    return 1/4*Theta(m1,m2,tau)**(-1/4)*(1+(743/4032+11/48*nu(m1,m2))*Theta(m1,m2,tau)**(-1/4)-1/5*pi*Theta(m1,m2,tau)**(-3/8)+(19583/254016+24401/193536*nu(m1,m2)+31/288*(nu(m1,m2))**2)*(Theta(m1,m2,tau))**(-1/2)+(-11891/53760+109/1920*nu(m1,m2))*pi*(Theta(m1,m2,tau))**(-5/8))
    #checked
def gamma(m1,m2,tau,r):
    #return G*(m1+m2)/(r*c**2) lowest order
    return x(m1,m2,tau)*(1+(1-nu(m1,m2)/3)*x(m1,m2,tau)+(1-65/12*nu(m1,m2))*(x(m1,m2,tau))**2)
    #checked
def phi(m1,m2,tau):
    return -1/nu(m1,m2)*(Theta(m1,m2,tau))**(5/8)*(1+(3715/8064+55/96*nu(m1,m2))*Theta(m1,m2,tau)**(-1/4)-3/4*pi*(Theta(m1,m2,tau))**(-3/8)+(9275495/14450688+284875/258048*nu(m1,m2)+1855/2048*(nu(m1,m2))**2)*(Theta(m1,m2,tau))**(-1/2)+(-38645/172032+65/2048*nu(m1,m2))*pi*(Theta(m1,m2,tau))**(-5/8)*np.log(Theta(m1,m2,tau)/(nu(m1,m2)*c**3/(5*G*(m1+m2))*(t_c))))
    #checked
def omega_s(m1,m2,tau,r):
    #return (G*(m1+m2)/r**3*(1+(-3+nu(m1,m2))*gamma(m1,m2,tau,r))+(6+41/4*nu(m1,m2)+(nu(m1,m2))**2)*(gamma(m1,m2,tau,r))**2)**(1/2)
    return c**3/(8*G*(m1+m2))*((Theta(m1,m2,tau))**(-3/8)+(743/2688+11/32*nu(m1,m2))*(Theta(m1,m2,tau))**(-5/8)-3*pi/10*(Theta(m1,m2,tau))**(-3/4)*+(1855099/14450688+59675/258048*nu(m1,m2)+371/2048*(nu(m1,m2))**2)*(Theta(m1,m2,tau))**(-7/8))#2PN Blanchet 1996
    #checked
def psi(m1,m2,r,tau,omega0):
    #return phi(m1,m2,tau)-2*G*omega_s(m1,m2,tau,r)/c**3*np.log(omega_s(m1,m2,tau,r)/omega0) #Maggiore
    return phi(m1,m2,tau)-2*(x(m1,m2,tau))**(3/2)*(1-nu(m1,m2)/2*x(m1,m2,tau))*np.log(omega_s(m1,m2,tau,r)/omega0) #Blanchet=Maggiore+correction: nu/2*x
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
def Hplus2(m1,m2,r,tau,iota,omega0):
    return pi*np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(psi(m1,m2,r,tau,omega0))*(-5/8-1/8*(np.cos(iota))**2)+np.cos(2*psi(m1,m2,r,tau,omega0))*(11/60+33/10*(np.cos(iota))**2+29/24*(np.cos(iota))**4-1/24*(np.cos(iota))**6+nu(m1,m2)*(353/36-3*(np.cos(iota))**2-251/72*(np.cos(iota))**4+5/24*(np.cos(iota))**6)+(nu(m1,m2))**2*(-49/12+9/2*(np.cos(iota))**2-7/24*(np.cos(iota))**4-5/24*(np.cos(iota))**6))+pi*np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(3*psi(m1,m2,r,tau,omega0))*(27/8*(1+(np.cos(iota))**2))+np.cos(4*psi(m1,m2,r,tau,omega0))*(118/15+16/5*(np.cos(iota))**2+86/15*(np.cos(iota))**4+16/15*(np.cos(iota))**6+nu(m1,m2)*(-262/9+16*(np.cos(iota))**2-166/9*(np.cos(iota))**4+16/3*(np.cos(iota))**6)+(nu(m1,m2))**2*(14+16*(np.cos(iota))**2-10/3*(np.cos(iota))**4-16/3*(np.cos(iota))**6))+np.cos(6*psi(m1,m2,r,tau,omega0))*(-81/40*(np.sin(iota))**4*(1+(np.cos(iota))**2)*(1-5*nu(m1,m2)+5*(nu(m1,m2))**2))+np.sin(iota)*(m1-m2)/(m1+m2)*np.sin(psi(m1,m2,r,tau,omega0))*(11/40+5*np.log(2)/4+(np.cos(iota))**2*(7/40+np.log(2)/4))+np.sin(iota)*(m1-m2)/(m1+m2)*np.sin(3*psi(m1,m2,r,tau,omega0))*((-189/40+27/4*np.log(3/2))*(1+(np.cos(iota))**2))

#pi*np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(psi(m1,m2,r,tau,omega0))*(-5/8-1/8*(np.cos(iota))**2)
#+np.cos(2*psi(m1,m2,r,tau,omega0))*(11/60+33/10*(np.cos(iota))**2+29/24*(np.cos(iota))**4-1/24*(np.cos(iota))**6+nu(m1,m2)*(353/36-3*(np.cos(iota))**2-251/72*(np.cos(iota))**4+5/24*(np.cos(iota))**6)+(nu(m1,m2))**2*(-49/12+9/2*(np.cos(iota))**2-7/24*(np.cos(iota))**4-5/24*(np.cos(iota))**6))
#+pi*np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(3*psi(m1,m2,r,tau,omega0))*(27/8*(1+(np.cos(iota))**2))
#+np.cos(4*psi(m1,m2,r,tau,omega0))*(118/15+16/5*(np.cos(iota))**2+86/15*(np.cos(iota))**4+16/15*(np.cos(iota))**6+nu(m1,m2)*(-262/9+16*(np.cos(iota))**2-166/9*(np.cos(iota))**4+16/3*(np.cos(iota))**6)+(nu(m1,m2))**2*(14+16*(np.cos(iota))**2-10/3*(np.cos(iota))**4-16/3*(np.cos(iota))**6))
#+np.cos(6*psi(m1,m2,r,tau,omega0))*(-81/40*(np.sin(iota))**4*(1+(np.cos(iota))**2)*(1-5*nu(m1,m2)+5*(nu(m1,m2))**2))
#+np.sin(iota)*(m1-m2)/(m1+m2)*np.sin(psi(m1,m2,r,tau,omega0))*(11/40+5*np.log(2)/4+(np.cos(iota))**2*(7/40+np.log(2)/4))
#+np.sin(iota)*(m1-m2)/(m1+m2)*np.sin(3*psi(m1,m2,r,tau,omega0))*((-189/40+27/4*np.log(3/2))*(1+(np.cos(iota))**2))

def Hcross2(m1,m2,r,tau,iota,omega0):
    return np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.cos(psi(m1,m2,r,tau,omega0))*(-9/20-3/2*np.log(2))+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.cos(3*psi(m1,m2,r,tau,omega0))*(189/20-27/2*np.log(3/2))-np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*(3*pi/4)*np.sin(psi(m1,m2,r,tau,omega0))+np.cos(iota)*np.sin(2*psi(m1,m2,r,tau,omega0))*(17/15+113/30*(np.cos(iota))**2-1/4*(np.cos(iota))**4+nu(m1,m2)*(143/9-245/18*(np.cos(iota))**2+5/4*(np.cos(iota))**4)+(nu(m1,m2))**2*(-14/3-35/6*(np.cos(iota))**2-5/4*(np.cos(iota))**4))+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(3*psi(m1,m2,r,tau,omega0))*(27*pi/4)+np.cos(iota)*np.sin(4*psi(m1,m2,r,tau,omega0))*(44/3+268/15*(np.cos(iota))**2+16/5*(np.cos(iota))**4+nu(m1,m2)*(-476/9+620/9*(np.cos(iota))**2-16*(np.cos(iota))**4)+(nu(m1,m2))**2*(68/3-116/3*(np.cos(iota))**2+16*(np.cos(iota))**4))+np.cos(iota)*np.sin(6*psi(m1,m2,r,tau,omega0))*(-81/20*(np.sin(iota))**4*(1-5*nu(m1,m2)+5*(nu(m1,m2))**2))

#+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.cos(3*psi(m1,m2,r,tau,omega0))*(189/20-27/2*np.log(3/2))
#-np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*(3*pi/4)*np.sin(psi(m1,m2,r,tau,omega0))
#+np.cos(iota)*np.sin(2*psi(m1,m2,r,tau,omega0))*(17/15+113/30*(np.cos(iota))**2-1/4*(np.cos(iota))**4+nu(m1,m2)*(143/9-245/18*(np.cos(iota))**2+5/4*(np.cos(iota))**4)+(nu(m1,m2))**2*(-14/3-35/6*(np.cos(iota))**2-5/4*(np.cos(iota))**4))
#+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(3*psi(m1,m2,r,tau,omega0))*(27*pi/4)
#+np.cos(iota)*np.sin(4*psi(m1,m2,r,tau,omega0))*(44/3+268/15*(np.cos(iota))**2+16/5*(np.cos(iota))**4+nu(m1,m2)*(-476/9+620/9*(np.cos(iota))**2-16*(np.cos(iota))**4)+(nu(m1,m2))**2*(68/3-116/3*(np.cos(iota))**2+16*(np.cos(iota))**4))
#+np.cos(iota)*np.sin(6*psi(m1,m2,r,tau,omega0))*(-81/20*(np.sin(iota))**4*(1-5*nu(m1,m2)+5*(nu(m1,m2))**2))

def Hplus52(m1,m2,r,tau,iota,omega0):
    return np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(psi(m1,m2,r,tau,omega0))*(1771/5120+1667/5120*(np.cos(iota))**2+217/9216*(np.cos(iota))**4-1/9216*(np.cos(iota))**6+nu(m1,m2)*(681/256+13/768*(np.cos(iota))**2-35/768*(np.cos(iota))**4+1/2304*(np.cos(iota))**6)+(nu(m1,m2))**2*(-3451/9216+673/3072*(np.cos(iota))**2-5/9216*(np.cos(iota))**4-1/3072*(np.cos(iota))**6))+pi*np.cos(2*psi(m1,m2,r,tau,omega0))*(19/3+3*(np.cos(iota))**2+2/3*(np.cos(iota))**4+nu(m1,m2)*(-16/3+14/3*(np.cos(iota))**2+2*(np.cos(iota))**4))+np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(3*psi(m1,m2,r,tau,omega0))*(3537/1024+22977/5120*(np.cos(iota))**2-15309/5120*(np.cos(iota))**4+729/5120*(np.cos(iota))**6+nu(m1,m2)*(-23829/1280+5529/1280*(np.cos(iota))**2+7749/1280*(np.cos(iota))**4-729/1280*(np.cos(iota))**6)+(nu(m1,m2))**2*(29127/5120-27267/5120*(np.cos(iota))**2-1647/5120*(np.cos(iota))**4+2187/5120*(np.cos(iota))**6))+np.cos(4*psi(m1,m2,r,tau,omega0))*(-16*pi/3*(np.sin(iota))**2*(1+(np.cos(iota))**2)*(1-3*nu(m1,m2)))+np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(5*psi(m1,m2,r,tau,omega0))*(-108125/9216+40625/9216*(np.cos(iota))**2+83125/9216*(np.cos(iota))**4-15625/9216*(np.cos(iota))**6+nu(m1,m2)*(8125/256-40625/2304*(np.cos(iota))**2-48125/2304*(np.cos(iota))**4-15625/2304*(np.cos(iota))**6)+(nu(m1,m2))**2*(-119375/9216+40625/3702*(np.cos(iota))**2+44375/9216*(np.cos(iota))**4+15625/3072*(np.cos(iota))**6))+(m1-m2)/(m1+m2)*np.cos(7*psi(m1,m2,r,tau,omega0))*(-117649/46080*(np.sin(iota))**5*(1+(np.cos(iota))**2)*(1-4*nu(m1,m2)+3*(nu(m1,m2))**2))+np.sin(2*psi(m1,m2,r,tau,omega0))*(-9/5+14/5*(np.cos(iota))**2+7/5*(np.cos(iota))**4+nu(m1,m2)*(96/5-8/5*(np.cos(iota))**2-28/5*(np.cos(iota))**4))+np.sin(iota)**2*(1+(np.cos(iota))**2)*np.sin(4*psi(m1,m2,r,tau,omega0))*(56/5-32*np.log(2)/3-nu(m1,m2)*(1193/30-32*np.log(2)))

#+pi*np.cos(2*psi(m1,m2,r,tau,omega0))*(19/3+3*(np.cos(iota))**2+2/3*(np.cos(iota))**4+nu(m1,m2)*(-16/3+14/3*(np.cos(iota))**2+2*(np.cos(iota))**4))
#+np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(3*psi(m1,m2,r,tau,omega0))*(3537/1024+22977/5120*(np.cos(iota))**2-15309/5120*(np.cos(iota))**4+729/5120*(np.cos(iota))**6+nu(m1,m2)*(-23829/1280+5529/1280*(np.cos(iota))**2+7749/1280*(np.cos(iota))**4-729/1280*(np.cos(iota))**6)+(nu(m1,m2))**2*(29127/5120-27267/5120*(np.cos(iota))**2-1647/5120*(np.cos(iota))**4+2187/5120*(np.cos(iota))**6))
#+np.cos(4*psi(m1,m2,r,tau,omega0))*(-16*pi/3*(np.sin(iota))**2*(1+(np.cos(iota))**2)*(1-3*nu(m1,m2)))
#+np.sin(iota)*(m1-m2)/(m1+m2)*np.cos(5*psi(m1,m2,r,tau,omega0))*(-108125/9216+40625/9216*(np.cos(iota))**2+83125/9216*(np.cos(iota))**4-15625/9216*(np.cos(iota))**6+nu(m1,m2)*(8125/256-40625/2304*(np.cos(iota))**2-48125/2304*(np.cos(iota))**4-15625/2304*(np.cos(iota))**6)+(nu(m1,m2))**2*(-119375/9216+40625/3702*(np.cos(iota))**2+44375/9216*(np.cos(iota))**4+15625/3072*(np.cos(iota))**6))
#+(m1-m2)/(m1+m2)*np.cos(7*psi(m1,m2,r,tau,omega0))*(-117649/46080*(np.sin(iota))**5*(1+(np.cos(iota))**2)*(1-4*nu(m1,m2)+3*(nu(m1,m2))**2))
#+np.sin(2*psi(m1,m2,r,tau,omega0))*(-9/5+14/5*(np.cos(iota))**2+7/5*(np.cos(iota))**4+nu(m1,m2)*(96/5-8/5*(np.cos(iota))**2-28/5*(np.cos(iota))**4))
#+np.sin(iota)**2*(1+(np.cos(iota))**2)*np.sin(4*psi(m1,m2,r,tau,omega0))*(56/5-32*np.log(2)/3-nu(m1,m2)*(1193/30-32*np.log(2)))

def Hcross52(m1,m2,r,tau,iota,omega0):
    return 6/5*(np.sin(iota))**2*np.cos(iota)*nu(m1,m2)+np.cos(iota)*np.cos(2*psi(m1,m2,r,tau,omega0))*(2-22/5*(np.cos(iota))**2+nu(m1,m2)*(-154/5+94/5*(np.cos(iota))**2))+(np.sin(iota))**2*np.cos(iota)*np.cos(2*psi(m1,m2,r,tau,omega0))*(-112/5-64/3*np.log(2)+nu(m1,m2)*(1193/15+64*np.log(2)))+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(psi(m1,m2,r,tau,omega0))*(-913/7680+1891/11520*(np.cos(iota))**2-7/4608*(np.cos(iota))**4+nu(m1,m2)*(1165/384-235/576*(np.cos(iota))**2+7/1152*(np.cos(iota))**4)+(nu(m1,m2))**2*(-1301/4608+301/2304*(np.cos(iota))**2-7/1536*(np.cos(iota))**4))+pi*np.cos(iota)*np.sin(2*psi(m1,m2,r,tau,omega0))*(34/3-8/3*(np.cos(iota))**2-nu(m1,m2)*(20/3-8*(np.cos(iota))**2))+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(3*psi(m1,m2,r,tau,omega0))*(12501/2560-12069/1280*(np.cos(iota))**2+1701/2560*(np.cos(iota))**4+nu(m1,m2)*(-19581/640-7821/320*(np.cos(iota))**2+1701/640*(np.cos(iota))**4)+(nu(m1,m2))**2*(18903/2560+11403/1280*(np.cos(iota))**2+5103/2560*(np.cos(iota))**4))+(np.sin(iota))**2*np.cos(iota)*np.sin(4*psi(m1,m2,r,tau,omega0))*(-32*pi/3*(1-3*nu(m1,m2)))+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(5*psi(m1,m2,r,tau,omega0))*(-101875/4608+6875/256*(np.cos(iota))**2+21875/4608*(np.cos(iota))**4+nu(m1,m2)*(66875/1152-44375/576*(np.cos(iota))**2+21875/1152*(np.cos(iota))**4)+(nu(m1,m2))**2*(-100625/4608+83125/2304*(np.cos(iota))**2-21875/1536*(np.cos(iota))**4))+(np.sin(iota))**5*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(7*psi(m1,m2,r,tau,omega0))*(117649/23040*(1-4*nu(m1,m2)+3*(nu(m1,m2))**2))

#+np.cos(iota)*np.cos(2*psi(m1,m2,r,tau,omega0))*(2-22/5*(np.cos(iota))**2+nu(m1,m2)*(-154/5+94/5*(np.cos(iota))**2))
#+(np.sin(iota))**2*np.cos(iota)*np.cos(2*psi(m1,m2,r,tau,omega0))*(-112/5-64/3*np.log(2)+nu(m1,m2)*(1193/15+64*np.log(2)))
#+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(psi(m1,m2,r,tau,omega0))*(-913/7680+1891/11520*(np.cos(iota))**2-7/4608*(np.cos(iota))**4+nu(m1,m2)*(1165/384-235/576*(np.cos(iota))**2+7/1152*(np.cos(iota))**4)+(nu(m1,m2))**2*(-1301/4608+301/2304*(np.cos(iota))**2-7/1536*(np.cos(iota))**4))
#+pi*np.cos(iota)*np.sin(2*psi(m1,m2,r,tau,omega0))*(34/3-8/3*(np.cos(iota))**2-nu(m1,m2)*(20/3-8*(np.cos(iota))**2))
#+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(3*psi(m1,m2,r,tau,omega0))*(12501/2560-12069/1280*(np.cos(iota))**2+1701/2560*(np.cos(iota))**4+nu(m1,m2)*(-19581/640-7821/320*(np.cos(iota))**2+1701/640*(np.cos(iota))**4)+(nu(m1,m2))**2*(18903/2560+11403/1280*(np.cos(iota))**2+5103/2560*(np.cos(iota))**4))
#+(np.sin(iota))**2*np.cos(iota)*np.sin(4*psi(m1,m2,r,tau,omega0))*(-32*pi/3*(1-3*nu(m1,m2)))
#+np.sin(iota)*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(5*psi(m1,m2,r,tau,omega0))*(-101875/4608+6875/256*(np.cos(iota))**2+21875/4608*(np.cos(iota))**4+nu(m1,m2)*(66875/1152-44375/576*(np.cos(iota))**2+21875/1152*(np.cos(iota))**4)+(nu(m1,m2))**2*(-100625/4608+83125/2304*(np.cos(iota))**2-21875/1536*(np.cos(iota))**4))
#+(np.sin(iota))**5*np.cos(iota)*(m1-m2)/(m1+m2)*np.sin(7*psi(m1,m2,r,tau,omega0))*(117649/23040*(1-4*nu(m1,m2)+3*(nu(m1,m2))**2))

def hplus(m1,m2,r,t,iota,omega0):
    return 2*G*redmass(m1,m2)*x(m1,m2,t)/(c**2*r)*(Hplus0(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**1/2*Hplus12(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))*Hplus1(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**3/2*Hplus32(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**2*Hplus2(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**5/2*Hplus52(m1,m2,r,t,iota,omega0))
    #checked
def hcross(m1,m2,r,t,iota,omega0):
    return 2*G*redmass(m1,m2)*x(m1,m2,t)/(c**2*r)*(Hcross0(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**1/2*Hcross12(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))*Hcross1(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**3/2*Hcross32(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**2*Hplus2(m1,m2,r,t,iota,omega0)+(x(m1,m2,t))**5/2*Hplus52(m1,m2,r,t,iota,omega0))
    #checked
t_c=1.1#1.5#0.8201#0.8#0.2#
NP=4507#3277+1#3360+1#3277+1#820+1#6144+1#
t=np.linspace(0,t_c-0.005,NP)

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
fncsv = 'GW151226_owndata25.csv'
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

print('integration constant Theta0: \n',nu(m1,m2)*c**3/(5*G*(m1+m2))*(t_c))