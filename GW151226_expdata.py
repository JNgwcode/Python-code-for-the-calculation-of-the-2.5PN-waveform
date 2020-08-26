# -*- coding: utf-8 -*-
"""

@author: Johann Ioannou-Nikolaides
"""
# Standard python library imports:

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.interpolation import shift

#-----------------------------------------------------------------------------


# Import data from csv files:

data = pd.read_csv("GW151226_data.csv").values
odata0 = pd.read_csv("GW151226_owndata0.csv").values
odat05 = pd.read_csv("GW151226_owndata05.csv").values
odat15 = pd.read_csv("GW151226_owndata15.csv").values
odat25 = pd.read_csv("GW151226_owndata25.csv").values
templ = pd.read_csv("GW151226_template.csv").values
time = pd.read_csv("GW151226_f_GW_time.csv").values

plottype = "png"
eventname = "GW151226"

#-----------------------------------------------------------------------------

# Calculate time stamp of the GW entering the detector math the computed waveform
# to the one taken from the 'BINARY BLACK HOLE SIGNALS IN LIGO OPEN DATA' tutorial 
# by the GW Open Science Center: 

t1=0
for i in range(len(time)):
    if time[i]>-1.1:
        print("\n The timestep -1.1 lies at: ",i)
        t1=i
        break
t2=0
for i in range(len(time)):
    if time[i]>0:
        print("\n The timestep 0 lies at: ",i)
        t2=i
        break
print('t2-t1: \n', t2-t1)

#-----------------------------------------------------------------------------



time=time.transpose()
time=time[0,:]
time=time.transpose()

plt.figure("GW 151226 strain")
#plt.style.use('seaborn-ticks')
plt.plot(time,templ, label="H1_template", color="b")
plt.xlim([-1,0.001])
#plt.plot(odat05[:,0],shift(odat05[:,1]*0.95,12, cval=np.NaN), label="own template 0.5PN", color="r")
#plt.xlim([-0.2,0.001])
#plt.plot(odat15[:,0],shift(odat15[:,1], 12, cval=np.NaN), label="own template 1.5PN", color="y")
#plt.plot(time[t1:t2],shift(odat15[:,1]*0.8, 21, cval=np.NaN), label="own template 1.5PN", color="y")#0.2=64716
plt.plot(time[t1:t2],shift(odat25[:,1]*(0.79),0, cval=np.NaN), label="own template 2.5PN", color = '#ffa500')
#plt.xlim([-0.2,0.001])
#plt.plot(odata0[:,0],shift(odata0[:,1]*(-0.95), 27, cval=np.NaN), label="own template 0PN", color="g")
plt.xlim([-1,0.001])
#plt.plot(odat15[:,0],odatay15, label="own template 1.5", color="g")
plt.xlabel("Time (s) since event")
plt.ylabel("strain (*10^(19))")
#plt.title("H1 whitened data around the event GW151226")
plt.legend(loc='upper left')
#plt.savefig(eventname+'_comparison_2515PNearlyfit.'+plottype)
plt.show()


#h1=hilbert(odat25[:,1]*(-0.79))
#h2=hilbert(templ[t1:t2])
#p1=np.angle(h1)
#p2=np.angle(h2)
#instantaneous_frequency = (np.diff(p1)/(2.0*np.pi)*odat25[:,1].shape[-1])
#print(p1,'\n',p2)
#p=np.unwrap(p2)-np.unwrap(p1)
#t3=t1+1
#plt.figure('Phase difference')
#plt.plot(time[t3:t2],instantaneous_frequency, label="H1_template", color="b")
##plt.plot(time[t1:t2],p, label="H1_template", color="b")
##plt.plot(time[t1:t2],p2, label="H1_template", color="b")
##plt.plot(time[t1:t2],p1, label="own_template", color="r")
#plt.xlabel("Time (s) since event")
#plt.ylabel("Instantaneous Phase difference")
#plt.show()

template=templ.transpose()
template=template[0,:]
data_extracted = template[t1:t2]
print("length of data extracted",len(data_extracted), "\nlength of odat: \n", len(odat15))
#
#ratio05 = np.divide(data_extracted,shift(odat05[:,1]*0.95,12, cval=np.NaN))
#ratio05=np.where(ratio05>5,5,ratio05)
#ratio05=np.where(ratio05<-5,-5,ratio05)
#meanratio05 = np.abs(np.mean(ratio05.reshape(-1, 5), axis=1)-1)
#
#ratio15 = np.divide(data_extracted,shift(odat15[:,1]*0.8, 21, cval=np.NaN))
#ratio15=np.where(ratio15>5,5,ratio15)
#ratio15=np.where(ratio15<-5,-5,ratio15)

ratio25 = np.divide(data_extracted,shift(odat25[:,1]*(0.79),0, cval=np.NaN))
ratio25=np.where(ratio25>3,3,ratio25)
ratio25=np.where(ratio25<-3,-3,ratio25)

#ratio15y = np.divide(data[:,1],odatay15)-1
#meanratio15 =np.abs(np.mean(ratio15.reshape(-1, 5), axis=1)-1)
#meanratio15y =np.abs(np.mean(ratio15y.reshape(-1, 8), axis=1))
meanratio25 =np.abs(np.mean(ratio25.reshape(-1, 6), axis=1)-1)

#meantime = np.mean(odat15[:,0].reshape(-1,5), axis=1)
meantime = np.mean(odat25[:,0].reshape(-1,6), axis=1)
#print("length of meantime :",len(meantime),len(meanratio15))

plt.figure("ratio plot")
#plt.plot(meantime,meanratio05, label="0.5PN waveform", color="b")
plt.plot(meantime,meanratio25, label="0.5PN waveform", color="b")
#plt.xlim([-0.04,0.001])
#plt.ylim([0,5])
#plt.plot(meantime,meanratio15, label="1.5PN waveform", color="y")
plt.xlim([-0.8,0.001])
#plt.ylim([0,5])
##plt.plot(meantime,meanratio15y, label="1.5pn waveform", color="g")
#
plt.xlabel("Time (s) since event")
plt.ylabel("ratio")
#plt.title("Ratio plot around the event GW151226")
plt.legend(loc='upper left')
#plt.savefig(eventname+'_ratio_plot25.'+plottype)
plt.show()

