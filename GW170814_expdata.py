# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 10:23:28 2020

@author: laugo
"""
#----------------------
# Import needed modules
#----------------------
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.interpolation import shift
import h5py
filename = "L-L1_LOSC_CLN_4_V1-1186741845-32.hdf5"
filename = "L-L1_GWOSC_4KHZ_R1-1186741846-32.hdf5"
#filename= "H-H1_LOSC_4_V2-1135136334-32.hdf5"
#filename = "H-H1_LOSC_CLN_4_V1-1186741845-32.hdf5"
#filename = "L-L1_LOSC_C00_4_V1-1186741845-32.hdf5"


# Open the files
data = h5py.File(filename, "r")
odata0 = pd.read_csv("GW170814_owndata0.csv").values
odat05 = pd.read_csv("GW170814_owndata05.csv").values
odat15 = pd.read_csv("GW170814_owndata15.csv").values
datasc = pd.read_csv("ScatterDataset.csv", sep=';', decimal=',').values
datasc= datasc[datasc[:,0].argsort()]
print(datasc)
plottype = "png"
eventname = "GW170814"


plt.figure("GW170814 strain")
#plt.plot(data[:,0],data[:,1], label="L1_data", color="g")
plt.plot(datasc[:,0]-0.52, datasc[:,1], '-o',color="b",label="L1 data")
plt.plot(odat15[:,0],shift(odat15[:,1], +10, cval=np.NaN), label="own template 1.5PN", color="y")
plt.plot(odat05[:,0],shift(odat05[:,1], +10, cval=np.NaN), label="own template 0.5PN", color="r")
plt.plot(odata0[:,0],shift(odata0[:,1], +5, cval=np.NaN), label="own template 0PN", color="g")
plt.xlim([-0.08,0.01])
plt.xlabel('Time (s) since event')
plt.ylabel('Strain')
#plt.title("L1 whitened data around the event GW170814")
plt.legend(loc='upper left')
plt.savefig(eventname+'_comparison_earlyfit.'+plottype)
plt.show()


# Explore the file
for key in data.keys():
    print(key)

#---------------------
# Read in strain data
#---------------------
strain = data['strain']['Strain'].value
ts = data['strain']['Strain'].attrs['Xspacing']
#-----------------------
# Print out some meta data
#-----------------------
print("\n\n")
metaKeys = data['meta'].keys()
meta = data['meta']
for key in metaKeys:
    print(key, meta[key].value)
#---------------------------
# Create a time vector
#---------------------------
gpsStart = meta['GPSstart'].value
duration = meta['Duration'].value #1186741861.53
gpsEnd   = gpsStart + duration
t_c=1186741861.53
print(gpsStart,gpsStart-1186741861.53,gpsEnd-1186741861.53)
time = np.arange(gpsStart-t_c, gpsEnd-t_c, ts)
#----------------------
# Plot the time series
#----------------------
numSamples = len(time)
plt.figure("GW170814 strain data")
plt.plot(time[0:numSamples], strain[0:numSamples]*10**(17),label="L1_data", color="g")
#plt.plot(odat15[:,0],shift(odat15[:,1]/30, +10, cval=np.NaN), label="own template 1.5PN", color="y")
plt.xlabel('Time (s) since event')
plt.ylabel('Strain(*10^(17))')
#plt.title("L1 whitened data around the event GW170814")
plt.legend(loc='upper left')
plt.show()


#print("\n len data:\n", data[:,3].size, odat15[:,1].size)
plottype = "png"
eventname = "GW170814"


#plt.figure("GW170814 strain")
#plt.plot(data[:,0],data[:,1], label="L1_data", color="g")
##plt.plot(data[:,0],data[:,3], label="H1_template", color="b")
##plt.plot(odat05[:,0],odat05[:,1], label="own template 0.5PN", color="r")
##plt.xlim([-0.1,0.001])
#plt.plot(odat15[:,0],shift(odat15[:,1], 12, cval=np.NaN), label="own template 1.5PN", color="y")
##plt.plot(time[t1:t2],shift(odat15[:,1], 15, cval=np.NaN), label="own template 1.5PN", color="y")
##plt.xlim([-0.2,0.001])
##plt.plot(odata0[:,0],odata0[:,1], label="own template 0PN", color="g")
##plt.xlim([-0.1,0.001])
##plt.plot(odat15[:,0],odatay15, label="own template 1.5", color="g")
#plt.xlabel("Time (s) since event")
#plt.ylabel("strain (*10^(19))")
##plt.title("H1 whitened data around the event GW151226")
#plt.legend(loc='upper left')
#plt.savefig(eventname+'_comparison_test.'+plottype)
#plt.show()