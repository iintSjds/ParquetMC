import os
import sys
import re
import glob
import numpy as np
from matplotlib import pyplot as plt

g=0.7*2.0*4.0
Omega=0.5

rs = None
Lambda = None
Mass2 = None
Beta = None
Charge2 = None
TotalStep = None
BetaStr = None
rsStr = None
ChargeStr = None
Mass2Str = None
LambdaStr = None

folder = "./"

FreqBin=None
FreqBinSize=None
AngleBin = None
ExtMomBin = None
AngleBinSize = None
ExtMomBinSize = None
TauBin=None
TauBinSize=None
Data = {}  # key: (order, channel)
DataWithAngle = {}  # key: (order, channel)
DataErr = {}  # key: (order, channel)

kF = 1.0

with open("config.dat","r") as file:
    line0 = file.readline()
    print line0.split(" ")[-1]

    #read Beta from delta.dat
    Beta = float(line0.split(" ")[1])  
    print ("Beta", Beta)
    omega_c=float(line0.split(" ")[-1])  

    line2 = file.readline()
    if FreqBin is None:
        FreqBin = np.fromstring(
            line2, sep=' ')
        FreqBinSize = len(FreqBin)
        print (FreqBinSize)
    #line3 = file.readline()
    #if AngleBin is None:
    #    AngleBin = np.fromstring(
    #       line3.split(":")[1], sep=' ')
    #    AngleBinSize = len(AngleBin)
    #    print (AngleBinSize)
    line4 = file.readline()
    if ExtMomBin is None:
        ExtMomBin = np.fromstring(
            line4, sep=' ')
        ExtMomBinSize = len(ExtMomBin)
        ExtMomBin /= kF
print ExtMomBin

delta=np.zeros(FreqBinSize*ExtMomBinSize)
delta=delta.reshape((FreqBinSize,ExtMomBinSize))
TauBin=np.zeros(FreqBinSize)
TauBinSize=len(TauBin)
for i in range (TauBinSize):
    TauBin[i]=i*Beta/TauBinSize


# read f
TauBin=[]
F=[]
with open("f.dat","r") as file:
    TauBin=np.fromstring(file.readline(),sep=' ')
    F=np.fromstring(file.readline(), sep=' ')
    F=F.reshape((TauBinSize,ExtMomBinSize))

delta0=-g/4.0/np.pi*sum(ExtMomBin**2*F[0])/ExtMomBinSize*ExtMomBin[-1]

print ExtMomBin.shape, F.shape
# taudep is tau dependent part; taudep=cosh(...)*F(tau,k)*k**2 with k integrated
taudep= np.cosh(Omega*(0.5*Beta-np.abs(TauBin))) * np.tensordot(ExtMomBin**2,F,axes=([0,1]))/ExtMomBinSize*ExtMomBin[-1]
delta1=np.tensordot(taudep,np.cos(np.tensordot(FreqBin,TauBin,axes=0)),axes=([0,1]))/TauBinSize*Beta
delta1=delta1*g*Omega/8.0/np.pi/np.sinh(0.5*Beta*Omega)

with open("delta0.dat","w") as file:
    for i in range(FreqBinSize):
        for k in range(ExtMomBinSize):
            file.write("{0}\t".format(delta0+delta1[i]))



