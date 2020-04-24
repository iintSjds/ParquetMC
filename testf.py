import os
import sys
import re
import glob
import math
import numpy as np
from pynufft import NUFFT_cpu, NUFFT_hsa
from matplotlib import pyplot as plt

Channel = [0, 1, 2, 3]
# Channel = [3]
ChanName = {0: "I", 1: "T", 2: "U", 3: "S"}
# 0: total, 1: order 1, ...
Order = [0, ]

#folder = "./Beta{0}_rs{1}_lambda{2}/".format(BetaStr, rsStr, Mass2Str)
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
#Nf = kF/2.0/np.pi**2
#Bubble = 0.0971916  # 3D, Beta=10, rs=1
#Step = None

def AngleIntegation(Data, l):
    # l: angular momentum
    shape = Data.shape[1:]
    Result = np.zeros(shape)
    for x in range(AngleBinSize):
        # Result += Data[x, ...] * \
        #     np.cos(l*AngleBin[x])*2.0*np.pi/AngleBinSize
        Result += Data[x, ...]*2.0/AngleBinSize
    return Result/2.0
    # return Result

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

#for order in Order:
   # for chan in Channel:
files = os.listdir(folder)
Num = 0
Norm = 0
Data0 = None
DataList = []
#FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)
FileName = "delta0.dat"

if(np.all(FreqBin[:-1] <= FreqBin[1:])):
    print ("fine")
else:
    print ("error")


#initialize delta
delta=np.zeros(FreqBinSize*ExtMomBinSize)
delta=delta.reshape((FreqBinSize,ExtMomBinSize))
cut=np.searchsorted(FreqBin,omega_c,side='right') 
delta[0:cut,:]=1.0
delta[cut:,:]=-0.1
modulus=math.sqrt(np.tensordot(delta[0:cut,:],delta[0:cut,:],axes=([0,1],[0,1])))
delta[0:cut,:]=delta[0:cut,:]/modulus
TauBin=np.zeros(FreqBinSize)
TauBinSize=len(TauBin)
for i in range (TauBinSize):
    TauBin[i]=i*Beta/TauBinSize

IterationType=1
loopcounter=0
lamu=0
shift=2.0
with open("f.dat","w") as file:
    for i in range(TauBinSize):
        file.write("{0} ".format(TauBin[i]))
        file.write("\n")
    for i in range(TauBinSize):
        for k in range(ExtMomBinSize):
            file.write("{0}\t".format( np.exp(-ExtMomBin[k]**2)*np.sin(2*np.pi*TauBin[i]) ))
