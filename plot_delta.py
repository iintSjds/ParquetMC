import os
import sys
import re
import glob
import numpy as np
from matplotlib import pyplot as plt

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

with open("inlist", "r") as file:
    line = file.readline()
    para = line.split(" ")
    BetaStr = para[1]
    Beta = float(BetaStr)
    rsStr = para[2]
    rs = float(rsStr)
    Mass2Str = para[3]
    Mass2 = float(Mass2Str)
    LambdaStr = para[4]
    Lambda = float(LambdaStr)
    ChargeStr = para[5]
    Charge2 = float(ChargeStr)
    TotalStep = float(para[7])

print rs, Beta, Mass2, Lambda, TotalStep
Miu=1.0

# 0: I, 1: T, 2: U, 3: S
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

kF = (9.0*np.pi/4.0)**(1.0/3.0)/rs
Nf = kF/2.0/np.pi**2
Bubble = 0.0971916  # 3D, Beta=10, rs=1
Step = None

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

TauBin=np.linspace(Beta/20.0,Beta*19.0/20.0,num=10)
TauBinSize=len(TauBin)



#for order in Order:
   # for chan in Channel:
files = os.listdir(folder)
Num = 0
Norm = 0
Data0 = None
DataList = []
#FileName = "vertex{0}_{1}_pid[0-9]+.dat".format(order, chan)
FileName = "delta0.dat"
for f in files:
    if re.match(FileName, f):
        print "Loading ", f
        Norm0 = -1
        d = None
        #try: 
        with open(folder+f, "r") as file:
        #file.open(folder+f, "r")
            print 1
            line0 = file.readline()
            Step = int(line0.split(":")[-1])/1000000
            # print "Step:", Step
            line1 = file.readline()
            # print line1
            Norm0 = float(line1.split(":")[-1])
            # print "Norm: ", Norm0
            line2 = file.readline()
            if FreqBin is None:
                FreqBin = np.fromstring(
                    line2.split(":")[1], sep=' ')
                FreqBinSize = len(FreqBin)
                print FreqBinSize
            line3 = file.readline()
            if AngleBin is None:
                AngleBin = np.fromstring(
                    line3.split(":")[1], sep=' ')
                AngleBinSize = len(AngleBin)
                print AngleBinSize
            
            line4 = file.readline()
            if ExtMomBin is None:
                ExtMomBin = np.fromstring(
                    line4.split(":")[1], sep=' ')
                ExtMomBinSize = len(ExtMomBin)
                ExtMomBin /= kF   
            d = np.loadtxt(folder+f)
            #if d is not None and Norm0 > 0:
                #if Data0 is None:
                    #Data0 = d
                #else:
                    # Data0 = d
                 #   Data0 += 
           
           
           ##Fourier Transformation of delta*GG 
            dd = d.reshape(
                  (FreqBinSize,AngleBinSize,ExtMomBinSize,2))/Norm0
           
            GG=np.zeros(FreqBinSize*ExtMomBinSize)
            gg=GG.reshape(
                    (FreqBinSize,ExtMomBinSize))

            for i0 in range (FreqBinSize):
                for j0 in range (ExtMomBinSize):
                    E=ExtMomBin[j0]*ExtMomBin[j0]-1.0
                    gg[i0][j0]=1.0 #/(FreqBin[i0]*FreqBin[i0]+E*E)
            
            
            print gg.shape, dd.shape  
            for i0 in range (AngleBinSize):
                for j0 in range (2):
                    dd[:,i0,:,j0]=np.multiply(gg,dd[:,i0,:,j0])
            print dd.shape
            
            Fourier=np.exp(1j*2*np.pi*np.tensordot(FreqBin,TauBin,0))
            print Fourier.shape
            print Fourier
           
            ff=np.tensordot(Fourier,dd,axes=([0,0]))
            print ff.shape
            ax1=plt.axes()
            plt.plot(TauBin,np.real(ff[:,0,0,0]),'.');
            plt.savefig('new_3.pdf')

            plt.show()

             #  Norm += Norm0
            #   dd = d.reshape(
            #       (AngleBinSize, ExtMomBinSize, 2))/Norm0
            #   dd = AngleIntegation(dd, 0)
            #   DataList.append(SpinMapping(dd))
                       
       # except:
        #    print "fail to load ", folder+f
        
#if Norm > 0 and Data0 is not None:
#        print "Total Weight: ", Data0[0]
#        Data0 /= Norm
#        Data0 = Data0.reshape((AngleBinSize, ExtMomBinSize, 2)
