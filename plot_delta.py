import os
import sys
import re
import glob
import math
import numpy as np
from pynufft import NUFFT_cpu, NUFFT_hsa
from matplotlib import pyplot as plt



#rs = None
#Lambda = None
#Mass2 = None
#Beta = None
#Charge2 = None
#TotalStep = None
#BetaStr = None
#rsStr = None
#ChargeStr = None
#Mass2Str = None
#LambdaStr = None
#
#with open("inlist", "r") as file:
#    line = file.readline()
#    para = line.split(" ")
#    BetaStr = para[1]
#    Beta = float(BetaStr)
#    rsStr = para[2]
#    rs = float(rsStr)
#    Mass2Str = para[3]
#    Mass2 = float(Mass2Str)
#    LambdaStr = para[4]
#    Lambda = float(LambdaStr)
#    ChargeStr = para[5]
#    Charge2 = float(ChargeStr)
#    TotalStep = float(para[7])
#
#print rs, Beta, Mass2, Lambda, TotalStep
#Miu=1.0
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
for loopcounter in range(100):

    gg=np.zeros((FreqBinSize,ExtMomBinSize))

    for i0 in range (FreqBinSize):
        for j0 in range (ExtMomBinSize):
            E=ExtMomBin[j0]*ExtMomBin[j0]-1.0
            gg[i0][j0]=1.0 /(FreqBin[i0]*FreqBin[i0]+E*E)
            #gg[i0][j0]=1.0 /FreqBin[i0] / (ExtMomBin[j0]**4+1)
    
    
    #print (gg.shape) 
    f_freq=np.multiply(gg,delta)
    #print (delta.shape)

    phase_shift=np.exp(-1j*np.pi/Beta*TauBin)
    ff=phase_shift[:,np.newaxis]*np.fft.fft(f_freq,axis=0)/Beta
    #ff=np.fft.fft(delta,axis=0)

   # print (ff.shape)
    
    plt.figure()
    ax1=plt.axes()
    plt.plot(TauBin,ff.real[:,0],'.')
    plt.savefig('new_3.pdf')
    plt.close()
    #plt.show()
    
    with open("f.dat","w") as file:
        for i in range(TauBinSize):
            file.write("{0} ".format(TauBin[i]))
        file.write("\n")
        for i in range(TauBinSize):
            for k in range(ExtMomBinSize):
                file.write("{0}\t".format(2*ff.real[i][k]))
    
    #if(loopcounter%5==0):
    #     IterationType=(IterationType+1)%2
    print IterationType
    myCmd='./feyncalc.exe >out.txt'
    #myCmd='python phonon.py>out.txt'
    os.system(myCmd)

    d = None
    for f in files:
        if re.match(FileName, f):
            #print ("Loading ")
            Norm0 = -1
            d = None
            #try: 
            with open(folder+f, "r") as file:
                #line1 = file.readline()
                # print (line1)
                #Norm0 = float(line1.split(":")[-1])
                # print ("Norm: ", Norm0)
                #file.open(folder+f, "r")
                d = np.loadtxt(folder+f)
                plt.figure()
                #plt.axvline(omega_c)
                #plt.plot(FreqBin,delta[:,0],'ko')
                plt.plot(ExtMomBin,delta[0,:],'r.')
                plt.savefig('delta_new.pdf')
                plt.close()
                #if d is not None and Norm0 > 0:
                    #if Data0 is None:
                        #Data0 = d
                    #else:
                        # Data0 = d
                    #   Data0 += 
            
   
#       print (TauBin)
    d=d.reshape((FreqBinSize,ExtMomBinSize))
    #test=np.sin((Beta*FreqBin-np.pi)/FreqBinSize)
# print (test)
# dd[:,0,0,0]=dd[:,0,0,0]+test
    #separate delta in to high and low frequency 

    
    #print (FreqBin)
    lamu=0.0
    if(IterationType==0):
        delta[cut:,:]=d[cut:,:]
    elif(IterationType==1): 
        d[0:cut,:]=d[0:cut,:]+shift*delta[0:cut,:]
        lamu=np.tensordot(d[0:cut,:],delta[0:cut,:],axes=([0,1],[0,1]))
        print lamu-shift
        modulus=math.sqrt(np.tensordot(d[0:cut,:],d[0:cut,:],axes=([0,1],[0,1])))
        print ("modulus:",modulus)
        delta[0:cut,:]=d[0:cut,:]/modulus  
        



    plt.figure()
    #plt.axvline(omega_c)
    #plt.plot(FreqBin,delta[:,0],'ko')
    plt.plot(ExtMomBin,delta[0,:],'r.')
    plt.savefig('delta_q.pdf')
    plt.close()
    #plt.show()

    plt.figure()
    #plt.axvline(omega_c)
    plt.plot(FreqBin[0:cut],delta[0:cut,0],'ko')
    plt.plot(FreqBin[cut:],delta[cut:,0]*lamu,'ko')
    #plt.plot(ExtMomBin,delta[0,:],'r.')
    plt.savefig('delta_w.pdf')
    plt.close()


    ##Fourier Transformation of delta*GG 

# om=np.random.randn(1512,1) #None-Uniform-Fourier-Transform
# Nd=(256,)
# Kd=(512,)
# Jd=(6,)
    
# time_data = np.zeros(256, )
# time_data[96:128+32] = 1.0

# NufftObj=NUFFT_cpu()
# NufftObj.plan(om,Nd,Kd,Jd)
# ff=NufftObj.forward(time_data)
        
# Fourier=np.exp(1j*2*np.pi*np.tensordot(FreqBin,TauBin,0))
# print Fourier.shape
# print Fourier
#  ff=np.tensordot(Fourier,dd,axes=([0,0]))
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

   
