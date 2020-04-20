import os
import sys
import numpy as np

Order=1
Beta=2000.0
Temp=1.0/Beta
rs=1.91916
Mass2=1.0
Lambda=1.0
charge2=1.0
MaxExtMom=7.0
TotalStep=10
seed=1453
Duplication=0
Omega_c=0.5

Freqsize=2000
ExtMomBinSize=80
Freq=np.pi*Temp
Mom=1e-7
step=(MaxExtMom-Mom)/(ExtMomBinSize-1)
with open("config.dat","w") as file:
    file.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(Order, Beta, rs, Mass2, Lambda, charge2, MaxExtMom, TotalStep, seed, Duplication,Omega_c))    
    for k in range(Freqsize):
        file.write("{0}\t".format(Freq))
        Freq=Freq+2*np.pi*Temp
    file.write("\n")
    for k in range(ExtMomBinSize):
        file.write("{0}\t".format(Mom))
        Mom=Mom+step
    file.write("\n")

    


        
