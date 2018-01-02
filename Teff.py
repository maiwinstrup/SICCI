import math
import Tgr
import Pd
import numpy as np
import matplotlib.pyplot as plt


def Teff(Ta,hs,hi,freq):
    #compute the effective temperature i.e the physical temperature at penetration depth
    #the temperature profile
    Gs,Gi=Tgr.Tgr(Ta,hs,hi)
    #C-shaped salinity profile? hold on, wait a bit
    
    #set up profile
    layer_number_snow=int(math.ceil(hs*100.0))
    if hs==0.0: layer_number_snow=0
    layer_number_ice=int(math.ceil(hi*100.0))
    if hi==0.0: layer_number_ice=0
    layer_number_all=layer_number_snow + layer_number_ice
    num=range(layer_number_all+2)
    
    si=np.zeros(layer_number_all+2)
    salinity=np.zeros(layer_number_all+2)
    density=np.zeros(layer_number_all+2)
    temperature=np.zeros(layer_number_all+2)
    scatter=0.07*np.ones(layer_number_all+2)
    scatter[0]=0.0
    rms=np.zeros(layer_number_all+2)
    d=0.01*np.ones(layer_number_all+2)
    wc=np.zeros(layer_number_all+2)

    Tb=np.zeros(90)
    emissivity=np.zeros(90)
    Tb_=np.zeros(90)
    emis_=np.zeros(90)
    Tv=np.zeros(90)
    Th=np.zeros(90)
    
    
    si[0:layer_number_snow]=0
    si[layer_number_snow:layer_number_all]=1
    salinity[0:layer_number_snow]=0.0
    density[1:layer_number_snow]=300.0
    density[layer_number_snow:layer_number_all]=920.0
    for s in range(layer_number_snow): temperature[s]=Ta+Gs*0.01*s
    for i in range(layer_number_snow,layer_number_all): temperature[i]=Ta+Gs*0.01*layer_number_snow+Gi*0.01*(i-layer_number_snow+1)
    for i in range(layer_number_snow,layer_number_all): salinity[i]=16.0/(i-layer_number_snow+1)**2 + 6.0
    
    frequency=np.array([1.4,6,10,18,23,36,50,89,118,150,183,243,325,448,664])
    for i in range(0,90):
    #ryd op i returparametre
         Tb[i],emissivity[i],Tb_[i],emis_[i],Tv[i],Th[i]=Pd.emod(num,temperature,rms,density,d,scatter,salinity,wc,si,i,freq)
         print i,Th,Tv, Tb,emissivity,Tb_,emis_
    plt.plot(Tv,'x',Th,'o',Tb,'x',Tb_,'o')
    plt.show()
    return Tv,Th
