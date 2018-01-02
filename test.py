import numpy as np
import Pd

def test(q,freq):
    Tphys=260.0
    Ka=0.5
    Ks=0.5
    Ke=Ka+Ks
    theta=0.707 
    d=1.0
    k=(2.0*3.14159)/(0.3/freq)
    e=1.4+0.01j
    Tsnow=Tphys*(Ka/(Ke-q*Ks))*(1.0-np.exp((-1.0*Ke+q*Ks)*d/np.cos(theta)))
    ka=Pd.HWabsorp(e,k)
    print ka
    return Tsnow, k
