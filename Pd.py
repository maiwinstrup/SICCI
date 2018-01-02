import cmath
import math
import numpy as np

def transang(ia,e1,e2):
    #conputes the transmission angle using Snellius
    ta=math.asin((math.sqrt(e1.real)/math.sqrt(e2.real))*math.sin(ia))
    return ta

def Vb(T,sal):
    # the volume of brine in sea ice Ulaby et al. 1986, E71
    # T: thermometric temperature [C]
    # sal: sea ice salinity [ppt]
    if (abs(T) < 0.1):
         volbrine = 0.001*sal*49.717
    else:
         volbrine = 0.001*sal*((-49.185/T)+0.532)
    return volbrine

def extloss(d,ang,ke):
    #computes the loss coefficient in a layer with extinction, ke, and thickness, d.
    loss=(ke*d/math.cos(ang))
    return abs(loss)

def absloss(d,ang,ka):
    #computes the absorption loss coefficient in a layer with the absorption coefficient ka, and thickness, d.
    absorp=(ka*d/math.cos(ang))
    return abs(absorp)

def eice_rn2p(eh,ei,vi):
    # Polder / Van Santen dielectric constant of sea ice with random brine needles
    # e.g. Shokr (1998) Arctic sea ice in the microwave c-band, IEEE TGRS 36(2), 463-78.
    measureA=1.2
    measureB=1.2
    # initial estimate: Vant et al. (1978), J. Appl. Phys. 49(3), 1264-80.
    emi=complex((3.05+0.72*vi),(0.024+0.33*vi))
    m=0
    # iterate to solve for emi: the dielectric const. of the mixture
    while measureA > 0.001 and measureB > 0.0001:
        f1=(ei-eh)/(ei+emi)
        f2=5.0*emi+ei
        est=eh+((vi/3.0)*f1*f2)
        measureA=abs(emi.real-est.real)
        measureB=abs(emi.imag-est.imag)
        emi=est
        m=m+1
        if (m>30): break
    return emi

def epsib0(T,N):
    # static dielectric constant of brine, se. Ulaby et al. 1986, E65a
    # T: thermometric temperature [C]
    # N: normality of brine
    eps=88.045-0.4147*T+6.295E-4*(T**2)+1.075E-5*(T**3)
    a=1.0-0.255*N+5.15E-2*(N**2)-6.89E-3*(N**3)
    se=eps*a
    return se

def relaxt(T,N):
    # relaxation time of brine, Ulaby et al. 1986, E67
    # T: thermometric temperature
    # N: normality of brine
    rel=(1.1109E-10)-(3.824E-12)*T+(6.938E-14)*(T**2)-(5.096E-16)*(T**3)
    b=1.0-(0.146E-2)*T*N-(4.89E-2)*N-(2.97E-2)*(N**2)+(5.64E-3)*(N**3)
    relax=rel*b
    return relax

def Sb(T):
    # Salinity of brine
    # Vant et al. 1978 J. Appl. Phys. 49(3),1264-1280 & Ulaby et al. 1986 E63
    # T: thermometric temperature of brine [C]
    if T<-1.8 and T>=-8.2:
        s=1.725-18.756*T-0.3964*(T**2)
    elif T<-8.2 and T>=-22.9:
        s=57.041-9.929*T-0.16204*(T**2)-0.002396*(T**3)
    elif T<-22.9 and T>=-36.8:
        s=242.94+1.5299*T+0.0429*(T**2)
    elif T<-36.8 and T>=-43.2:
        s=508.18+14.535*T+0.2018*(T**2)
    else:
        s=34.2
    return s

def Nsw(Ssw):
    # normality of sea water or brine Ulaby et al. 1986, E20
    # Ssw: salinity of brine or sea water
    N = 0.91141*Ssw*(1.707E-2+1.205E-5*Ssw+4.058E-9*(Ssw**2))
    return N

def condbrine(T,N):
    # conductivity of brine Ulaby et al. 1986, vol. III E68,E69
    # T: thermometric temperature [C]
    # N: normality of brine
    D=25.0-T
    sig=N*(10.39-2.378*N+0.683*(N**2)-0.135*(N**3)+(1.01E-2)*(N**4))
    c=1.0-(1.96E-2)*D+(8.08E-5)*(D**2)-N*D*(3.02E-5)+(3.92E-5)*D+N*((1.72E-5)-(6.58E-6)*D)
    sigma=c*sig
    return sigma

def ebrine2(freq,T,Si):
    # permittivity and loss of brine Ulaby et al. 1986 E64a
    # freq: em frequency
    volb=Vb(T,Si)
    salb=Sb(T)
    N=Nsw(salb)
    conb=condbrine(T,N)
    relax=relaxt(T,N)
    eb0=epsib0(T,N)

    f=freq*1.0E9
    epsiwoo=4.9
    e0=8.854E-12

    eb=epsiwoo+((eb0-epsiwoo)/(1.0+((relax*f)**2)))
    ebi=((relax*f*(eb0-epsiwoo))/(1.0+((relax*f)**2)))+(conb/(6.28*f*e0))
    ebr=complex(eb,ebi)
    return ebr

def eice_s2p(e1,e2,v):
    #improved born approximation by C. Matzler (1998). J. Appl. Phys. 83(11),6111-7
    #Polder/VanSanten mixing formulae for spheical inclusions
    #effective dielectric constant of medium consisting of e1 and e2
    #e1: dielectric constant of background
    #e2: dielectric constant of sherical inclusion
    #v: fraction of inclusions

    eeff=0.25*(2.0*e1-e2+3.0*v*(e2-e1)+cmath.sqrt((2.0*e1-e2+3.0*v*(e2-e1))**2.0+8.0*e1*e2))
    return eeff

def epicesal(T,freq,sal):
    # Dielectric constant of pure ice: Matzler, 1998, Microwave properties of ice and snow, in:
    # B. Smitht et al. (eds.) solar system ices, 241-257, Kluwer.
    # T: thermometric temperature in K
    # freq: Frequency in GHz
    epi=3.1884+9.1E-4*(T-273.15)

    # The Hufford model for the imagenary part.
    if T == 0: T=273.15
    theta=300.0/T
    alpha=(0.00504+0.0062*theta)*math.exp(-22.1*theta)
    beta=(0.502-0.131*theta/(1+theta))*1e-4+(0.542e-6*((1+theta)/(theta+0.0073))**2)

    epii=(alpha/freq) + (beta*freq)

    epui=complex(epi,epii)

    eh=complex(1.0,0.0)
    v=Vb((T-273.15),sal)
    eb=ebrine2(freq,(T-273.15),sal)
    e=eice_rn2p(epui,eb,v)
    if sal > 0.0: e=epui
    return epui

def e_wet_snow(k,Dens,wc):
    # beregner dielektrisitetskonstanten af vad sne
    fws=k/20.94
    denws=Dens/1000.0
    wc=wc*100.0
    A1ws=1.0
    A2ws=1.0
    B1ws=0.0
    Aws=1.0+1.83*denws+0.02*A1ws*(wc**1.015)+B1ws
    Bws=0.073*A1ws
    Cws=0.073*A2ws
    f0=9.07 #relaxation frequency

    ews1=Aws+((Bws*wc)/(1.0+((fws/f0)**2)))
    ews2=(Cws*(fws/f0)*wc)/(1.0+((fws/f0)**2))
    e_ws=complex(ews1,ews2)
    return e_ws

def FT_(ia,ta,e1,e2):
    #computes the Transmission coefficient (HH), using angles and dielectrics
    nn1=cmath.sqrt(e1.real)
    nn2=cmath.sqrt(e2.real)
    T_=(2.0*nn2*cmath.cos(ia))/(nn2*cmath.cos(ta)+nn1*cmath.cos(ta))
    return abs(T_)

def FTll(ia,ta,e1,e2):
    #computes the transmission coefficient (VV) using the angle and dielectrics
    nn1=cmath.sqrt((e1.real))
    nn2=cmath.sqrt((e2.real))
    Tll=(2.0*nn2*cmath.cos(ia))/(nn2*cmath.cos(ta)+nn1*cmath.cos(ia))
    return abs(Tll)

def FRll(ia,ta,e1,e2):
    #computes the Fresnell reflection coefficient for vertical polarisation
    nn1=math.sqrt(e1.real)
    nn2=math.sqrt(e2.real)
    r=abs((nn1*math.cos(ia)-nn2*math.cos(ta))/(nn1*math.cos(ia)+nn2*math.cos(ta)))**2
    return r

def FR_(ia,ta,e1,e2):
    #computes the Fresnell reflection coefficient for horisontal polarisation
    nn1=math.sqrt(e1.real)
    nn2=math.sqrt(e2.real)
    r=abs((nn2*math.cos(ia)-nn1*math.cos(ta))/(nn2*math.cos(ia)+nn1*math.cos(ta)))**2
    return r

def HWabsorp(e,k):
    kk=cmath.sqrt(e)
    #print e, k*abs(cmath.sqrt(kk)), abs(kk.imag)
    #a=2.0*k*abs(kk)*abs(kk.imag)
    a=2.0*k*abs(kk.imag)
    return a

def iborn_s2p_sc(e1,e2,eeff,v,k,pcc):
    # improved born approximation by C. Matzler (1998). J. Appl. Phys. 83(11),6111-7
    #scattering coefficient of a collection of spherical
    #particles with correlation length pc using improved born approximation
    ss=(3.0*(pcc**3.0)*(k**4)/32.0)*v*(1.0-v)*abs(((e2-e1)*(2.0*eeff+e1))/(2.0*eeff+e2))**2
    return ss

def iborn_s2p_bc(e1,e2,eeff,v,k,pcc):
    # improved born approximation by C. Matzler (1998). J. Appl. Phys. 83(11),6111-7
    #backscattering coefficient of a collection of spherical
    #particles with correlation length pcc using improved born approximation
    R=3.0*pcc/4.0
    ss=((R**3.0)*(k**4.0/3.0))*v*(1.0-v)*abs(((e2-e1)*(2.0*eeff+e1))/(2.0*eeff+e2))**2.0
    return ss


def emod(num,TK,rms,Dens,d,pcc,sal,wc,si,ia,freq):
    d2rf=0.01745329
    k=(2.0*3.14159)/(0.3/freq)
    eair=1.0
    ang_ia=ia*d2rf
    n=len(num)

    #temperature in Kelvin or Celcius
    T=TK-273.15

    #correlation length in meters
    pcc=0.001*pcc

    #initialize arrays
    emipar=np.array(3)
    v=np.zeros(n)
    eh=np.array(np.zeros(n),dtype=complex)
    ep=np.array(np.zeros(n),dtype=complex)
    e=np.array(np.zeros(n),dtype=complex)
    ang=np.zeros(n)
    T_=np.zeros(n)
    Rll=np.zeros(n)
    R_=np.zeros(n)
    ks=np.zeros(n)
    bc=np.zeros(n)
    ka=np.zeros(n)
    ke=np.zeros(n)
    loss=np.ones(n)
    aloss=np.ones(n)

    #volume of scatters in snow
    for idx in range(n):
        if si[idx] == 0:
            v[idx]=Dens[idx]/926.0
        else:
            v[idx]=Vb(T[idx],sal[idx])

    #dielectric constant of background and particles
    for idx in range(n):
        if si[idx] == 0:
            eh[idx]=1.0
            ep[idx]=epicesal(TK[idx],freq,sal[idx])
        else:
            eh[idx]=epicesal(TK[idx],freq,0)
            ep[idx]=ebrine2(freq,T[idx],sal[idx])

    #effective dielectric constant using a two-phase PVS form. for spheres
    for idx in range(n):
        e[idx]=eice_s2p(eh[idx],ep[idx],v[idx])
        if (wc[idx] > 0.0001 and si[idx] == 0):
            e[idx]=e_wet_snow(k,Dens[idx],wc[idx]);

    #estimate the angle
    for idx in range(n):
        if idx==0:
            ang[idx]=abs(transang(ang_ia,eair,e[idx]))
        else:
            ang[idx]=abs(transang(ang[idx-1],e[idx-1],e[idx]))

    #Fresnell transmission coefficient
    for idx in range(n):
        if idx==0:
            T_[idx]=FTll(ang_ia,ang[idx],eair,e[idx])
            Rll[idx]=FRll(ang_ia,ang[idx],eair,e[idx])
            R_[idx]=FR_(ang_ia,ang[idx],eair,e[idx])
        else:
            T_[idx]=FTll(ang[idx-1],ang[idx],e[idx-1],e[idx])
            Rll[idx]=FRll(ang[idx-1],ang[idx],e[idx-1],e[idx])
            R_[idx]=FR_(ang[idx-1],ang[idx],e[idx-1],e[idx])

    #scattering and absorption coefficients and Loss
    for idx in range(n):
        ks[idx]=iborn_s2p_sc(eh[idx],ep[idx],e[idx],v[idx],k,pcc[idx])
        bc[idx]=iborn_s2p_bc(eh[idx],ep[idx],e[idx],v[idx],k,pcc[idx])
        ka[idx]=HWabsorp(e[idx],k)
        ke[idx]=ks[idx]+ka[idx]
        loss[idx]=extloss(d[idx],ang[idx],ke[idx])
        aloss[idx]=absloss(d[idx],ang[idx],ka[idx])

    bc=bc*d/np.cos(ang)

    Tbll,emissivityll = rtransf2(n,TK,Rll,d,8.68*ke,2.7)
    Tb_,emissivity_ = rtransf2(n,TK,R_,d,8.68*ke,2.7)


    nelements=len(d)
    ka=ka[1:nelements-2]
    ks=ks[1:nelements-2]
    d=d[1:nelements-2]
    ang=ang[1:nelements-2]
    Rll=Rll[1:nelements-2]
    R_=R_[1:nelements-2]
    TK=TK[1:nelements-2]
    #print "ka,ks,Rll,R_,ang,TK"
    #print ka,ks,Rll,R_,ang,TK

    Tbv = mulhut(8.68*ka,8.68*ks,d,ang,Rll,TK,271.35,2.7)
    Tbh = mulhut(8.68*ka,8.68*ks,d,ang,R_, TK,271.35,2.7)
    #Tbv = mulhut(ka,ks,d,ang,Rll,TK,271.35,2.7)
    #Tbh = mulhut(ka,ks,d,ang,R_, TK,271.35,2.7)
    return Tbll, emissivityll, Tb_, emissivity_, Tbv, Tbh

def mulhut(Ka,Ks,d,theta,r,Tphys,Tgrd,Tsky):
    #Lemmetyinen et al. 2010. Multiple-layer adaption of HUT snow emission model. IEEE TGRS 48(7), 2781-2794.
    #input: Ka,Ks,d,theta,Tgrd,Tsky,r,Tphys
    #output: Tb

    #print Ka[0], Ks[0], d[0], theta[0], r[0], Tphys[0]
    #empirical forward scattering factor, probably frequency dependent etc
    #unclear if it is related with the mie assymetry parameter
    #original HUT value is 0.96
    q=0.96
    #extinction
    Ke=Ks+Ka
    #loss
    loss=np.exp((Ke-q*Ks)*d/np.cos(theta))
    #bottom reflection, transmissivity
    rgrd=0.0
    tgrd=1.0-rgrd
    #transmissivity
    t=1.0-r
    #number of layers
    nlayers=len(Tphys)
    #define reflectivity array
    S=np.zeros(nlayers)
    #define A matrix
    A=np.identity(2*nlayers)
    #define b vector
    b=np.zeros(2*nlayers)
    #define effective temperature array
    Tsnow=np.zeros(nlayers)

    #the effective temperature of the layer
    #Tsnow[0]=0.0
    #Tsnow=Tphys*(Ka/(Ke-q*Ks))*(1.0-np.exp((-1.0*Ke+q*Ks)*d/np.cos(theta)))
    for n in range(0,nlayers-1):
        Tsnow[n]=Tphys[n]*((Ka[n])/(Ke[n]-q*Ks[n]))*(1.0-np.exp((-1.0*Ke[n]+q*Ks[n])*d[n]/np.cos(theta[n])))

    #top layer, is it really surface reflectivity squared? p.2783: ..from rn-1 and reduced to rn.
    S[0]=1.0/(1.0-(r[0]*r[1]/(loss[0]**2)))
    #intermediate layers
    for n in range(1,nlayers-1):
        S[n]=1.0/(1.0-(r[n]*r[n+1]/(loss[n]**2)))
    #bottom layer
    S[nlayers-1]=1.0/(1.0-(r[nlayers-1]*rgrd/(loss[nlayers-1]**2)))

    #set up the A matrix
    #top layer
    A[0,2]=S[0]*t[1]/loss[0]
    A[1,2]=S[0]*t[1]*r[0]/loss[0]**2
    #intermediate layers
    for i in range(1,nlayers-1):
        A[2*i,2*i-1]=S[i]*t[i]*r[i+1]/loss[i]**2
        A[2*i+1,2*i-1]=S[i]*t[i]/loss[i]
        A[2*i,2*i+2]=S[i]*t[i+1]/loss[i]
        A[2*i+1,2*i+2]=S[i]*t[i+1]*r[i]/loss[i]**2
    #bottom layer
    A[2*nlayers-2,2*nlayers-3]=S[nlayers-1]*t[nlayers-1]*rgrd/loss[nlayers-1]**2
    A[2*nlayers-1,2*nlayers-3]=S[nlayers-1]*t[nlayers-1]/loss[nlayers-1]

    #set up the b vector
    #top layer
    b[0]=S[0]*(Tsnow[0]*(1.0+r[1]/loss[0])+Tsky*t[0]*r[1]/loss[0]**2)
    b[1]=S[0]*(Tsnow[0]*(1.0+r[0]/loss[0])+Tsky*t[0]/     loss[0])
    #intermediate layers
    for i in range(1,nlayers-1):
        b[2*i]=  S[i]*(Tsnow[i]*(1.0+r[i+1]/loss[i]))
        b[2*i+1]=S[i]*(Tsnow[i]*(1.0+r[i]/  loss[i]))
    #bottom layers
    b[2*nlayers-2]=S[nlayers-1]*(Tsnow[nlayers-1]*(1.0+rgrd/loss[nlayers-1])+Tgrd*tgrd/loss[nlayers-1])
    b[2*nlayers-1]=S[nlayers-1]*(Tsnow[nlayers-1]*(1.0+r[nlayers-1]/loss[nlayers-1])+Tgrd*tgrd*r[nlayers-1]/loss[nlayers-1]**2)

    #copy A into A' (Ap) and assign b to left column in A' (Ap)
    Ap=np.copy(A)
    Ap[:,0]=b[:]

    #compute the brightness temperature
    #Tb=(1.0-r[0])*np.linalg.det(Ap)/np.linalg.det(A)
    (Apsign,Aplogdet)=np.linalg.slogdet(Ap)
    (Asign,Alogdet)=np.linalg.slogdet(A)
    taeller=abs(Apsign*np.exp(Aplogdet))
    naevner=Asign*np.exp(Alogdet)
    if naevner == 0: naevner=1.0
    Tb=(1.0-r[0])*taeller/naevner
    return Tb

def rtransf2(N,T,R,path,extinction,Tsky):
    #computes the individual contributions for each layer,N, using a first-order non-coherent radiative transfer model
    #after Burke et al 1979: comparison of 2.8-... JGR 84(C1), 287-294. NB: correction of the bending angle (theta_0) compared to the article.
    #N: number of layers
    #T: thermometric temperature [K]
    #R: the reflection coefficient
    #loss: loss coefficient
    #Tsky=2.7

    #define arrays
    tf_tmp=np.zeros(N+2)
    tf_tmp[N+1]=271.35
    rf_tmp=np.zeros(N+2)
    afa_tmp=np.ones(N+2)
    afa1_tmp=np.ones(N+2)
    afa2_tmp=np.ones(N+2)
    Teff_tmp=np.zeros(N+2)
    Teff_tmp[N+1]=271.35
    r_tmp=np.zeros(N+2)
    path_tmp=np.zeros(N+2)
    extinction_tmp=np.zeros(N+2)
    #the emissivity is 1 minus the reflectivity of the surface :: 
    #actually it is corrected to the reflectivity of the system in the end
    emissivity=1.0-R[0]
    #setting infinite lower and upper halfspace for computational reasons
    #the physical temperature of the lower half-space 
    T_tmp=np.zeros(N+2)
    T_tmp[N+1]=271.35
    #the reflection at the lower boundary
    R_tmp=np.zeros(N+2)
    R_tmp[N+1]=0.5
    #filling in the blanks in the temporary arrays
    for idx in range(1,N+1):
        T_tmp[idx]=T[idx-1]
        R_tmp[idx]=R[idx-1]
        path_tmp[idx]=path[idx-1]
        extinction_tmp[idx]=extinction[idx-1]
    #computing the attenuation factor from layers above
    for i in range(1,N+2):
        afa1_tmp[i-1] = np.prod((1.0 - R_tmp[0:i]))
    for i in range(1,N+2):
        afa2_tmp[i-1] = np.prod(np.exp(-1.0 * np.sum(extinction_tmp[1:i-1]*path_tmp[1:i-1])))
    #afa2_tmp[1] = 0.5
    afa_tmp=(afa1_tmp*afa2_tmp)
    afa_tmp[0] = 0.0
    afa_tmp[N+1] = afa_tmp[N]*((1.0 - R_tmp[N+1]))
    #computing the layer self emission in terms of effective temperature
    for j in range(0,N+1):
        tf_tmp[j] = T_tmp[j] * (1.0 - np.exp(-1.0*extinction_tmp[j]*path_tmp[j])) * (1.0 + R_tmp[j+1] *np.exp(-1.0*extinction_tmp[j]*path_tmp[j]))
    tf_tmp[0] = 0.0
    #computing the layer contributions to the self emission temperature from a point above the surface
    for idx in range(0,N+2):
        r_tmp[idx]    =  R_tmp[idx] * afa_tmp[idx]
        Teff_tmp[idx] = tf_tmp[idx] * afa_tmp[idx]
        print afa_tmp[idx]
    Tb=np.sum(r_tmp)*Tsky+np.sum(Teff_tmp)
    Tb0=np.sum(r_tmp)*0.0+np.sum(Teff_tmp)
    Tb100=np.sum(r_tmp)*100.0+np.sum(Teff_tmp)
    emissivity=1.0-(Tb100-Tb0)/100.0
    return np.array([Tb, emissivity])
