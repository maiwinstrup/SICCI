#Two layer incoherent emissivity model after Ulaby et al 1981.
#Volume 1, chapter 4., p. 243

#             air
#-------------------------------------
#:::::::::::::snow::::::::::::::::::::
#-------------------------------------
#|||||||||||||ice|||||||||||||||||||||

from scipy.interpolate import interp1d
import numpy as np
import math

def ice_emissivity(wave,myf):
    #isemissiviteten er vigtig men svaer at faa paa plads, lige nu er det myis der ligner mest, men det bliver kompenseret af en alt for hoej Teff.
    #kanal mean            std            mean e[260K]    e+std           e-std
    #6v    249.75824       5.6748425      0.96060863      0.98243495      0.93878231
    #6h    225.39771       7.7126605      0.86691426      0.89657834      0.83725018
    #10v   247.02053       7.1304733      0.95007895      0.97750385      0.92265406
    #10h   224.34160       8.8444326      0.86285231      0.89686935      0.82883526
    #18v   242.74486       11.378086      0.93363406      0.97739593      0.88987219
    #18h   223.34190       11.973182      0.85900729      0.90505799      0.81
    #37v   229.07616       19.711571      0.88106215      0.95687588      0.80524841
    #37h   214.43431       18.533517      0.82474735      0.89603011      0.75346459
    #89v   215.14091       20.977036      0.82746503      0.90814593      0.74678412
    #89h   205.14236       19.954230      0.78900907      0.86575611      0.71226203
    x=   [1,6, 10, 18, 37, 89, 100]
    #Eppler et al. 1992: table 4.1
    efyv=[0.935,0.900,0.924,0.941,0.960,0.955,0.926]
    efyh=[0.850,0.840,0.876,0.888,0.910,0.913,0.886]
    #Eppler et al. 1992: table 4.1
    emyv=[0.926,0.925,0.890,0.850,0.787,0.764,0.680]
    emyh=[0.865,0.820,0.817,0.780,0.635,0.706,0.650]
    fxfyv=interp1d(x,efyv,kind='linear')
    fxfyh=interp1d(x,efyh,kind='linear')
    fxmyv=interp1d(x,emyv,kind='linear')
    fxmyh=interp1d(x,emyh,kind='linear')
    ice_emis_v=myf*fxmyv(wave)+(1.0-myf)*fxfyv(wave)
    ice_emis_h=myf*fxmyh(wave)+(1.0-myf)*fxfyh(wave)
    return ice_emis_v,ice_emis_h

def e_snow_dry(k,dens_snow):
    fds=1.0E9*k/20.94
    mfcgs=dens_snow/1000.  
    e1ds=(1.+0.51*mfcgs)**3
    vi=dens_snow/916.
    e2i=(91.5-3.15)/(fds/7.23E3)
    e2ds=((3.*vi*e2i)*((e1ds)**2)*(2.*e1ds+1.0))/((3.15+2.0*e1ds)*(3.15+2.0*((e1ds)**2)))
    return e1ds,e2ds

def e_snow_wet(k,dens_snow,vol_water):
    fws=k/20.94
    denws=dens_snow/1000.
    A1ws=1.0
    A2ws=1.0
    B1ws=0.0
    Aws=1.0+1.83*denws+0.02*A1ws*(vol_water**1.015)+B1ws
    Bws=0.073*A1ws
    Cws=0.073*A2ws
    f0=9.07 #relaxation frequency
    ews1=Aws+((Bws*vol_water)/(1.0+((fws/f0)**2)))
    ews2=(Cws*(fws/f0)*vol_water)/(1.0+((fws/f0)**2))
    return ews1,ews2

def angle(e1,e2,t1):
    t2=math.asin((math.sqrt(e1)/math.sqrt(e2))*math.sin(t1))
    return t2

def rough_reflec_h(e1,e2,t1,rms,k):
    #check her!
    # rough surface reflection horisontal (Ulaby et al. 11.20-21)
    Gh=((math.cos(t1)-math.sqrt(e2-(math.sin(t1))**2))/(math.cos(t1)+math.sqrt(e2-(math.sin(t1))**2)))**2
    h=4.0*(k**2)*(rms**2)
    angexp=math.exp(-h*(math.cos(t1))**2)
    return Gh*angexp

def rough_reflec_v(e1,e2,t1,rms,k):
    #check her!
    # rough surface reflection vertical (Ulaby et al. 11.20-21)
    Gv=((e2*math.cos(t1)-math.sqrt(e2-(math.sin(t1))**2))/(e2*math.cos(t1)+math.sqrt(e2-(math.sin(t1))**2)))**2
    h=4.0*(k**2)*(rms**2)
    angexp=math.exp(-h*(math.cos(t1))**2)
    return Gv*angexp

def Kscat_ray(eb,ep,k,r,vs):
    # scattering coefficient using the Rayleigh scattering function (no dense media)
    N=(vs*3.0)/(4.0*3.14*(r**3))
    ec=ep/eb
    Kf=(ec-1)/(ec+2)
    return N*4.0*3.14*(r**6)*(k**4)*(Kf**2)

def Kscat(eb,ep,k,r,vs):
    # the snow scattering coefficeint using the improved born approximation
    eice=ep
    pci=r
    vfi=vs
    epseff = (2.0-eice+3.0*vfi*(eice-1.0) + math.sqrt((2.0-eice+3.0*vfi*(eice-1.0))**2.0+8.0*eice))/4.0
    return (3.0/32.0)*((0.001*pci)**3.0)*(k**4.0)*vfi*(1.0-vfi)*abs((2.0*epseff+1.0)*(eice-1.0)/(2.0*epseff+eice))**2.0

def Kabs(vs,k,e2_snow):
    # absorbtion coefficient
    wl=2*3.14/k
    alpha=(k)*vs*math.sqrt(e2_snow)
    return 2.0*alpha

def Loss(d,t1,Ke):
    # the scattering and absorbtion loss
    return math.exp(Ke*d/math.cos(t1))

def ice_property(myf):
    my_sal=0.5 #ppt
    fy_sal=8.0 #ppt
    my_dens=720 #kg/m3
    my_dens=926 #kg/m3
    salinity=my_sal*myf+(1-myf)*fy_sal
    density=my_dens*myf+(1-myf)*fy_dens
    return salinity,density

def emice_amsr(r,snow_depth,iceT,myf):
    #rms_snow=0.0028 #[rms/mm]
    #rms_ice=0.005 #[rms/mm]
    rms_snow=0.0 #[rms/mm]
    rms_ice=0.0 #[rms/mm]
    snowT=iceT*0.98 #[K]
    #r=0.10 #[mm]
    dens_snow=320 #[kg/m3]
    vol_water=0 #[%]
    ia=53.0
    #frequencies=np.array([6.93, 10.65, 18.70, 23.80, 36.50, 50.30, 52.80, 89.00])
    Tv6, Th6,  ev6, eh6 = emimod(ia, 6.93,rms_snow,rms_ice,r,snowT,dens_snow,snow_depth,vol_water,iceT,myf)
    Tv10,Th10, ev10, eh10 = emimod(ia,10.65,rms_snow,rms_ice,r,snowT,dens_snow,snow_depth,vol_water,iceT,myf)
    Tv18,Th18, ev18, eh18 = emimod(ia,18.70,rms_snow,rms_ice,r,snowT,dens_snow,snow_depth,vol_water,iceT,myf)
    Tv23,Th23, ev23, eh23 = emimod(ia,23.80,rms_snow,rms_ice,r,snowT,dens_snow,snow_depth,vol_water,iceT,myf)
    Tv36,Th36, ev36, eh36 = emimod(ia,36.50,rms_snow,rms_ice,r,snowT,dens_snow,snow_depth,vol_water,iceT,myf)
    Tv89,Th89, ev89, eh89 = emimod(ia,89.00,rms_snow,rms_ice,r,snowT,dens_snow,snow_depth,vol_water,iceT,myf)
    Tb=np.array([Tv6,Th6,Tv10,Th10,Tv18,Th18,Tv23,Th23,Tv36,Th36,Tv89,Th89])
    return Tb

def emimod(incidence_angle,frequency,snow_roughness,ice_roughness,snow_correlation_length,snow_physical_temp,snow_density,snow_depth,snow_water_content,ice_physical_temp,multiyear_ice_frac):
    #Snow and ice parameters
    #;;;;;;;;;;;;;;;;;;;;;;;
    t1=0.01744*incidence_angle #[deg]
    GHZ=frequency #[GHz]
    rms_snow=snow_roughness #[rms/mm]
    rms_ice=ice_roughness #[rms/mm]
    r=snow_correlation_length #[mm]
    T_snow=snow_physical_temp #[K]
    dens_snow=snow_density #[kg/m3]
    d=snow_depth #[m]
    vol_water=snow_water_content #[%]
    T_ice=ice_physical_temp #[K]
    myf=multiyear_ice_frac

    e_air=1.0 #dielectric c of air
    ep=3.15 #dielectric c of fresh ice
    #e_ice=3.5 #dielectric c of ice
    e_ice=2.5+(1.0-myf)
    k=20.94*GHZ
    vs=dens_snow/910.0 #volume of scatters

    #call emissivity function
    ev_fy,eh_fy=ice_emissivity(GHZ,myf)

    #Brightness temperature of the ice layer
    #can be extended with a physical model
    Tb_ice_v=T_ice*ev_fy
    Tb_ice_h=T_ice*eh_fy

    #The dielectric constant of snow
    [e1_snow, e2_snow]=e_snow_dry(k,dens_snow)  #for dry snow
    if vol_water > 0.0: [e1_snow, e2_snow]=e_snow_wet(k,dens_snow,100.0*vol_water)  #for wet snow
    #The refraction angle
    t2=math.asin((math.sqrt(e_air)/math.sqrt(e1_snow))*math.sin(t1)) 
    #The reflection at the surface
    G1h=rough_reflec_h(e_air,e1_snow,t1,rms_snow,k)
    G1v=rough_reflec_v(e_air,e1_snow,t1,rms_snow,k)
    G2h=rough_reflec_h(e1_snow,e_ice,t2,rms_ice,k)
    G2v=rough_reflec_v(e1_snow,e_ice,t2,rms_ice,k)
    #The scattering coefficient
    Ks=Kscat(e1_snow,ep,k,r,vs)
    #The absorbtion coefficient
    Ka=Kabs(vs,k,e2_snow)
    #The extinction coefficient
    Ke=Ks+Ka
    #The scattering albedo
    a=Ks/Ke
    #The loss
    L2=Loss(d,t2,Ke)

    #horisontally polarised brightness temperature
    F1h=(1.0-G1h)/(1.0-(G1h*G2h/L2**2))
    F2h=(1.0+(G2h/L2))*(1.0-(1.0/L2))*(1.0-a)*T_snow+(((1.0-G2h)*Tb_ice_h)/L2)
    Tbh=F1h*F2h

    #horisontally polarised emissivity
    Eh=(1.0+(G2h/L2))*(1.0-(1.0/L2))*(1.0-a)+((1.0-G2h)/L2)
    emissivityh=F1h*Eh

    #vertically polarised brightness temperature
    F1v=(1.0-G1v)/(1.0-(G1v*G2v/L2**2))
    F2v=(1.0+(G2v/L2))*(1.0-(1.0/L2))*(1.0-a)*T_snow+(((1.0-G2v)*Tb_ice_v)/L2)
    Tbv=F1v*F2v

    #vertically polarised emissivity
    Ev=(1.0+(G2v/L2))*(1.0-(1.0/L2))*(1.0-a)+((1.0-G2v)/L2)
    emissivityv=F1v*Ev

    return Tbv, Tbh, emissivityv, emissivityh

