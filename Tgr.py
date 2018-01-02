def Tgr(Ta,hs,hi):
    #Nakawo & Sinha, 1981, Growth rate and salinity profile of FY. J. Glaciol. 
    Tm=271.3
    ks=0.31
    ki=2.1
    Gs=(Tm-Ta)/((ks/ki)*hi+hs)
    Gi=(Tm-Ta)/(hi+(ki/ks)*hs)
    if hs == 0: Gs=0.0
    if hi == 0: Gi=0.0
    return Gs,Gi
