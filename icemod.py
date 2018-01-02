import numpy as np
import reg_mod13 as rm

def Ts_amsr(Tb):
   Tb6v = np.float(Tb[0]) #6v
   Tb6h = np.float(Tb[1]) #6h
   Tb10v = np.float(Tb[2]) #10v
   Tb10h = np.float(Tb[3]) #10h
   Tb19v = np.float(Tb[4]) #19v
   Tb19h = np.float(Tb[5]) #19h
   Tb23v = np.float(Tb[6]) #23v
   Tb23h = np.float(Tb[7]) #23h
   Tb37v = np.float(Tb[8]) #37v
   Tb37h = np.float(Tb[9]) #37h
   Tb89v = np.float(Tb[10]) #89v
   Tb89h = np.float(Tb[11]) #89h

   #Lise's equations
   SDsim = 1.7701 + 0.017462*Tb6v - 0.02801*Tb19v + 0.0040926*Tb37v
   Tsi6 = 1.144*Tb6v - 0.815/SDsim - 27.08
   Tsi10 = 1.101*Tb10v - 0.999/SDsim - 14.28
   Tsi = 0.5*(Tsi6 + Tsi10)

   Teff6v  = 0.888  * Tsi + 30.245
   Teff10v = 0.901  * Tsi + 26.569
   Teff19v = 0.920  * Tsi + 21.536
   Teff23v = 0.932  * Tsi + 18.417
   Teff37v = 0.960  * Tsi + 10.902
   Teff50v = 0.989  * Tsi + 2.959
   Teff89v = 1.0604 * Tsi - 16.384
   Teff=np.array([Teff6v,Teff10v,Teff19v,Teff23v,Teff37v,Teff50v,Teff50v,Teff89v])
   return Teff

def ev_ice(Tb,IST,SD):
   Tb6v = np.float(Tb[0]) #6v
   Tb6h = np.float(Tb[1]) #6h
   Tb10v = np.float(Tb[2]) #10v
   Tb10h = np.float(Tb[3]) #10h
   Tb19v = np.float(Tb[4]) #19v
   Tb19h = np.float(Tb[5]) #19h
   Tb23v = np.float(Tb[6]) #23v
   Tb23h = np.float(Tb[7]) #23h
   Tb37v = np.float(Tb[8]) #37v
   Tb37h = np.float(Tb[9]) #37h
   Tb89v = np.float(Tb[10]) #89v
   Tb89h = np.float(Tb[11]) #89h

   GR0610v=(Tb10v-Tb6v)/(Tb10v+Tb6v)
   GR0618v=(Tb19v-Tb6v)/(Tb19v+Tb6v)
   GR0618h=(Tb19h-Tb6h)/(Tb19h+Tb6h)
   GR0636v=(Tb37v-Tb6v)/(Tb37v+Tb6v)

   #Lise's equations
   SDsim = 1.7701 + 0.017462*Tb6v - 0.02801*Tb19v + 0.0040926*Tb37v

   Tsi6 = 1.144*Tb6v - 0.815/SDsim - 27.08
   Tsi10 = 1.101*Tb10v - 0.999/SDsim - 14.28
   Tsi = 0.5*(Tsi6 + Tsi10)

   Teff6v  = 0.888  * Tsi + 30.245
   Teff10v = 0.901  * Tsi + 26.569
   Teff19v = 0.920  * Tsi + 21.536
   Teff23v = 0.932  * Tsi + 18.417
   Teff37v = 0.960  * Tsi + 10.902
   Teff50v = 0.989  * Tsi + 2.959
   Teff89v = 1.0604 * Tsi - 16.384

   #Rasmus equation for ice thickness
   ITsim= 4.81691061767 - 47.3872663622*GR0610v - 52.8250933515*GR0618v + 49.9619055273*GR0618h + 33.4020392591*GR0636v + 0.142250481385*Tb6h - 0.115797387369*Tb19h - 0.0385743856297*Tb37v
   ITsim=max(0.3,ITsim)
   ITsim=min(3.5,ITsim)

   #Rasmus forward model
   Tbsim=rm.reg_mod13(IST,SD,ITsim)

   e6v=Tbsim[0]/Teff6v
   e6v=max(0.0,e6v)
   e6v=min(1.0,e6v)

   e10v=Tbsim[2]/Teff10v
   e10v=max(0.0,e10v)
   e10v=min(1.0,e10v)

   e18v=Tbsim[4]/Teff19v
   e18v=max(0.0,e18v)
   e18v=min(1.0,e18v)

   e37v=Tbsim[6]/Teff37v
   e37v=max(0.0,e37v)
   e37v=min(1.0,e37v)

   e89v=Tbsim[8]/Teff89v
   e89v=max(0.0,e89v)
   e89v=min(1.0,e89v)
   ev=np.array([e6v,e10v,e18v,e18v,e37v,e37v,e37v,e89v])

   return ev


def eh_ice(Tb,IST,SD):
   Tb6v = np.float(Tb[0]) #6v
   Tb6h = np.float(Tb[1]) #6h
   Tb10v = np.float(Tb[2]) #10v
   Tb10h = np.float(Tb[3]) #10h
   Tb19v = np.float(Tb[4]) #19v
   Tb19h = np.float(Tb[5]) #19h
   Tb23v = np.float(Tb[6]) #23v
   Tb23h = np.float(Tb[7]) #23h
   Tb37v = np.float(Tb[8]) #37v
   Tb37h = np.float(Tb[9]) #37h
   Tb89v = np.float(Tb[10]) #89v
   Tb89h = np.float(Tb[11]) #89h

   GR0610v=(Tb10v-Tb6v)/(Tb10v+Tb6v)
   GR0618v=(Tb19v-Tb6v)/(Tb19v+Tb6v)
   GR0618h=(Tb19h-Tb6h)/(Tb19h+Tb6h)
   GR0636v=(Tb37v-Tb6v)/(Tb37v+Tb6v)

   #Lise's equations
   SDsim = 1.7701 + 0.017462*Tb6v - 0.02801*Tb19v + 0.0040926*Tb37v

   Tsi6 = 1.144*Tb6v - 0.815/SDsim - 27.08
   Tsi10 = 1.101*Tb10v - 0.999/SDsim - 14.28
   Tsi = 0.5*(Tsi6 + Tsi10)

   Teff6v  = 0.888  * Tsi + 30.245
   Teff10v = 0.901  * Tsi + 26.569
   Teff19v = 0.920  * Tsi + 21.536
   Teff23v = 0.932  * Tsi + 18.417
   Teff37v = 0.960  * Tsi + 10.902
   Teff50v = 0.989  * Tsi + 2.959
   Teff89v = 1.0604 * Tsi - 16.384

   #Rasmus equation for ice thickness
   ITsim= 4.81691061767 - 47.3872663622*GR0610v - 52.8250933515*GR0618v + 49.9619055273*GR0618h + 33.4020392591*GR0636v + 0.142250481385*Tb6h - 0.115797387369*Tb19h - 0.0385743856297*Tb37v
   ITsim=max(0.3,ITsim)
   ITsim=min(3.5,ITsim)

   #Rasmus forward model
   Tbsim=rm.reg_mod13(IST,SD,ITsim)

   e6h=Tbsim[1]/Teff6v
   e6h=max(0.0,e6h)
   e6h=min(1.0,e6h)

   e10h=Tbsim[3]/Teff10v
   e10h=max(0.0,e10h)
   e10h=min(1.0,e10h)

   e18h=Tbsim[5]/Teff19v
   e18h=max(0.0,e18h)
   e18h=min(1.0,e18h)

   e37h=Tbsim[7]/Teff37v
   e37h=max(0.0,e37h)
   e37h=min(1.0,e37h)

   e89h=Tbsim[9]/Teff89v
   e89h=max(0.0,e89h)
   e89h=min(1.0,e89h)
   eh=np.array([e6h,e10h,e18h,e18h,e37h,e37h,e37h,e89h])

   return eh


def icemod(Tb,IST,SD):
   Tb6v = np.float(Tb[0]) #6v
   Tb6h = np.float(Tb[1]) #6h
   Tb10v = np.float(Tb[2]) #10v
   Tb10h = np.float(Tb[3]) #10h
   Tb19v = np.float(Tb[4]) #19v
   Tb19h = np.float(Tb[5]) #19h
   Tb23v = np.float(Tb[6]) #23v
   Tb23h = np.float(Tb[7]) #23h
   Tb37v = np.float(Tb[8]) #37v
   Tb37h = np.float(Tb[9]) #37h
   Tb89v = np.float(Tb[10]) #89v
   Tb89h = np.float(Tb[11]) #89h
 
   GR0610v=(Tb10v-Tb6v)/(Tb10v+Tb6v)
   GR0618v=(Tb19v-Tb6v)/(Tb19v+Tb6v)
   GR0618h=(Tb19h-Tb6h)/(Tb19h+Tb6h)
   GR0636v=(Tb37v-Tb6v)/(Tb37v+Tb6v)

   ITsim= 4.81691061767 - 47.3872663622*GR0610v - 52.8250933515*GR0618v + 49.9619055273*GR0618h + 33.4020392591*GR0636v + 0.142250481385*Tb6h - 0.115797387369*Tb19h - 0.0385743856297*Tb37v
   ITsim=max(0.3,ITsim)
   ITsim=min(3.5,ITsim)
   SDsim = 1.7701 + 0.017462*Tb6v - 0.02801*Tb19v + 0.0040926*Tb37v
   Tsi6 = 1.144*Tb6v - 0.815/SDsim - 27.08
   Tsi10 = 1.101*Tb10v - 0.999/SDsim - 14.28 
   Tsi = 0.5*(Tsi6 + Tsi10)
   Tsi = np.sqrt(Tsi6*Tsi10)
   Tsi=max(200.0,Tsi)
   Tsi=min(273.15,Tsi)

   Teff6v  = 0.888  * Tsi + 30.245
   Teff10v = 0.901  * Tsi + 26.569
   Teff19v = 0.920  * Tsi + 21.536
   Teff23v = 0.932  * Tsi + 18.417
   Teff37v = 0.960  * Tsi + 10.902
   Teff50v = 0.989  * Tsi + 2.959
   Teff89v = 1.0604 * Tsi - 16.384
   Teff=np.array([Teff6v,Teff10v,Teff19v,Teff23v,Teff37v,Teff50v,Teff89v])

   Tbsim=rm.reg_mod13(IST,SD,ITsim)

   e6v=Tbsim[0]/Teff6v
   e6v=max(0.0,e6v)
   e6v=min(1.0,e6v)

   e6h=Tbsim[1]/Teff6v
   e6h=max(0.0,e6h)
   e6h=min(1.0,e6h)

   e10v=Tbsim[2]/Teff10v
   e10v=max(0.0,e10v)
   e10v=min(1.0,e10v)

   e10h=Tbsim[3]/Teff10v
   e10h=max(0.0,e10h)
   e10h=min(1.0,e10h)

   e18v=Tbsim[4]/Teff19v
   e18v=max(0.0,e18v)
   e18v=min(1.0,e18v)

   e18h=Tbsim[5]/Teff18v
   e18h=max(0.0,e18h)
   e18h=max(1.0,e18h)

   e37v=Tbsim[6]/Teff37v
   e37v=max(0.0,e37v)
   e37v=min(1.0,e37v)

   e37h=Tbsim[7]/Teff37v
   e37h=max(0.0,e37h)
   e37h=min(1.0,e37h)

   e89v=Tbsim[8]/Teff89v
   e89v=max(0.0,e89v)
   e89v=min(1.0,e89v)

   e89h=Tbsim[9]/Teff89v
   e89h=max(0.0,e89h)
   e89h=min(1.0,e89h)

   e=np.array([e6v, e6h, e10v, e10h, e18v, e18h, e37v, e37h, e89v, e89h])

   return IST, SD, Tbsim, e, Teff, Tsi
