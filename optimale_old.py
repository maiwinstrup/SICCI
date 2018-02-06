#!/usr/bin/python
# -*- coding: ascii -*-

import numpy as np
#import ssmirtm2 as ssmi
#import smmr
import amsrrtm2
#import emimod
import icemod_old as icemod
import bootstrap_f_old as bootstrap_f


def optimal(Tb):
    #input parametre aabent vand :
    #maalte straalingstemperature Tb [K], 1xn vektor
    #kovariansmatricen Se for Tb nxn matrice
    #foerste gaet paa fysiske vaerdier P0, 1xm vektor
    #kovarians matricen Sp for de fysiske vaerdier, mxm matrice

    #kovariansmatricen for de estimerede straalingstemperatur vaerdier S
    Se=np.matrix([[0.5000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
                  [0, 0.5000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
                  [0, 0, 0.5000, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
                  [0, 0, 0, 0.5000, 0, 0, 0, 0, 0, 0, 0, 0],\
                  [0, 0, 0, 0, 0.5000, 0, 0, 0, 0, 0, 0, 0],\
                  [0, 0, 0, 0, 0, 0.5000, 0, 0, 0, 0, 0, 0],\
                  [0, 0, 0, 0, 0, 0, 0.5000, 0, 0, 0, 0, 0],\
                  [0, 0, 0, 0, 0, 0, 0, 0.5000, 0, 0, 0, 0],\
                  [0, 0, 0, 0, 0, 0, 0, 0, 0.5000, 0, 0, 0],\
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5000, 0, 0],\
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000, 0],\
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0000]])

    #kovariansmatricen for de estimerede fysiske vardier P
    Sp=np.matrix([[15.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
                  [0.0, 5.0, 0.0, 0.0, 0.0, 0.0],\
                  [0.0, 0.0, 5.0, 0.0, 0.0, 0.0],\
                  [0.0, 0.0, 0.0, 0.3, 0.0, 0.0],\
                  [0.0, 0.0, 0.0, 0.0, 0.4, 0.0],\
                  [0.0, 0.0, 0.0, 0.0, 0.0, 0.1]])

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

    #Lise's equations for snowdepth
    SDsim = 1.7701 + 0.017462*Tb6v - 0.02801*Tb19v + 0.0040926*Tb37v

    #vector med fysiske vaerdier
    Ts=275.0 #sst, overfladetemperatur, effektiv temperatur
    #Ts_amsr=Ts*np.ones(8) #sst, overfladetemperatur, effektiv temperatur 6-89GHz
    V=3.6 #default water vapour
    W=5.0 #default wind speed
    L=0.08 #default integrated liquid water 
    #c_ice=1.0 is open water and sic=2.0 is ice
    cf=bootstrap_f.bootstrap_f(Tb19v,Tb37v)
    c_ice=cf+1.0 #openwater eller hva? vi mangler en isemissivitetsfunktion
    sd=SDsim
    ev_ice=np.ones(8) #emissivitet vertikal polarisation 6-89GHz over is!
    eh_ice=np.ones(8) #emissivitet horisontal polarisation 6-89GHz over is!
    Ts_amsr=icemod.Ts_amsr(Tb)
    P0 = np.array([Ts,V,W,L,c_ice,sd]) #initial vaerdier for aabentvandsatmosfaeren 
    #print icemod.ev_ice(Tb,P0[0],P0[5])

    #eksempel paa kald af amsr-funktionen
    #amsr_Tb = amsrrtm2.amsr(V, W, L, Ti, Ts, Ts_amsr, Ts_amsr, c_ice, ev_ice, eh_ice)
    #Tb = np.array([178.0, 111.5, 183.8, 117.3, 204.9,\
    #               146.0, 220.6, 175.9, 225.3, 183.7, 204.3, 139.7])
    #Tb = np.array([157.0, 77.5, 163.8, 82.3, 178.9,\
    #               97.0, 188.6, 111.9, 204.3, 130.7, 238.3, 179.7])

    #oevre og nedre graenser til de fysiske parametre
    L=np.array([200.0,300.0,0.0,40.0,0.0,50.0,0.0,1.0,0.9,2.1,0.0,1.0])
    len_p = len(P0)
    len_tb = len(Tb)
    melements = (len_p,len_tb)
    M = np.matrix(np.ones(melements))
    #P0 = np.array([Ts,V,W,L,c_ice,sd])
    n=0
    #start while loop
    while n<7:
        n=n+1
        #kald forward model med seneste gaet. samtidig reference til adjoint model.
        if n==1:
            Ta=amsrrtm2.amsr(P0[1],P0[2],P0[3],P0[0],P0[0],Ts_amsr,Ts_amsr,P0[4],icemod.ev_ice(Tb,P0[0],P0[5]),icemod.eh_ice(Tb,P0[0],P0[5]))
        else:
            Ta=amsrrtm2.amsr(P[1],P[2],P[3],P[0],P[0],Ts_amsr,Ts_amsr,P[4],icemod.ev_ice(Tb,P[0],P[5]),icemod.eh_ice(Tb,P[0],P[5]))
        #start for loop over parametre
        for i in range(len_p):
            #beregn partiell afledede ->kald forwardmodel med pertubationer for hver parameter, pertubationer 1%
            if n==1:
                c_ice=P0[4]
                #M matricen "the adjoint" foeste iteration bruger initialvaerdierne
                M[i,:] = (amsrrtm2.amsr((0.01*(i==1)*P0[1])+P0[1],(0.01*(i==2)*P0[2])+P0[2],(0.01*(i==3)*P0[3])+\
                         P0[3],(0.01*(i==0)*P0[0])+P0[0],(0.01*(i==0)*P0[0])+P0[0],Ts_amsr,Ts_amsr,(0.01*(i==4)*c_ice+c_ice),\
                         icemod.ev_ice(Tb,(0.01*(i==0)*P0[0])+P0[0],(0.01*(i==5)*P0[5])+P0[5]),icemod.eh_ice(Tb,(0.01*(i==0)*P0[0])+P0[0],(0.01*(i==5)*P0[5])+P0[5])) - \
                         amsrrtm2.amsr(P0[1],P0[2],P0[3],P0[0],P0[0],Ts_amsr,Ts_amsr,c_ice,icemod.ev_ice(Tb,P0[0],P0[5]),icemod.eh_ice(Tb,P0[0],P0[5]))) / \
                         (((i==0)*0.01)*P0[0]+((i==1)*0.01)*P0[1]+((i==2)*0.01)*P0[2]+((i==3)*0.01)*P0[3]+((i==4)*0.01)*P0[4]+((i==5)*0.01)*P0[5])
                #efterfoelgende iterationer
            else:
                c_ice=P[4]
                #M matricen "the adjoint" paa iterationspunktet
                M[i,:] = (amsrrtm2.amsr((i==1)*P[1]*0.01+P[1],(i==2)*P[2]*0.01+P[2],(i==3)*P[3]*0.01+\
                         P[3],(i==0)*P[0]*0.01+P[0],(i==0)*P[0]*0.01+P[0],Ts_amsr,Ts_amsr,(0.01*(i==4)*c_ice+c_ice),icemod.ev_ice(Tb,(0.01*(i==0)*P[0])+\
                         P[0],(0.01*(i==5)*P[5])+P[5]),icemod.eh_ice(Tb,(0.01*(i==0)*P[0])+P[0],(0.01*(i==5)*P[5])+P[5])) - \
                         amsrrtm2.amsr(P[1],P[2],P[3],P[0],P[0],Ts_amsr,Ts_amsr,c_ice,icemod.ev_ice(Tb,P[0],P[5]),icemod.eh_ice(Tb,P[0],P[5]))) / \
                         (((i==0)*0.01)*P[0]+((i==1)*0.01)*P[1]+((i==2)*0.01)*P[2]+((i==3)*0.01)*P[3]+((i==4)*0.01)*P[4]+((i==5)*0.01)*P[5])
        M = M.T
        #beregn S
        S = (Sp.I + M.T*Se.I*M).I
        #check konvergens kriteria f.eks. n>5
        #hvis n>4 (antallet af iterationer) print S diagonal + alt andet
        if n>6:
            P_final = P
            S_diag = np.diag(S)
            cost=np.sqrt(sum((Tb-Ta)**2))
        elif n==1:
            P0 = np.matrix(P0).T
            Tb = np.matrix(Tb).T
            Ta = np.matrix(Ta).T
            #estimer fysiske vaerdier ud fra foerste gaet
            P = P0 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P0))
            #Det virker strengt taget ikke til at vaere noedvendigt med fysiske begraensninger
            if (P[0] < L[0]): P[0] = L[0]+0.1
            if (P[0] > L[1]): P[0] = L[1]
            if (P[1] < L[2]): P[1] = L[2]+0.1
            if (P[1] > L[3]): P[1] = L[3]
            if (P[2] < L[4]): P[2] = L[4]+0.1
            if (P[2] > L[5]): P[2] = L[5]
            if (P[3] < L[6]): P[3] = L[6]+0.01
            if (P[3] > L[7]): P[3] = L[7]
            if (P[4] < L[8]): P[4] = L[8]+0.1
            if (P[4] > L[9]): P[4] = L[9]
            if (P[5] < L[10]): P[5] = L[10]+0.1
            if (P[5] > L[11]): P[5] = L[11]
            P_1 = P
            M = M.T
            P0 = np.array(P0.T).reshape(len_p)
            Tb = np.array(Tb.T).reshape(len_tb)
            Ta = np.array(Ta.T).reshape(len_tb)
            P = np.array(P.T).reshape(len_p)
        else:
            P0 = np.matrix(P0).T
            Tb = np.matrix(Tb).T
            Ta = np.matrix(Ta).T
            #2. 3. 4. iteration frem mod fysiske vaerdier der faar Tb til at passe
            P = P_1 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P_1))
            if (P[0] < L[0]): P[0] = L[0]+0.1
            if (P[0] > L[1]): P[0] = L[1]
            if (P[1] < L[2]): P[1] = L[2]+0.1
            if (P[1] > L[3]): P[1] = L[3]
            if (P[2] < L[4]): P[2] = L[4]+0.1
            if (P[2] > L[5]): P[2] = L[5]
            if (P[3] < L[6]): P[3] = L[6]+0.01
            if (P[3] > L[7]): P[3] = L[7]
            if (P[4] < L[8]): P[4] = L[8]+0.1
            if (P[4] > L[9]): P[4] = L[9]
            if (P[5] < L[10]): P[5] = L[10]+0.1
            if (P[5] > L[11]): P[5] = L[11]
            P_1 = P
            M = M.T
            P0 = np.array(P0.T).reshape(len_p)
            Tb = np.array(Tb.T).reshape(len_tb)
            Ta = np.array(Ta.T).reshape(len_tb)
            P = np.array(P.T).reshape(len_p)
            #tilbage til while loop start
    P_final[4]=P_final[4]-1.0
    return P_final, S_diag, cost


#def optimale(Tb):
    ##input parametre aabent vand :
    ##maalte straalingstemperature Tb [K], 1xn vektor
    ##kovariansmatricen Se for Tb nxn matrice
    ##foerste gaet paa fysiske vaerdier P0, 1xm vektor
    ##kovarians matricen Sp for de fysiske vaerdier, mxm matrice

    ##kovariansmatricen for de estimerede straalingstemperatur vaerdier S
    #Se=np.matrix([[0.0900, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
                  #[0, 0.1089, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
                  #[0, 0, 0.2209, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
                  #[0, 0, 0, 0.2916, 0, 0, 0, 0, 0, 0, 0, 0],\
                  #[0, 0, 0, 0, 0.2304, 0, 0, 0, 0, 0, 0, 0],\
                  #[0, 0, 0, 0, 0, 0.2116, 0, 0, 0, 0, 0, 0],\
                  #[0, 0, 0, 0, 0, 0, 0.2025, 0, 0, 0, 0, 0],\
                  #[0, 0, 0, 0, 0, 0, 0, 0.1936, 0, 0, 0, 0],\
                  #[0, 0, 0, 0, 0, 0, 0, 0, 0.2025, 0, 0, 0],\
                  #[0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1600, 0, 0],\
                  #[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1600, 0],\
                  #[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1600]])

    ##kovariansmatricen for de estimerede fysiske vardier P
    #Sp=np.matrix([[5.0, 0.0, 0.0, 0.0],\
                  #[0.0, 5.0, 0.0, 0.0],\
                  #[0.0, 0.0, 5.0, 0.0],\
                  #[0.0, 0.0, 0.0, 5.0]])

    ##vector med fysiske vaerdier
    #Ts=279.5 #sst, overfladetemperatur, effektiv temperatur
    #Ts_amsr=Ts*np.ones(8) #sst, overfladetemperatur, effektiv temperatur 6-89GHz
    #V=3.6 #default water vapour
    #W=5.0 #default wind speed
    #L=0.08 #default integrated liquid water 
    #c_ice=0.0 #openwater eller hva? vi mangler en isemissivitetsfunktion
    #ev_ice=np.ones(8) #emissivitet vertikal polarisation 6-89GHz over is!
    #eh_ice=np.ones(8) #emissivitet horisontal polarisation 6-89GHz over is!
    #P0 = np.array([Ts,V,W,L]) #initial vaerdier for aabentvandsatmosfaeren 

    ##eksempel paa kald af amsr-funktionen
    ##amsr_Tb = amsrrtm2.amsr(V, W, L, Ti, Ts, Ts_amsr, Ts_amsr, c_ice, ev_ice, eh_ice)
    ##Tb = np.array([178.0, 111.5, 183.8, 117.3, 204.9,\
    ##               146.0, 220.6, 175.9, 225.3, 183.7, 204.3, 139.7])
    ##Tb = np.array([157.0, 77.5, 163.8, 82.3, 178.9,\
    ##               97.0, 188.6, 111.9, 204.3, 130.7, 238.3, 179.7])

    ##oevre og nedre graenser til de fysiske parametre
    #L=np.array([270.0,320.0,0.0,40.0,0.0,50.0,0.0,1.0])
    #len_p = len(P0)
    #len_tb = len(Tb)
    #melements = (len_p,len_tb)
    #M = np.matrix(np.ones(melements))
    #n=0
    ##start while loop
    #while n<5:
        #n=n+1
        ##kald forward model med seneste gaet. samtidig reference til adjoint model.
        #if n==1:
            #Ta=amsrrtm2.amsr(P0[1],P0[2],P0[3],P0[0],P0[0],Ts_amsr,Ts_amsr,c_ice,ev_ice,eh_ice)
        #else:
            #Ta=amsrrtm2.amsr(P[1],P[2],P[3],P[0],P[0],Ts_amsr,Ts_amsr,c_ice,ev_ice,eh_ice)
        ##start for loop over parametre
        #for i in range(len_p):
            ##beregn partiell afledede ->kald forwardmodel med pertubationer for hver parameter, pertubationer 1%
            #if n==1:
                ##M matricen "the adjoint" foeste iteration bruger initialvaerdierne
                #M[i,:] = (amsrrtm2.amsr((0.01*(i==1)*P0[1])+P0[1],(0.01*(i==2)*P0[2])+P0[2],(0.01*(i==3)*P0[3])+\
                         #P0[3],(0.01*(i==0)*P0[0])+P0[0],(0.01*(i==0)*P0[0])+P0[0],Ts_amsr,Ts_amsr,c_ice,ev_ice,eh_ice) - \
                         #amsrrtm2.amsr(P0[1],P0[2],P0[3],P0[0],P0[0],Ts_amsr,Ts_amsr,c_ice,ev_ice,eh_ice)) / \
                         #(((i==0)*0.01)*P0[0]+((i==1)*0.01)*P0[1]+((i==2)*0.01)*P0[2]+((i==3)*0.01)*P0[3])
                ##efterfoelgende iterationer
            #else:
                #M[i,:] = (amsrrtm2.amsr((i==1)*P[1]*0.01+P[1],(i==2)*P[2]*0.01+P[2],(i==3)*P[3]*0.01+\
                         #P[3],(i==0)*P[0]*0.01+P[0],(i==0)*P[0]*0.01+P[0],Ts_amsr,Ts_amsr,c_ice,ev_ice,eh_ice) - \
                         #amsrrtm2.amsr(P[1],P[2],P[3],P[0],P[0],Ts_amsr,Ts_amsr,c_ice,ev_ice,eh_ice)) / \
                         #(((i==0)*0.01)*P[0]+((i==1)*0.01)*P[1]+((i==2)*0.01)*P[2]+((i==3)*0.01)*P[3])
        #M = M.T
        ##beregn S
        #S = (Sp.I + M.T*Se.I*M).I
        ##check konvergens kriteria f.eks. n>5
        ##hvis n>4 (antallet af iterationer) print S diagonal + alt andet
        #if n>4:
            #P_final = P
            #S_diag = np.diag(S)
        #elif n==1:
            #P0 = np.matrix(P0).T
            #Tb = np.matrix(Tb).T
            #Ta = np.matrix(Ta).T
            ##estimer fysiske vaerdier ud fra foerste gaet
            #P = P0 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P0))
            ##Det virker strengt taget ikke til at vaere noedvendigt med fysiske begraensninger
            #if (P[0] < L[0]): P[0] = L[0]+0.1
            #if (P[0] > L[1]): P[0] = L[1]
            #if (P[1] < L[2]): P[1] = L[2]+0.1
            #if (P[1] > L[3]): P[1] = L[3]
            #if (P[2] < L[4]): P[2] = L[4]+0.1
            #if (P[2] > L[5]): P[2] = L[5]
            #if (P[3] < L[6]): P[3] = L[6]+0.1
            #if (P[3] > L[7]): P[3] = L[7]
            #P_1 = P
            #M = M.T
            #P0 = np.array(P0.T).reshape(len_p)
            #Tb = np.array(Tb.T).reshape(len_tb)
            #Ta = np.array(Ta.T).reshape(len_tb)
            #P = np.array(P.T).reshape(len_p)
        #else:
            #P0 = np.matrix(P0).T
            #Tb = np.matrix(Tb).T
            #Ta = np.matrix(Ta).T
            ##2. 3. 4. iteration frem mod fysiske vaerdier der faar Tb til at passe
            #P = P_1 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P_1))
            #if (P[0] < L[0]): P[0] = L[0]+0.1
            #if (P[0] > L[1]): P[0] = L[1]
            #if (P[1] < L[2]): P[1] = L[2]+0.1
            #if (P[1] > L[3]): P[1] = L[3]
            #if (P[2] < L[4]): P[2] = L[4]+0.1
            #if (P[2] > L[5]): P[2] = L[5]
            #if (P[3] < L[6]): P[3] = L[6]+0.1
            #if (P[3] > L[7]): P[3] = L[7]
            #P_1 = P
            #M = M.T
            #P0 = np.array(P0.T).reshape(len_p)
            #Tb = np.array(Tb.T).reshape(len_tb)
            #Ta = np.array(Ta.T).reshape(len_tb)
            #P = np.array(P.T).reshape(len_p)
            ##tilbage til while loop start
    #return P_final, S_diag


#def optimal_ice(Tb):
    ##input parametre is :
    ##maalte straalingstemperature Tb [K], 1xn vektor
    ##kovariansmatricen Se for Tb nxn matrice
    ##foerste gaet paa fysiske vaerdier P0, 1xm vektor
    ##kovarians matricen Sp for de fysiske vaerdier, mxm matrice

    ##kovariansmatricen for de estimerede straalingstemperatur vaerdier S
    #Se=np.matrix([[0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0, 0, 0, 0],\
                  #[0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0, 0, 0],\
                  #[0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0, 0],\
                  #[0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0],\
                  #[0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0],\
                  #[0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0],\
                  #[0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0],\
                  #[0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0],\
                  #[0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05],\
                  #[0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1],\
                  #[0, 0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15],\
                  #[0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2]])

    ##kovariansmatricen for de estimerede fysiske vardier P
    #Sp=np.matrix([[0.5, 0.1, 0.1, 0.0],\
                  #[0.1, 0.5, 0.0, 0.0],\
                  #[0.0, 0.0, 15.0, 0.1],\
                  #[0.0, 0.1, 0.1, 0.5]])

    ##vector med fysiske vaerdier
    #r=0.1
    #snow_depth=0.1
    #iceT=261.0
    #myf=0.01
    #P0 = np.array([r,snow_depth,iceT,myf]) #initial vaerdier for isen 

    ##eksempel paa kald af emimod-funktionen
    ##tbv,tbh = emimod.emimodia(ia,freq,snow_rough,ice_rough,snow_corr,snowT,snow_dens,snow_depth,snow_water,iceT,myf)
    ##tbv,tbh = emimod.emimodia(53,6,0.0028,0.005,0.12,snowT,300,snow_depth,0,iceT,myf)
    ##Tb=emimod.emice(snowT,snow_depth,iceT,myf)

    ##oevre og nedre graenser til de fysiske parametre
    #L=np.array([0.0,0.55,0.0,1.0,190.0,273.15,0.0,1.0])
    #len_p = len(P0)
    #len_tb = len(Tb)
    #melements = (len_p,len_tb)
    #M = np.matrix(np.ones(melements))
    #n=0
    ##start while loop
    #while n<5:
        #n=n+1
        ##kald forward model med seneste gaet. samtidig reference til adjoint model.
        #if n==1:
            #Ta=emimod.emice_amsr(P0[0],P0[1],P0[2],P0[3])
        #else:
            #Ta=emimod.emice_amsr(P[0],P[1],P[2],P[3])
        ##start for loop over parametre
        #for i in range(len_p):
            ##beregn partiell afledede ->kald forwardmodel med pertubationer for hver parameter, pertubationer 1%
            #if n==1:
                ##M matricen "the adjoint" foeste iteration bruger initialvaerdierne
                #M[i,:] = (emimod.emice_amsr((0.01*(i==0)*P0[0])+P0[0],(0.01*(i==1)*P0[1])+P0[1],(0.01*(i==2)*P0[2])+P0[2],(0.01*(i==3)*P0[3])+P0[3]) -\
                         #(emimod.emice_amsr(P0[0],P0[1],P0[2],P0[3]))) /\
                         #((0.01*(i==0)*P0[0])+(0.01*(i==1)*P0[1])+(0.01*(i==2)*P0[2])+(0.01*(i==3)*P0[3]))
                ##efterfoelgende iterationer
            #else:
                #M[i,:] = (emimod.emice_amsr((0.01*(i==0)*P[0])+P[0],(0.01*(i==1)*P[1])+P[1],(0.01*(i==2)*P[2])+P[2],(0.01*(i==3)*P[3])+P[3]) -\
                         #(emimod.emice_amsr(P[0],P[1],P[2],P[3]))) /\
                         #((0.01*(i==0)*P[0])+(0.01*(i==1)*P[1])+(0.01*(i==2)*P[2])+(0.01*(i==3)*P[3]))
        #M = M.T
        ##beregn S
        #S = (Sp.I + M.T*Se.I*M).I
        ##check konvergens kriteria f.eks. n>5
        ##hvis n>4 (antallet af iterationer) print S diagonal + alt andet
        #if n>4:
            #P_final = P
            #S_diag = np.diag(S)
        #elif n==1:
            #P0 = np.matrix(P0).T
            #Tb = np.matrix(Tb).T
            #Ta = np.matrix(Ta).T
            ##estimer fysiske vaerdier ud fra foerste gaet
            #P = P0 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P0))
            ##Det virker strengt taget ikke til at vaere noedvendigt med fysiske begraensninger, for aabent vand
            #if (P[0] < L[0]): P[0] = L[0]+0.1
            #if (P[0] > L[1]): P[0] = L[1]
            #if (P[1] < L[2]): P[1] = L[2]+0.01
            #if (P[1] > L[3]): P[1] = L[3]
            #if (P[2] < L[4]): P[2] = L[4]+0.1
            #if (P[2] > L[5]): P[2] = L[5]
            #if (P[3] < L[6]): P[3] = L[6]+0.01
            #if (P[3] > L[7]): P[3] = L[7]
            #P_1 = P
            #M = M.T
            #P0 = np.array(P0.T).reshape(len_p)
            #Tb = np.array(Tb.T).reshape(len_tb)
            #Ta = np.array(Ta.T).reshape(len_tb)
            #P = np.array(P.T).reshape(len_p)
        #else:
            #P0 = np.matrix(P0).T
            #Tb = np.matrix(Tb).T
            #Ta = np.matrix(Ta).T
            ##2. 3. 4. iteration frem mod fysiske vaerdier der faar Tb til at passe
            #P = P_1 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P_1))
            #if (P[0] < L[0]): P[0] = L[0]+0.1
            #if (P[0] > L[1]): P[0] = L[1]
            #if (P[1] < L[2]): P[1] = L[2]+0.01
            #if (P[1] > L[3]): P[1] = L[3]
            #if (P[2] < L[4]): P[2] = L[4]+0.1
            #if (P[2] > L[5]): P[2] = L[5]
            #if (P[3] < L[6]): P[3] = L[6]+0.01
            #if (P[3] > L[7]): P[3] = L[7]
            #P_1 = P
            #M = M.T
            #P0 = np.array(P0.T).reshape(len_p)
            #Tb = np.array(Tb.T).reshape(len_tb)
            #Ta = np.array(Ta.T).reshape(len_tb)
            #P = np.array(P.T).reshape(len_p)
            ##tilbage til while loop start
    #return P_final, S_diag

#def optimal_e50(Tb):
    ##input parametre is :
    ##maalte straalingstemperature Tb [K], 1xn vektor
    ##kovariansmatricen Se for Tb nxn matrice
    ##foerste gaet paa fysiske vaerdier P0, 1xm vektor
    ##kovarians matricen Sp for de fysiske vaerdier, mxm matrice

    ##kovariansmatricen for de estimerede straalingstemperatur vaerdier S
    #Se=np.matrix([[32.203837,       35.812594,       38.760240,       37.464913,       50.689187,       42.778272,       58.008934,       51.726425,       69.792600,       59.173215,       69.817865,       64.864224],\
                  #[35.812594,       59.485133,       44.253562,       62.952884,       60.797946,       68.276922,       69.753896,       75.752940,       86.327074,       84.218645,       73.381330,       70.005190],\
                  #[38.760240,       44.253562,       50.843649,       51.554061,       75.104245,       67.202734,       89.283167,       82.562600,       113.06214,       99.261512,       108.67905,       101.63358],\
                  #[37.464913,       62.952884,       51.554061,       78.223989,       81.229877,       96.116579,       97.199100,       108.24073,       127.37195,       127.40176,       112.05488,       108.50751],\
                  #[50.689187,       60.797946,       75.104245,       81.229877,       129.46084,       122.72271,       160.65676,       153.68886,       215.34336,       194.77497,       198.12912,       186.15222],\
                  #[42.778272,       68.276922,       67.202734,       96.116579,       122.72271,       143.35708,       153.13969,       168.69856,       210.17655,       208.66716,       186.18606,       180.61640],\
                  #[58.008934,       69.753896,       89.283167,       97.199100,       160.65676,       153.13969,       202.37678,       194.20282,       276.19449,       250.64699,       256.11883,       241.03267],\
                  #[51.726425,       75.752940,       82.562600,       108.24073,       153.68886,       168.69856,       194.20282,       206.44278,       268.73769,       260.72077,       245.53194,       236.92715],\
                  #[69.792600,       86.327074,       113.06214,       127.37195,       215.34336,       210.17655,       276.19449,       268.73769,       388.54602,       356.42801,       367.92349,       346.83117],\
                  #[59.173215,       84.218645,       99.261512,       127.40176,       194.77497,       208.66716,       250.64699,       260.72077,       356.42801,       343.49126,       332.44020,       320.22830],\
                  #[69.817865,       73.381330,       108.67905,       112.05488,       198.12912,       186.18606,       256.11883,       245.53194,       367.92349,       332.44020,       440.03604,       414.99232],\
                  #[64.864224,       70.005190,       101.63358,       108.50751,       186.15222,       180.61640,       241.03267,       236.92715,       346.83117,       320.22830,       414.99232,       398.17128]])

    ##kovariansmatricen for de estimerede fysiske vardier P
    #Sp=np.matrix([[5.0, 0.0, 0.0, 0.0],\
                  #[0.0, 5.0, 0.0, 0.0],\
                  #[0.0, 0.0, 5.0, 0.0],\
                  #[0.0, 0.0, 0.0, 5.0]])

    ##vector med fysiske vaerdier
    #snowT=240.0
    #snow_depth=0.1
    #r=0.10
    #myf=0.5
    #P0 = np.array([snowT,snow_depth,r,myf]) #initial vaerdier for isen 

    ##eksempel paa kald af emimod-funktionen
    ##tbv,tbh = emimod.emimodia(ia,freq,snow_rough,ice_rough,snow_corr,snowT,snow_dens,snow_depth,snow_water,iceT,myf)
    ##tbv,tbh = emimod.emimodia(53,6,0.0028,0.005,0.12,snowT,300,snow_depth,0,iceT,myf)
    ##Tb=emimod.emice(snowT,snow_depth,r,myf)

    ##oevre og nedre graenser til de fysiske parametre
    #L=np.array([190.0,273.15,0.0,1.0,0.04,0.9,0.0,1.0])
    #len_p = len(P0)
    #len_tb = len(Tb)
    #melements = (len_p,len_tb)
    #M = np.matrix(np.ones(melements))
    #n=0
    ##start while loop
    #while n<5:
        #n=n+1
        ##kald forward model med seneste gaet. samtidig reference til adjoint model.
        #if n==1:
            #Ta=emimod.emice_e(P0[0],P0[1],P0[2],P0[3])
        #else:
            #Ta=emimod.emice_e(P[0],P[1],P[2],P[3])
        ##start for loop over parametre
        #for i in range(len_p):
            ##beregn partiell afledede ->kald forwardmodel med pertubationer for hver parameter, pertubationer 1%
            #if n==1:
                ##M matricen "the adjoint" foeste iteration bruger initialvaerdierne
                #M[i,:] = (emimod.emice_e((0.01*(i==0)*P0[0])+P0[0],(0.01*(i==1)*P0[1])+P0[1],(0.01*(i==2)*P0[2])+P0[2],(0.01*(i==3)*P0[3])+P0[3]) -\
                         #(emimod.emice_e(P0[0],P0[1],P0[2],P0[3]))) /\
                         #((0.01*(i==0)*P0[0])+(0.01*(i==1)*P0[1])+(0.01*(i==2)*P0[2])+(0.01*(i==3)*P0[3]))
                ##efterfoelgende iterationer
            #else:
                #M[i,:] = (emimod.emice_e((0.01*(i==0)*P[0])+P[0],(0.01*(i==1)*P[1])+P[1],(0.01*(i==2)*P[2])+P[2],(0.01*(i==3)*P[3])+P[3]) -\
                         #(emimod.emice_e(P[0],P[1],P[2],P[3]))) /\
                         #((0.01*(i==0)*P[0])+(0.01*(i==1)*P[1])+(0.01*(i==2)*P[2])+(0.01*(i==3)*P[3]))
        #M = M.T
        ##beregn S
        #S = (Sp.I + M.T*Se.I*M).I
        ##check konvergens kriteria f.eks. n>5
        ##hvis n>4 (antallet af iterationer) print S diagonal + alt andet
        #if n>4:
            #P_final = P
            #S_diag = np.diag(S)
        #elif n==1:
            #P0 = np.matrix(P0).T
            #Tb = np.matrix(Tb).T
            #Ta = np.matrix(Ta).T
            ##estimer fysiske vaerdier ud fra foerste gaet
            #P = P0 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P0))
            ##Det virker strengt taget ikke til at vaere noedvendigt med fysiske begraensninger
            #if (P[0] < L[0]): P[0] = L[0]+0.1
            #if (P[0] > L[1]): P[0] = L[1]
            #if (P[1] < L[2]): P[1] = L[2]+0.01
            #if (P[1] > L[3]): P[1] = L[3]
            #if (P[2] < L[4]): P[2] = L[4]+0.1
            #if (P[2] > L[5]): P[2] = L[5]
            #if (P[3] < L[6]): P[3] = L[6]+0.01
            #if (P[3] > L[7]): P[3] = L[7]
            #P_1 = P
            #M = M.T
            #P0 = np.array(P0.T).reshape(len_p)
            #Tb = np.array(Tb.T).reshape(len_tb)
            #Ta = np.array(Ta.T).reshape(len_tb)
            #P = np.array(P.T).reshape(len_p)
        #else:
            #P0 = np.matrix(P0).T
            #Tb = np.matrix(Tb).T
            #Ta = np.matrix(Ta).T
            ##2. 3. 4. iteration frem mod fysiske vaerdier der faar Tb til at passe
            #P = P_1 + S*(M.T*Se.I*(Tb - Ta) + Sp.I*(P0 - P_1))
            #if (P[0] < L[0]): P[0] = L[0]+0.1
            #if (P[0] > L[1]): P[0] = L[1]
            #if (P[1] < L[2]): P[1] = L[2]+0.01
            #if (P[1] > L[3]): P[1] = L[3]
            #if (P[2] < L[4]): P[2] = L[4]+0.1
            #if (P[2] > L[5]): P[2] = L[5]
            #if (P[3] < L[6]): P[3] = L[6]+0.01
            #if (P[3] > L[7]): P[3] = L[7]
            #P_1 = P
            #M = M.T
            #P0 = np.array(P0.T).reshape(len_p)
            #Tb = np.array(Tb.T).reshape(len_tb)
            #Ta = np.array(Ta.T).reshape(len_tb)
            #P = np.array(P.T).reshape(len_p)
            ##tilbage til while loop start
    ##print P[0],P[1],P[2],P[3]
    #ev,eh=emimod.emice_e50(P[0],P[1],P[2],P[3])
    #e50=np.array([ev,eh,P[0],P[1],P[2],P[3]])
    #return e50
