#!/usr/bin/python

import csv
import numpy as np
import optimale
import time

def runows(ncf):
    tic =time.clock()
      
    i=0
    for row_met in met:
	#print row_met[:]
	#including the atmosphere:
        #amsr_channels
        Tb6v = ncf['tb6v'][:].flatten()
        #Tb6v = np.float(row_met[0]) #6v
        Tb6h = np.float(row_met[1]) #6h
        Tb10v = np.float(row_met[2]) #10v
        Tb10h = np.float(row_met[3]) #10h
        Tb19v = np.float(row_met[4]) #19v
        Tb19h = np.float(row_met[5]) #19h
        Tb23v = np.float(row_met[6]) #23v
        Tb23h = np.float(row_met[7]) #23h
        Tb37v = np.float(row_met[8]) #37v
        Tb37h = np.float(row_met[9]) #37h
        Tb89v = np.float(row_met[10]) #89v
        Tb89h = np.float(row_met[11]) #89h

        Tb = np.array([Tb6v,Tb6h,Tb10v,Tb10h,Tb19v,Tb19h,Tb23v,Tb23h,Tb37v,Tb37h,Tb89v,Tb89h])
        P_final, S_diag, cost = optimale.optimal(Tb)
        #print P_final 
        #the string to be printed to the file
        oic = str(P_final[0])+" "+str(S_diag[0])+" "+str(P_final[1])+" "+str(S_diag[1])+" "\
             +str(P_final[2])+" "+str(S_diag[2])+" "+str(P_final[3])+" "+str(S_diag[3])+" "\
             +str(P_final[4])+" "+str(S_diag[4])+" "+str(P_final[5])+" "+str(S_diag[5])+" "+str(cost)+"\n"
        ofile.write(oic)
        i=i+1
    ifile.close()
    ofile.close()
    toc=time.clock()
    print toc-tic

#ifile = open('AMSR_SIC10_2008.txt','rt')
#ofile = open('OE_AMSR_SIC0_2008.txt','a')
#met = csv.reader(ifile, delimiter = ' ')
met = f.variables
runows(met)

#def runice():
    #tic =time.clock()
    #ifile = open('AMSR_SIC1_2008.txt','rt')
    #ofile = open('OE_AMSR_icemodel.txt','a')
    #met = csv.reader(ifile, delimiter = ' ')
    #i=0
    #for row_met in met:
	##print row_met[:]
	##including the atmosphere:
        ##amsr_channels
        #Tb6v = np.float(row_met[0]) #6v
        #Tb6h = np.float(row_met[1]) #6h
        #Tb10v = np.float(row_met[2]) #10v
        #Tb10h = np.float(row_met[3]) #10h
        #Tb19v = np.float(row_met[4]) #19v
        #Tb19h = np.float(row_met[5]) #19h
        #Tb23v = np.float(row_met[6]) #23v
        #Tb23h = np.float(row_met[7]) #23h
        #Tb37v = np.float(row_met[8]) #37v
        #Tb37h = np.float(row_met[9]) #37h
        #Tb89v = np.float(row_met[10]) #89v
        #Tb89h = np.float(row_met[11]) #89h

        #Tb = np.array([Tb6v,Tb6h,Tb10v,Tb10h,Tb19v,Tb19h,Tb23v,Tb23h,Tb37v,Tb37h,Tb89v,Tb89h])
        #P_final, S_diag = optimale.optimal(Tb)
        #print P_final 
        ##the string to be printed to the file
        #oic = str(P_final[0])+" "+str(S_diag[0])+" "+str(P_final[1])+" "+str(S_diag[1])+" "\
             #+str(P_final[2])+" "+str(S_diag[2])+" "+str(P_final[3])+" "+str(S_diag[3])+"\n"
        #ofile.write(oic)
        #i=i+1
    #ifile.close()
    #ofile.close()
    #toc=time.clock()
    #print toc-tic

#def runemi():
    #tic =time.clock()
    #ifile = open('AMSR_SIC1_2008.txt','rt')
    #ofile = open('OE_AMSR_emis.txt','a')
    #met = csv.reader(ifile, delimiter = ' ')
    #i=0
    #for row_met in met:
	##print row_met[:]
	##including the atmosphere:
        ##amsr_channels
        #Tb6v = np.float(row_met[0]) #6v
        #Tb6h = np.float(row_met[1]) #6h
        #Tb10v = np.float(row_met[2]) #10v
        #Tb10h = np.float(row_met[3]) #10h
        #Tb19v = np.float(row_met[4]) #19v
        #Tb19h = np.float(row_met[5]) #19h
        #Tb23v = np.float(row_met[6]) #23v
        #Tb23h = np.float(row_met[7]) #23h
        #Tb37v = np.float(row_met[8]) #37v
        #Tb37h = np.float(row_met[9]) #37h
        #Tb89v = np.float(row_met[10]) #89v
        #Tb89h = np.float(row_met[11]) #89h

        #Tb = np.array([Tb6v,Tb6h,Tb10v,Tb10h,Tb19v,Tb19h,Tb23v,Tb23h,Tb37v,Tb37h,Tb89v,Tb89h])
        #e = optimale.optimal_e50(Tb)
        #print e[0],e[1] 
        ##the string to be printed to the file
        #oic = str(e[0])+" "+str(e[1])+" "+str(e[2])+" "+str(e[3])+" "\
             #+str(e[4])+" "+str(e[5])+"\n"
        #ofile.write(oic)
        #i=i+1
    #ifile.close()
    #ofile.close()
    #toc=time.clock()
    #print toc-tic
