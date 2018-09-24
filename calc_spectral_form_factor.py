#! /usr/local/intelpython3/bin/python
import numpy as np
import matplotlib.pylab as plt
import subprocess
import csv

##############
XVF="X"
ver="0001"
t="0.0"
D="0.1"
##############
Ndiv=2000 # number of division
ini_tau=0E0
fin_tau=1E14
Dtau = (fin_tau - ini_tau) / float(Ndiv)
##############
ini_sample=1 # which SV 
num_samples=1
##############


ALL_SV_FILE="EIGENS/" + XVF + "MD_t" + t +"D" + D + "_SV_" + ver
#ALL_SV = np.loadtxt(ALL_SV_FILE)

sv_file = open( ALL_SV_FILE, 'r' ) 
#lines = sv_file.readlines()[ini_sample: ini_sample+num_samples]
#num_svs = len( lines[0].split() )
#print( num_svs, float(lines[0].split()[-1] ) )

Xtau = []
SFF_sample = []
nline=0 # runs from 0 to num_samples-1
for line in sv_file.readlines()[ini_sample: ini_sample+num_samples]:
    tmp = 0. + 0.j
    ntau=0 # runs from  0 to Ndiv
    for tau in np.arange(ini_tau, fin_tau, Dtau):
        Xtau.append( tau )
        SFF_sample.append( [] )
        n=0
        for sv in line.split(): 
            tmp += np.exp( 1.j * float(sv) * tau )
            n += 1
        SFF_sample[ntau].append( (np.abs(tmp) / float(n))**2  )
        ntau += 1
    #print( SFF_sample[nline] )
    nline += 1

#print( SFF_sample )
SFF=[]
for spl in SFF_sample:
    SFF.append( np.average( spl ) )

np.savetxt( "test.dat", SFF )


plt.xscale("log")
plt.yscale("log")
plt.plot(Xtau, SFF)
plt.show()


    




