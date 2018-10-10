#! /usr/local/intelpython3/bin/python
import numpy as np
import matplotlib.pylab as plt
import subprocess
import csv
from scipy import optimize


##############
XVF="X"
MASS="0.0"
t="2.0" # "t" in the note
D="0.0" # "\Delta" in the note
ver="0000"
##############
Ndiv=500 # number of division
ini_tau=1E-1
fin_tau=1E3
Dlogtau = (np.log(fin_tau) - np.log(ini_tau)) / float(Ndiv)

##############
ini_sample=1 # which SV 
num_samples=1
##############


SV_FILE="SV/sv_" + XVF + "_M" + MASS + "t" + t +"D" + D + "_" + ver
SFF_FILE="SFF_" + XVF + "_M" + MASS + "t" + t +"D" + D + "_" + ver
num_lines = sum(1 for line in open(SV_FILE) )

SFF_sample=np.zeros(Ndiv*num_lines,dtype=float).reshape(Ndiv,num_lines) # [ [sample1, sample2, ...., sampleN], tau=ini..fin]
with open( SV_FILE, 'r' ) as sv_data:
    nline=-1
    for line in sv_data.readlines():
        nline += 1
        ntau=-1
        for logtau in np.arange(np.log(ini_tau), np.log(fin_tau), Dlogtau):
            ntau += 1
            tmp=0.+0.j
            n=0
            for sv in np.array(line.split(), dtype=float):
                n += 1
                tmp += np.exp( 1.j * np.log(sv) * np.exp(logtau) )
            SFF_sample[ntau,nline] = (np.abs(tmp) / float(n))**2  


SFF = np.zeros(Ndiv*2,dtype=float).reshape(Ndiv,2)
n=0
for logtau in np.arange(np.log(ini_tau), np.log(fin_tau), Dlogtau):
    SFF[n,0]=np.exp(logtau)
    SFF[n,1]=np.average(SFF_sample[n,:])
    n += 1

np.savetxt( SFF_FILE, SFF )

### fitting
def f(x,a,b):
    return a*x**b
parameter_initial=[1., 1.]

x1 = SFF[:,0][ ( SFF[:,0] > 20 ) & ( SFF[:,0] < 60 ) ]
y1 = SFF[:,1][ ( SFF[:,0] > 20 ) & ( SFF[:,0] < 60 ) ]
popt, pcov = optimize.curve_fit(f, x1, y1 , p0=parameter_initial)

y=f(SFF[:,0],*popt)

plt.xscale("log")
plt.yscale("log")
plt.plot(SFF[:,0],SFF[:,1])
#plt.plot(x1,f(x1,*popt), label='fitting')
plt.plot(SFF[:,0],y,label='x^'+f"{popt[1]:.2f}")
plt.legend()
plt.title( "SFF of " + XVF + XVF + ": M=" + MASS + ", t=" + t + ", $\Delta$=" + D )
plt.show()


    




