#! /usr/local/intelpython3/bin/python
import numpy as np

ver="0000"
D="1.0"
t="0.0"

XFILE="SV/sv" + ver + "_Xmat_t" + t +"D" + D 
VFILE="SV/sv" + ver + "_Vmat_t" + t +"D" + D 
FFILE="SV/sv" + ver + "_Fmat_t" + t +"D" + D 

print( XFILE )

Xa = np.loadtxt(XFILE)
Va = np.loadtxt(VFILE)
Fa = np.loadtxt(FFILE)
##
Xmax1 = np.log(Xa[:,0] )
Xmax2 = np.log(Xa[:,1])
##
Vmax1 = np.log(Va[:,0] )
Vmax2 = np.log(Va[:,1])
##
Fmax1 = np.log(Fa[:,0] )
Fmax2 = np.log(Fa[:,1])

np.set_printoptions(precision=4)
print("XX max1:", np.round(np.average(Xmax1),4), "+/-", np.round(np.cov(Xmax1),4) )
print("XX max2:", np.round(np.average(Xmax2),4), "+/-", np.round(np.cov(Xmax2),4) )
print("VV max1:", np.round(np.average(Vmax1),4), "+/-", np.round(np.cov(Vmax1),4) )
print("VV max2:", np.round(np.average(Vmax2),4), "+/-", np.round(np.cov(Vmax2),4) )
print("FF max1:", np.round(np.average(Fmax1),4), "+/-", np.round(np.cov(Fmax1),4) )
print("FF max2:", np.round(np.average(Fmax2),4), "+/-", np.round(np.cov(Fmax2),4) )




