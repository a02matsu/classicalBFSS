#! /usr/local/intelpython3/bin/python
import numpy as np

M="0.0"
t="0.0"
D="0.5"
ver="0001"

XFILE="SV/sv_X_M" + M + "t" + t +"D" + D + "_" + ver
VFILE="SV/sv_V_M" + M + "t" + t +"D" + D + "_" + ver
FFILE="SV/sv_F_M" + M + "t" + t +"D" + D + "_" + ver

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




