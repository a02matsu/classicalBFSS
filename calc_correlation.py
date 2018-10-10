#! /usr/local/intelpython3/bin/python
import numpy as np
import matplotlib.pylab as plt
import subprocess
import csv

NMAT=4
DIM=9
MASS=0.0
matrix_size = DIM*NMAT*NMAT
##############
XVF="X"
ver="0000"
delay_time="0.0"  # delay time ("t" in the note)
Delta="2.0"  # range of integration

deltaT=1.0  # time step of the input data (in OUTPUT)
##############
INPUT_FILE="OUTPUT/M_" + XVF + XVF + "t" + delay_time + "D" + Delta + "_" + ver
SV_FILE="SV/sv_" + XVF + "_t" + delay_time + "D" + Delta + "_" + ver

# Count lines (the first line is a comment)
num_lines = sum(1 for line in open(INPUT_FILE) )
num_delay = int( float(delay_time) / deltaT )
num_Delta = int( float(Delta) / deltaT ) 
num_data = int( ( num_lines - num_delay) / num_Delta )

#print( num_lines, num_delay, num_Delta, num_data)

############
matrices = open( INPUT_FILE, 'r' ) 
n=0
counter=0
tmp1=np.zeros(num_data*(2*matrix_size),dtype=float).reshape(num_data, (2*matrix_size))
for line in matrices.readlines()[1:]:
    counter += 1
    tmp1[n] += np.array(line.split(), dtype=float)[1:] * deltaT
    if counter == num_Delta:
        counter = 0
        n += 1
matrices.close()
############
matrices = open( INPUT_FILE, 'r' ) 
n=0
counter=0
skip=0
tmp2=np.zeros(num_data*2*matrix_size,dtype=float).reshape(num_data, 2*matrix_size)
for line in matrices.readlines()[1:]:
    skip += 1
    if skip > num_delay:
        counter += 1
        tmp2[n] += np.array(line.split(), dtype=float)[1:] * deltaT
        if counter == num_Delta:
            counter = 0
            n += 1
matrices.close()

ave1=np.zeros(num_data*matrix_size, dtype=complex).reshape(num_data, matrix_size)
ave2=np.zeros(num_data*matrix_size, dtype=complex).reshape(num_data, matrix_size)
for n in np.arange(0,2*matrix_size-1,2):
    ave1[:,int(n/2)] = tmp1[:,n] + tmp1[:,n+1]*1.j
    ave2[:,int(n/2)] = tmp2[:,n] + tmp2[:,n+1]*1.j

print( ave1[0,:], ave1[-1,:] )

matrix=np.zeros(matrix_size*matrix_size, dtype=complex).reshape(matrix_size,matrix_size)
sv_data=np.zeros(matrix_size*num_data, dtype=float).reshape(num_data,matrix_size)
for d in np.arange(0,num_data-1):
    for i in np.arange(0,matrix_size):
        for j in np.arange(0,matrix_size):
            matrix[i,j] = np.conj(ave1[d,i]) * ave2[d,j]
    sv_data[d,:] = np.linalg.svd(matrix, full_matrices=True, compute_uv=False)

header="# NMAT= " + str(NMAT) + " M= ", str(MASS), " Delta= ", Delta, " t= ", delay_time
#np.savetxt( SV_FILE, header )
np.savetxt( SV_FILE, sv_data )


