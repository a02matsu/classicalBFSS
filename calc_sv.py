#! /usr/local/intelpython3/bin/python
import numpy as np
import matplotlib.pylab as plt
import subprocess
import csv

NMAT=4
DIM=9
matrix_size = DIM*NMAT*NMAT
##############
#XVF="X"
ver="0001"
MASS="0.0"
delay_time="0.0"  # delay time ("t" in the note)
Delta="0.5"
##############
for XVF in ["X","V","F"]:
    INPUT_FILE="OUTPUT/M" + XVF + XVF + "_M" + MASS + "t" + delay_time + "D" + Delta + "_" + ver
    SV_FILE="SV/sv" + "_" + XVF + "_M" + MASS + "t" + delay_time + "D" + Delta + "_" + ver

    # Count lines (the first line is a comment)
    num_lines = sum(1 for line in open(INPUT_FILE) )
    #num_delay = int( float(delay_time) / deltaT )
    #num_Delta = int( float(Delta) / deltaT ) 
    #num_data = int( ( num_lines - num_delay) / num_Delta )

    ############
    with open( INPUT_FILE, 'r' ) as matrices:
        counter=0
        vec=np.zeros(matrix_size*matrix_size,dtype=complex)
        sv=np.zeros((num_lines-1)*matrix_size,dtype=float).reshape(num_lines-1, matrix_size)

#tmp1=np.zeros((2*matrix_size)*(2*matrix_size),dtype=float).reshape(num_data, (2*matrix_size))

        n=0
        for line in matrices.readlines()[1:]:
            counter = 0
            parity = 0
            #print( len(np.array(line.split(), dtype=float)[1:]))
            for ele in np.array(line.split(), dtype=float)[1:]:
                if parity == 0:
                    real = ele
                    parity = 1
                else:
                    vec[counter] = real + ele*1.j
                    #print( "::",counter, vec[counter] )
                    counter += 1
                    parity = 0
            #print( "::", vec[0], vec[matrix_size-1], vec[2*matrix_size-2] )
            #print( "-----------" )
    
            matrix = vec.reshape(matrix_size, matrix_size)
            #print( matrix )
            #print( "-----------" )

            #print( matrix[0,0], matrix[1,1], matrix[2,2], matrix[-1,-1] )
            sv[n] = np.linalg.svd(matrix, full_matrices=True, compute_uv=False)
            n += 1
        

    np.savetxt( SV_FILE, sv )


