"""
Main routine for solving a stationary 1D Schroedinger equation
"""

import reader
import interpolate as ip
import solver
import plotter
import numpy as np
#Deleted a bunch of useless commands and maybe some usefull here but simple ones for merge
#and to make everything simple and understandable again


def main():
    #get data
    mass, xmin, xmax, xnum, first_val, last_val, inttype, numinterpol, pot, inputpath = reader.getdata()

    #interpolate data
    int_fct = ip.interpol(pot[:, 0], pot[:, 1], inttype)
    potential = ip.potential_grid(xmin, xmax, xnum, int_fct)

    hamiltonian = solver.hamilton_operator(potential, mass, xnum)

    #find Eigenvalues and Eigenvectors
    energies, eigenvec = (solver.diag_solver(hamiltonian, xnum, first_val, last_val))

    #Norming
    wavefuncs = np.zeros((xnum-2, last_val - first_val + 3))
    print(wavefuncs, np.shape(wavefuncs))
    expvalues = np.zeros((last_val - first_val + 2, 2))
    wavefuncs[:,0] = potential[1:xnum - 1, 0]
    for ii in range (0,last_val - first_val + 2):
        wavefuncs[:,ii+1] = solver.QM_Norming(eigenvec[:,ii], potential[1:xnum - 1, 0])
        print(solver.Integrate(wavefuncs[:,ii]**2,potential[1:xnum - 1, 0]))
        expvalues[ii,:] = solver.QM_position_info(wavefuncs[:,ii], potential[1:xnum - 1, 0])
    print(expvalues)
    #Expectationvalue and Variance
    #Save calculated data
    reader.savedata(inputpath, potential, energies, wavefuncs, expvalues)

if __name__ == '__main__':
    main()
