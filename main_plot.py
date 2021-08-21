#!/usr/bin/env python3
"""Main routine for plotting the results of the solver"""

from modules import plotter

_DESCRIPTION = "Plots given data from a 1-dimension qm-system."

def main():
    '''
    Main routine for plotting the results of the solver
    '''
    energies, wavefuncs, expvalues, potential, path = plotter.read_data()
    plotter.ploty(energies, wavefuncs, expvalues, potential, path)

if __name__ == '__main__':
    main()
