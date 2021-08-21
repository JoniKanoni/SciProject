'''
Main routine for plotting the results of the solver
'''

from modules import plotter

def main():
    '''
    Main routine for plotting the results of the solver
    '''
    energies, wavefuncs, expvalues, potential = plotter.read_data()
    plotter.ploty(energies, wavefuncs, expvalues, potential)

if __name__ == '__main__':
    main()
