'''
Module that reads data
'''
import os.path
import numpy as np




def getdata():
    '''
    Reads data from input file
    Args:
        None yet
    Returns:
        masse:           mass of the particle as float
        xmin:           first considered x-value for interpolation as float
        xmax:           last considered x-value for interpolation as float
        xnum:           Number of Points used for interpolation
        first, last:    first and last desired eigenvalue as int
        interpolation:  Desired interpolation type as string
        numinterpol:    Number of given x-values with corresponding potential
        potential:      2D Array of the Potential ([:,1]) and corresponding
                        x-values([:,0])
    '''
    try: 
        inputpath = "schroodinger.inp"
        open(inputpath)  
    except: 
        inputpath = input("give data: pls give path danke schön wunder bar ")
        

    with open(inputpath, "r") as data:
        masse = data.readline()
        masse = float(masse.replace(' # mass', ''))

        minmax = data.readline()
        minmax = minmax.replace('# xMin xMax nPoint', '')
        minmax = minmax.split()
        xmin, xmax, xnum = minmax
        xmin = float(xmin)
        xmax = float(xmax)
        xnum = int(xnum)
        firstlast = data.readline()
        firstlast = firstlast.replace('# first and last eigenvalue to print', '')
        first, last = firstlast.split()
        first = int(first)
        last = int(last)
        interpolation = data.readline()
        interpolation = interpolation.replace('# interpolation type', '')
        interpolation = interpolation.strip()
        numinterpol = data.readline()
        numinterpol = numinterpol.replace('# nr. of interpolation points and xy declarations', '')
        numinterpol = numinterpol.strip()
        numinterpol = float(numinterpol)
        potential = np.loadtxt(data.readlines())

    return masse, xmin, xmax, xnum, first, last, interpolation, numinterpol, potential

masse, xmin, xmax, xnum, first, last, interpolation, numinterpol, potential =   getdata()
print(masse)


'''

def savedata( ):


    np.savetxt(   , 'potential.dat'    )
    np.savetxt(   , 'energies.dat'    )
    np.savetxt(   , 'wavefuncs.dat'    )
    np.savetxt(   , 'expvalues.dat'    )

'''