'''
Module reads data
'''

import os.path
import numpy as np
#import scipy



def getdata(path):
    '''
    Reads data from input file
    '''
    inputpath = os.path.join("c:/", "users/", "boett/", "onedrive/", "dokumente/", "github/",   path, "schrodinger.inp") #pylint: disable=line-too-long
    #data = open(inputpath, "r")
    with open(inputpath, "r") as data:
        masse = data.readline()
        masse = float(masse.replace(' # mass',''))

        minmax = data.readline()
        minmax = minmax.replace('# xMin xMax nPoint','')
        minmax = minmax.split()
        xmin, xmax, xnum = minmax
        xmin = float(xmin)
        xmax = float(xmax)
        xnum = float(xnum)
        firstlast = data.readline()
        firstlast = firstlast.replace('# first and last eigenvalue to print','')
        first, last = firstlast.split()

        interpolation = data.readline()
        interpolation = interpolation.replace('# interpolation type','')
        numberinterpol = data.readline()
        numberinterpol = numberinterpol.replace('# nr. of interpolation points and xy declarations','') #pylint: disable=line-too-long
        rest = np.loadtxt(data.readlines())

    return xmin, xmax, xnum, masse, rest, interpolation, numberinterpol, first, last

print(getdata("sciproject"))
