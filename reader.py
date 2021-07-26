'''
Module reads data
'''
import time
start_time = time.time()

import os.path
import numpy as np

def getdata(path):
    '''
    Reads data from input file
    '''
    inputpath = os.path.join("c:/", "users/", "boett/", "onedrive/", "dokumente/", "github/",   path, "schrodinger.inp") #pylint: disable=line-too-long
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
        first= float(first)
        last = float(last)
        interpolation = data.readline()
        interpolation = interpolation.replace('# interpolation type','')
        interpolation = interpolation.strip()
        numinterpol = data.readline()
        numinterpol = numinterpol.replace('# nr. of interpolation points and xy declarations','')
        numinterpol = numinterpol.strip()
        numinterpol = float(numinterpol)
        rest = np.loadtxt(data.readlines())

    return masse, xmin, xmax, xnum, first, last, interpolation, numinterpol, rest

print(getdata("sciproject"))
print("--- %s seconds ---" % (time.time() - start_time))
