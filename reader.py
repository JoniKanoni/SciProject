import numpy as np
import scipy 
import os.path





def getdata(path):
    '''
    Reads data from input file
    '''
    
    inputpath = os.path.join("c:/", "users/", "boett/", "onedrive/", "dokumente/", "github/",   path, "schrodinger.inp")
    ei = open(inputpath, "r")

    #for line in ei:
    #    print(line) 
    masse = ei.readline()
    masse = float(masse.replace(' # mass',''))

    minmax = ei.readline()
    minmax = minmax.replace('# xMin xMax nPoint','')
    minmax = minmax.split()
    xmin, xmax, xnum = minmax
    xmin = float(xmin)
    firstlast = ei.readline()

    interpolation = ei.readline()

    numberinterpol = ei.readline()

    rest = ei.readlines()
    
    
    ei.close()
    print(xmin)
    return  

getdata("sciproject")