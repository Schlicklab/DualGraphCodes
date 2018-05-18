#!/usr/loca/bin/python
# Swati Jain (S.J.) - python file containig generic classes and functions used by different scripts

import sys
from copy import deepcopy
import numpy.linalg as LA
from decimal import *
import os
from numpy import *

#class to contain information about dual graphs - 05/17/2018
class DualGraph:
    
    def __init__(self,a,b,c,d):
        self.vertices = a
        self.adjMatrix = b
        self.graphID = c
        self.eigenvalues = d
    
    # function to see if eigen values are the same or not
    def match(self,d):
        if self.eigenvalues == d:
            return self.graphID
        else:
            return "NA"

    # to print eigen values of a dual graph - 05/18/2017
    def printEigen(self,filename):
        file = open(filename,'a')
        file.write(">%d"%self.graphID)
        file.write("\n")
        file.close()
        printEigenValues(self.eigenvalues,filename)

# function to print matrices - 05/17/2018
def printMat(Matrix,filename):
    
    file = open(filename,'a')
    file.write("\n")
    for i in range(0,len(Matrix)):
        for j in range(0,len(Matrix)):
            file.write("\t%d"%Matrix[i][j]),
        file.write("\n")
    file.write("\n")
    file.close()

#function to calculate eigen values of a given adjacency matrix - 05/17/2018
def calcEigenValues(adjMatrix):

    laplacian = []
    numVertices = len(adjMatrix)
    for i in range(0,numVertices):
        degree=0
        tempArray = []
        for j in range(0,numVertices):
            if i != j: # ignore self loops
                tempArray.append(-adjMatrix[i][j])
                degree += adjMatrix[i][j]
            else:
                tempArray.append(0)
        tempArray[i] += degree
        laplacian.append(tempArray)

    eigen = sort(LA.eigvals(laplacian))
    decimalArray = []
    decimalPlace = Decimal("0.00000001")
    for i in eigen:
        if isinstance(i, complex): # 04/20/2018 - S.J.
            #print "I have a complex eigen value"
            decimalArray.append(Decimal(str(i.real)).quantize(decimalPlace))
        else:
            decimalArray.append(Decimal(str(i)).quantize(decimalPlace))
    return decimalArray


#function to print eigen values in 0.8 precision that are in the form of a decimalArray - 05/17/2018
def printEigenValues(eigen,filename):

    file = open(filename,'a')
    decimalPlace = Decimal("0.00000001")
    for i in range(0,len(eigen)): # print eigen values upto .8 precision without negative signs
        if str(eigen[i])[0] == "-":
            eigen[i] = Decimal(str(eigen[i])[1:]).quantize(decimalPlace)
        file.write('{0:.8f}'.format(eigen[i]))
        file.write("\n")
    file.close()
