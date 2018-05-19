#!/usr/loca/bin/python
# Swati Jain (S.J.) - python file containig generic classes and functions used by different scripts

import sys
from copy import deepcopy
import numpy.linalg as LA
from decimal import *
import os
from numpy import *

#class to contain information about dual graphs - 05/17/2018
# # copied and modified from 2-find_graphID.py
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
        file.write(">%s"%self.graphID)
        file.write("\n")
        file.close()
        printEigenValues(self.eigenvalues,filename)

    # set functions - 05/18/2018
    def setVertices(self,v):
        self.vertices = v

    def setadjMatrix(self,a):
        self.adjMatrix = a

    def setGraphID(self,g):
        self.graphID = g

    def setEigen(self,e):
        self.eigenvalues = e


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
#copied and modified from dualGraphs.py
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
# copied and modified from calcEigenvals.py
def printEigenValues(eigen,filename):

    file = open(filename,'a')
    decimalPlace = Decimal("0.00000001")
    for i in range(0,len(eigen)): # print eigen values upto .8 precision without negative signs
        if str(eigen[i])[0] == "-":
            eigen[i] = Decimal(str(eigen[i])[1:]).quantize(decimalPlace)
        file.write('{0:.8f}'.format(eigen[i]))
        file.write("\n")
    file.close()


# function to load all dual graph IDs and eigen values from a given file and create dual graph instances - 05/18/2018
# copied and modified from dualGraphs.py
# creates the dual graph instances (with vertices, graphID, and eigenvalues set) and adds them to the Graphs list.
def loadEigenvalues(Graphs,num_vertices,file):
    
    decimalPlace = Decimal("0.00000001")
    emptyArray = []
    tArray = []
    graph_num = 0
    f = open(file,'r')
    for line in f:
        if(len(line) > 0 and line[0] == '>'): # reading in the graph id
            key = line[1:-1]
        elif(len(line) > 0 and line[0] != '>'): # reaing in the eigen values
            tArray.append(Decimal(str(line[:-1])).quantize(decimalPlace))
        if len(tArray) == num_vertices: # when enough eigen values are read
            Graphs.append(DualGraph(num_vertices,emptyArray,key,tArray)) # creating a dual graph instance and adding it to the graph list
            graph_num+=1
            tArray = []
    f.close()
    print "Read eigan values for %d dual graphs from file %s"%(graph_num,file)


# function to load all adj matrices from given files into the Graphs list - 05/18/2018
# copied and modified from calcEigenvals.py
# needs the Graphs list to be already initialized with dual graph instances
# Assuming that the adjacency matrices are in the same order as eigen values file that was used to initialize the dual graphs in Graphs
def loadAdjMatrices(Graphs,num_vertices,file):

    tempAdjMatrix = []
    tempArray = []
    graph_num=0
    size=0
    f = open(file,'r')
    for line in f: # read the adjacency matrices
        if line != '\n':
            size +=  1
            tempArray = []
            for x in line.split():
                tempArray.append(int(x))
            tempAdjMatrix.append(tempArray)
            if size == num_vertices: # when enough rows are read, then add the adjacency matrix to the Graphs list
                Graphs[graph_num].setadjMatrix(tempAdjMatrix)
                graph_num += 1
                #re-initialize the variables
                tempAdjMatrix = []
                size = 0
    f.close()
    print "Read adjacency matrices for %d dual graphs from file %s"%(graph_num,file)


