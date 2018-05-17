#!/usr/local/bin/python
# Swati Jain (S.J.) - script for enumerating dual graph using the rules set in various papers

import sys
from copy import deepcopy
import numpy.linalg as LA
from decimal import *
import os
from numpy import *

#class to contain information about new dual graphs
class DualGraph:
    
    def __init__(self,a,b):
        self.graphID = a
        self.values = b
        
    def match(self,b):
        if self.values == b:
            return self.graphID
        else:
            return "NA"

# function to print adj matrices
def printAdjMat(adjMatrix):
    
    print "Adjacency matrix"
    for i in range(0,numVertices):
        for j in range(0,numVertices):
            print adjMatrix[i][j],
        print

# function to determine if the laplacian of the dual graph has already been enumerated or not
def determineUnique(adjMatrix):
    
    # calculate the laplacian and eigen values
    laplacian = []
    for i in range(0,numVertices):
        degree=0
        tempArray = []
        for j in range(0,numVertices):
            tempArray.append(-adjMatrix[i][j])
            degree += adjMatrix[i][j]
        tempArray[i] += degree
        laplacian.append(tempArray)
    
    eigen = sort(LA.eigvals(laplacian))
    decimalArray = []
    decimalPlace = Decimal("0.00000001")
    for i in eigen:
        if isinstance(i, complex): # 04/20/2018 - S.J.
            print "I have a complex eigen value"
            decimalArray.append(Decimal(str(i.real)).quantize(decimalPlace))
        else:
            decimalArray.append(Decimal(str(i)).quantize(decimalPlace))

    # determine if this graph in already there or not
    id = "NA"
    for g in range(0,len(UniqueGraphs)):
        id = UniqueGraphs[g].match(decimalArray)
        if id != "NA":
            print "Graph already found %d"%id
            return
    if id == "NA":
        print "New graph being added with id %d"%(len(UniqueGraphs)+1)
        UniqueGraphs.append(DualGraph(len(UniqueGraphs)+1,decimalArray))
        printAdjMat(adjMatrix)
        print "Eigen values:"
        for i in range(0,len(decimalArray)): # print eigen values upto .8 precision without negative signs
            if str(decimalArray[i])[0] == "-":
                decimalArray[i] = Decimal(str(decimalArray[i])[1:]).quantize(decimalPlace)
            print '{0:.8f}'.format(decimalArray[i])

# recursive function to enumerate dual graphs. Each new vertex is connected to each previous vertex one at a time with 1,2, and 3 edges (again one at a time)
# Right now this misses graphs where the new vertex can be connected to more than one previous vertices
def enumerate(adjMatrix,curVertex):
    
    if(curVertex == numVertices): #recursion ending condition
        print "Full graph connected"
        #printAdjMat(adjMatrix)
        determineUnique(adjMatrix)
        return

    print "Adding vertex %d to the graph"%curVertex
    for prev in range(0,curVertex): #range over all previous vertices
        print "Connecting it to Vertex %d"%prev
        for edges in range(1,4): #range over all possible number of edges between two vertices
            print "With %d edges"%edges
            curAdjMatrix = deepcopy(adjMatrix) #using the original adjacency matrix
            curAdjMatrix[prev][curVertex] = edges
            curAdjMatrix[curVertex][prev] = edges
            degree_cur=0
            degree_prev=0
            for i in range(0,numVertices):
                degree_prev += curAdjMatrix[prev][i]
                degree_cur += curAdjMatrix[curVertex][i]
            if degree_cur > 4 or degree_prev > 4:
                print "Graph not accoring to rules. Not going forward"
            else:
                print "Going forward with enumeration"
                enumerate(curAdjMatrix,curVertex+1)


numVertices = int(sys.argv[1])
UniqueGraphs = []
adjMatrix_init = []

for i in range(0,numVertices): # initializing the adjacency matrix
    tempArray = []
    for j in range(0,numVertices):
        tempArray.append(0.0000)
    adjMatrix_init.append(tempArray)

enumerate(adjMatrix_init,1)
