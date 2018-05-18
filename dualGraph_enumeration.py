#!/usr/local/bin/python
# Swati Jain (S.J.) - script for enumerating dual graph using the rules set in various papers

import sys
from copy import deepcopy
import numpy.linalg as LA
from decimal import *
import os
from numpy import *
from ClassesFunctions import *

# function to determine if a adjMatrix follows the rules for dual graphs or not - 05/17/2018
def followRules(adjMatrix):

    # to keep track of number of vertices with degree two, three, or four
    degree = []
    numTwo = 0
    numThree = 0
    numFour = 0
    numV = len(adjMatrix)
    
    #counting degree of each vertex
    for i in range(0,numV):
        edges = 0
        for j in range(0,numV):
            edges+=adjMatrix[i][j]
        if edges == 1: # if there is only one edge, then a self loop has to be added
            adjMatrix[i][i] = 2
            edges +=2
            numThree+=1
        elif edges == 2:
            numTwo+=1
        elif edges == 3:
            numThree+=1
        elif edges == 4:
            numFour+=1
        degree.append(edges)

    if numThree != 0 and numThree == 2: #if there are vertices with degree 3 present, there can only be 2 vertices with degree 3
        print "This graph follows the rules"
        for i in range(0,numV):
            if degree[i] == 2: # if any vertex has only 2, then make it 4 by adding self loop
                adjMatrix[i][i] = 2
                degree[i]+=2
        return True
    elif numThree == 0 and numTwo != 0: # if there are no vertices with degree 3, then there has to be atleast one vertex with degree 2
        print "This graph follows the rules"
        firstTwo = False
        for i in range(0,numV):
            if degree[i] == 2 and not firstTwo: # leave the first occurance of degree 2 vertex as is (as we need one and only one vertex of degree 2)
                firstTwo = True
            elif degree[i] == 2 and firstTwo: # changes subsequent vertices of degree 2 to degree 4 by adding self loops
                adjMatrix[i][i] = 2
                degree[i]+=2
        return True
    else:
        print "This graph does not follow the rules"
        return False


# function to determine if the laplacian of the dual graph has already been enumerated or not - - 05/16/2018
def determineUnique(adjMatrix):
    
    # calculate the laplacian and eigen values
    decimalArray = calcEigenValues(adjMatrix)

    # determine if this graph in already there or not
    id = "NA"
    for g in range(0,len(UniqueGraphs)):
        id = UniqueGraphs[g].match(decimalArray)
        if id != "NA":
            print "Graph already found %d"%id
            return
    if id == "NA":
        print "New graph being added with id %d"%(len(UniqueGraphs)+1)
        UniqueGraphs.append(DualGraph(len(adjMatrix),adjMatrix,len(UniqueGraphs)+1,decimalArray))
        printMat(adjMatrix,adjMatFile)
        UniqueGraphs[len(UniqueGraphs)-1].printEigen(eigenFile)


# recursive function to enumerate dual graphs. Each new vertex is connected to each previous vertex one at a time with 1,2, and 3 edges (again one at a time)
# Right now this misses graphs where the new vertex can be connected to more than one previous vertices - 05/16/2018
def enumerate(adjMatrix,curVertex):
    
    if(curVertex == numVertices): #recursion ending condition
        print "Full graph connected"
        flag = False
        flag = followRules(adjMatrix) # add self loops and see if this matrix follows the rules of dual graphs or not
        if flag:
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
adjMatFile = "V%dAdjDG_new"%numVertices
eigenFile = "%dEigen_new"%numVertices

for i in range(0,numVertices): # initializing the adjacency matrix
    tempArray = []
    for j in range(0,numVertices):
        tempArray.append(0.0000)
    adjMatrix_init.append(tempArray)

enumerate(adjMatrix_init,1)
