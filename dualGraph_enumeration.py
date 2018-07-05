#!/usr/local/bin/python
# Swati Jain (S.J.) - script for enumerating dual graph using the rules set in various papers

import sys
from copy import deepcopy
import numpy.linalg as LA
from decimal import *
import os
from numpy import *
from ClassesFunctions import *

# function to read the eigen values and adjacency matrices for dual graphs with one less vertex - 05/18/2018
def readOnelessDualGraphs():

    prevAdjFile = "V%dAdjDG_new"%(numVertices-1) # output files
    prevEigenFile = "%dEigen_new"%(numVertices-1)

    #reading the eigen values, this will create the dual graph class instances and initialize the graph ID, number of vertices, and eigen values
    loadEigenvalues(OnelessGraphs,numVertices-1,prevEigenFile)
    loadAdjMatrices(OnelessGraphs,numVertices-1,prevAdjFile) #reading the adjacency matrices

    #print the dual graphs for testing
    #testFile = "Test_eigen"
    #testFile2 = "Test_Adj"
    #for i in range(0,len(OnelessGraphs)):
    #    printMat(OnelessGraphs[i].adjMatrix,testFile2)
    #    OnelessGraphs[i].printEigen(testFile)


# function to generate all allowed edge combinations by which the new vertex can be connected to onelessgraphs read previously
def genEdgesCombo(edgesLeft,startVertex,prevEdges):

    # recursion termination condition 1
    if edgesLeft == 0: # if not edges are left to be added, add 0's for remaning vertices and add this combo to edgesCombo list
        edgesCombo.append(prevEdges)
        return

    # recursion termination condition 2
    if startVertex == numVertices-1: # this combination not possible with this number of vertices
        return

    for j in range(0,edgesLeft+1):
        if edgesLeft-j != 4: #because there can't be 4 edges between two vertices
            curEdges = deepcopy(prevEdges)
            curEdges[startVertex]=edgesLeft-j
            genEdgesCombo(j,startVertex+1,curEdges)


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


# function to determine if the laplacian of the dual graph has already been enumerated or not - 05/16/2018
def determineUnique(adjMatrix):
    
    # calculate the laplacian and eigen values
    decimalArray = calcEigenValues(adjMatrix)

    # determine if this graph in already there or not
    id = "NA"
    for g in range(0,len(UniqueGraphs)):
        id = UniqueGraphs[g].match(decimalArray,adjMatrix) #06/11/2018 - adding the second argument of adjMatrix to determine unique graphs
        if id != "NA":
            print "Graph already found %d"%id
            return
    if id == "NA":
        print "New graph being added with id %d"%(len(UniqueGraphs)+1)
        UniqueGraphs.append(DualGraph(len(adjMatrix),adjMatrix,len(UniqueGraphs)+1,decimalArray))
        printMat(adjMatrix,adjMatFile)
        UniqueGraphs[len(UniqueGraphs)-1].printEigen(eigenFile)


# recursive function to enumerate dual graphs. Each new vertex is connected to each previous vertex one at a time with 1,2, and 3 edges (again one at a time)
# This misses graphs where the new vertex can be connected to more than one previous vertices - 05/16/2018
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


# function to enumerate graphs with numVertices by adding a new vertex to all onelessgraphs by using edge combinations in edgesCombo
# for each dual graph starting point, this checks the compatibility of each edge combination, and generates unique graphs that follows the rules
def enumerateAll(): # TODO

    for g in OnelessGraphs:
        print "Adding new vertex to starting graph %s"%g.graphID
        
        initialAdjMat = [] # adjacency matrix of the new graph but all new entries are 0
        edgeConn = []
        for i in range(0,numVertices-1):
            deg=0
            tempArray=[]
            for j in range(0,numVertices-1):
                if i != j: # not counting self loops
                    tempArray.append(g.adjMatrix[i][j])
                    deg+=g.adjMatrix[i][j]
                else:
                    tempArray.append(0)
            tempArray.append(0) # for the extra vertex
            initialAdjMat.append(tempArray) # adding the extra column for the new vertex to the Adj matrix
            edgeConn.append(deg) # sroting the connections of the previous vertices
        tempArray = []
        for i in range(0,numVertices): # adding the extra row for the new vertex to the Adj matrix
            tempArray.append(0)
        initialAdjMat.append(tempArray)

        # for each edge combination previously generated
        for edges in edgesCombo:
            print "Checking if this edge combination is compatible with this starting graph"
            compatible = True
            for j in range(0,numVertices-1):
                print "Vertex %d: Edges Present: %d Edges to be added: %d"%(j,edgeConn[j],edges[j])
                if (edgeConn[j] + edges[j]) > 4: # if total edges will become more than 4 then this edge combination is not compatible with this graph
                    compatible = False
                    print "This edge combination is not compatible"
                    break

            if compatible:
                print "This edge combination is compatible. Adding vertex and creating adjacency matrix"
                curAdjMatrix = deepcopy(initialAdjMat)
                for j in range(0,numVertices-1): # adding connections for the new vertex to the adjacency matrix
                    curAdjMatrix[j][numVertices-1] = edges[j]
                    curAdjMatrix[numVertices-1][j] = edges[j]

                # check if the new adjacency matrix follows the rules and if unique, add it to the set of new unique graphs
                check = False
                check = followRules(curAdjMatrix) # add self loops and see if this matrix follows the rules of dual graphs or not
                if check:
                    determineUnique(curAdjMatrix)


# The main function
numVertices = int(sys.argv[1])

UniqueGraphs = [] # New unique dual graphs that will be enumerated
OnelessGraphs = [] # Dual graphs with one vertex less that will be used as starting points
edgesCombo = [] # all allowed edges combination by which the new vertex can be added to onelessgraphs

adjMatFile = "V%dAdjDG_new"%numVertices # output files for newly enumerated graphs
eigenFile = "%dEigen_new"%numVertices

# following lines were used to generate the graphs with 2 vertices (as there is only one previous vertex to connect to)
#adjMatrix_init = []
#for i in range(0,numVertices): # initializing the adjacency matrix
#    tempArray = []
#    for j in range(0,numVertices):
#        tempArray.append(0.0000)
#    adjMatrix_init.append(tempArray)
#enumerate(adjMatrix_init,1)

print "Reading starting graphs with %d vertices"%(numVertices-1)
readOnelessDualGraphs() # read the eigen values and adjacency matrices for dual graphs with one less vertex

print "Generating all allowed edge combinations by which the new vertex with be added to starting graphs read before"
curEdgeCombo=[]
for i in range(0,numVertices-1):
    curEdgeCombo.append(0)
genEdgesCombo(1,0,curEdgeCombo) # generate all allowed edge combinations for 1 egdes
for i in range(0,numVertices-1):
    curEdgeCombo[i] = 0
genEdgesCombo(2,0,curEdgeCombo) # generate all allowed edge combinations for 2 egdes
for i in range(0,numVertices-1):
    curEdgeCombo[i] = 0
genEdgesCombo(3,0,curEdgeCombo) # generate all allowed edge combinations for 3 egdes
for i in range(0,numVertices-1):
    curEdgeCombo[i] = 0
genEdgesCombo(4,0,curEdgeCombo) # generate all allowed edge combinations for 4 egdes
print "Total number of possible edge combinations generated: %d"%len(edgesCombo)
#print the edge combinations for testing
#for i in range(0,len(edgesCombo)):
#    for j in range(0,numVertices-1):
#        print edgesCombo[i][j],
#    print

# Enumerate all graphs
enumerateAll()
