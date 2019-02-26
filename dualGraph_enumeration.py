#!/usr/local/bin/python
# Swati Jain (S.J.) - script for enumerating dual graph using the rules set in various papers
# 07/06/2018 - takes three integers as command line arguments, 1: number of vertices for which graphs have to be enumerated, 2: number of vertices of subgraph 1, 3: number of vertices of subgraph 2. All graphs with sub 1 vertices and all graphs with sub 2 vertices will be combined to generate new graphs. Assumption: sub 1 vertex >= sub 2 vertex, therefore sub 2 graphs will be added to sub 1 graph

import sys
from copy import deepcopy
import numpy.linalg as LA
from decimal import *
import os
import time
from numpy import *
from ClassesFunctions import *

# function to read the eigen values and adjacency matrices for dual graphs - 05/18/2018
# changes to make it read any vertex number dual graphs that will be used as subgraphs - 07/06/2018
def readSubDualGraphs(vertex,Graphs):

    #print "Reading starting graphs with %d vertices"%(vertex)
    
    prevAdjFile = "V%dAdjDG"%(vertex) # output files
    prevEigenFile = "%dEigen"%(vertex)

    #reading the eigen values, this will create the dual graph class instances and initialize the graph ID, number of vertices, and eigen values
    loadEigenvalues(Graphs,vertex,prevEigenFile)
    loadAdjMatrices(Graphs,vertex,prevAdjFile) #reading the adjacency matrices

    #print the dual graphs for testing
    #testFile = "Test_eigen"
    #testFile2 = "Test_Adj"
    #for i in range(0,len(Graphs)):
    #    printMat(Graphs[i].adjMatrix,testFile2)
    #    Graphs[i].printEigen(testFile)


# function to generate all allowed edge combinations by which the a new vertex can be connected to sub 1 graphs read previously
# changes to work with any number of vertices in sub 1 - 07/06/2018
def genEdgesCombo(edgesLeft,startVertex,prevEdges):

    # recursion termination condition 1
    if edgesLeft == 0: # if not edges are left to be added, add 0's for remaning vertices and add this combo to edgesCombo list
        edgesCombo.append(prevEdges)
        return

    # recursion termination condition 2
    if startVertex == numVertices_sub1: # this combination not possible with this number of vertices
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
        else:
            return False # 09/17/2018 - if any other degree, then that does not follow the rules
        degree.append(edges)

    if numThree != 0 and numThree == 2: #if there are vertices with degree 3 present, there can only be 2 vertices with degree 3
        #print "This graph follows the rules"
        for i in range(0,numV):
            if degree[i] == 2: # if any vertex has only 2, then make it 4 by adding self loop
                adjMatrix[i][i] = 2
                degree[i]+=2
        return True
    elif numThree == 0 and numTwo != 0: # if there are no vertices with degree 3, then there has to be atleast one vertex with degree 2
        #print "This graph follows the rules"
        firstTwo = False
        for i in range(0,numV):
            if degree[i] == 2 and not firstTwo: # leave the first occurance of degree 2 vertex as is (as we need one and only one vertex of degree 2)
                firstTwo = True
            elif degree[i] == 2 and firstTwo: # changes subsequent vertices of degree 2 to degree 4 by adding self loops
                adjMatrix[i][i] = 2
                degree[i]+=2
        return True
    else:
        #print "This graph does not follow the rules"
        return False


# function to determine if the laplacian of the dual graph has already been enumerated or not - 05/16/2018
def determineUnique(adjMatrix):
    
    #global dict_time
    
    global total_new_graphs
    # calculate the laplacian and eigen values
    decimalArray = calcEigenValues(adjMatrix)
    eigen_tuple=tuple(decimalArray) # 09/14/2018 - implementation of the search using a python dictionary

    # 09/14/2018 - implementation of the search using a python dictionary
    Graphs=UniqueGraphs.get(eigen_tuple,[])
    flag=False
    #start_time_5=time.clock()
    for g in Graphs:
	    flag=checkIsomorphism(g.adjMatrix,adjMatrix)
            if flag: # if this is already present, just return
                #return
                break

    #end_time_5=time.clock()
    #dict_time+=end_time_5-start_time_5

    # code will come here if the current graph is not present, so needs to be added
    if not flag:
        total_new_graphs+=1
        new_g=DualGraph(len(adjMatrix),adjMatrix,total_new_graphs,decimalArray)
        new_g.printEigen(eigenFile)
        printMat(new_g.adjMatrix,adjMatFile)
        Graphs.append(new_g)
        UniqueGraphs[eigen_tuple]=Graphs


# function to enumerate graphs with numVertices by adding a new vertex/or graph with sub 2 number of vertices to all graphs with sub 1 number of vertices by using edge combinations in edgesCombo - 07/06/2018
# for each dual graph starting point, this checks the compatibility of each edge combination, and generates unique graphs that follows the rules
# 07/06/2018 - the three arguments are: the adjacency matrix (without self loops) of the sub 2 graph that will be added to all sub 1 graphs, list of total number of edges to be added to sub 1 graph's vertices, list of indices in the edgeCombo array for edges between each vertex of sub 2 and sub 1 graphs (each corresponding to the same entry in the list of total number of edges)
def enumerateAllSub1(toAddSubMat,numEdgestoAdd,edgesComboIndices):

    for g in Sub1_Graphs:

	#if g.graphID != "13137":
	#	continue
        
        print "Adding to starting graph %s"%g.graphID
        
        initialAdjMat = [] # adjacency matrix of the new graph but all new entries are 0
        edgeConn = []
        for i in range(0,numVertices_sub1):
            deg=0
            tempArray=[]
            for j in range(0,numVertices_sub1):
                if i != j: # not counting self loops
                    tempArray.append(g.adjMatrix[i][j])
                    deg+=g.adjMatrix[i][j]
                else:
                    tempArray.append(0)
            edgeConn.append(deg) # storing the connections of the previous vertices
            for j in range(0,numVertices_sub2):
                tempArray.append(0) # for the extra vertices
            initialAdjMat.append(tempArray) # adding the extra column for the new vertex to the Adj matrix
    
        for j in range(0,numVertices_sub2): # adding the extra rows for the new vertices to the Adj matrix
            tempArray = []
            for i in range(0,numVertices_sub1):
                tempArray.append(0)
            for k in range(0,numVertices_sub2):
                tempArray.append(toAddSubMat[j][k]); # copying the adj matrix of the sub 2 graph to be added (provided as an argument)
            initialAdjMat.append(tempArray)

        # for each edge combination previously generated
        for e in range(0,len(numEdgestoAdd)):
            #print "Checking if this edge combination is compatible with this starting graph"
            compatible = True
            for j in range(0,numVertices_sub1):
                #print "Vertex %d: Edges Present: %d Edges to be added: %d"%(j,edgeConn[j],numEdgestoAdd[e][j])
                if (edgeConn[j] + numEdgestoAdd[e][j]) > 4: # if total edges for any one vertex will become more than 4 then this edge combination is not compatible with this graph
                    compatible = False
                    #print "This edge combination is not compatible"
                    break

            if compatible:
                #print "This edge combination is compatible. Adding vertex and creating adjacency matrix"
                curAdjMatrix = deepcopy(initialAdjMat)
                for k in range(0,numVertices_sub2):
                    for j in range(0,numVertices_sub1): # adding connections for the new vertices to the adjacency matrix
                        curAdjMatrix[numVertices_sub1+k][j] = edgesCombo[edgesComboIndices[e][k]][j];
                        curAdjMatrix[j][numVertices_sub1+k] = edgesCombo[edgesComboIndices[e][k]][j];

                # check if the new adjacency matrix follows the rules and if unique, add it to the set of new unique graphs
                check = False
                check = followRules(curAdjMatrix) # add self loops and see if this matrix follows the rules of dual graphs or not
                if check:
                    determineUnique(curAdjMatrix)

# 09/20 function to generate and check graphs that are only compatible with the second graph
def generateGraphs_0920(V2,ListIndices,edges_left_sub1,initialAdjMat):

    global first
    #global rules_time
    #global determine_time
    global total_new_graphs
    
    if V2 == numVertices_sub2: # all vertices are checked
        
        if first: # to avoid getting un conencted graphs
            first = False
            return

        check = False
        #start_time_3=time.clock()
        check = followRules(initialAdjMat)
        #end_time_3=time.clock()
        #rules_time+=end_time_3-start_time_3
        if check:
            #start_time_4=time.clock()
            determineUnique(initialAdjMat)
            
            #decimalArray = calcEigenValues(initialAdjMat)
            #eigen_tuple=tuple(decimalArray) # 09/14/2018 - implementation of the search using a python dictionary
            #Graphs=UniqueGraphs.get(eigen_tuple,[])
            #flag=False
            #for g in Graphs:
            #    flag=checkIsomorphism(g.adjMatrix,initialAdjMat)
            #    if flag: # if this is already present, just return
            #        return
    
            # code will come here if the current graph is not present, so needs to be added
            #total_new_graphs+=1
            #new_g=DualGraph(len(initialAdjMat),initialAdjMat,total_new_graphs,decimalArray)
            #new_g.printEigen(eigenFile)
            #printMat(new_g.adjMatrix,adjMatFile)
            #Graphs.append(new_g)
            #UniqueGraphs[eigen_tuple]=Graphs
            
            #end_time_4=time.clock()
            #print end_time_4-start_time_4
            #determine_time+=end_time_4-start_time_4
        return

    for e in range(0,len(ListIndices[V2])): # for each edge combinationa compatible with vertex V2 of sub 2 graph

        compatible = True
        edge = ListIndices[V2][e]
        for V1 in range(0,numVertices_sub1):

            if edgesCombo[edge][V1] > edges_left_sub1[V1]: # to see if this edge combinationa agrees with this vertex of the sub 1 graph
                compatible = False # not compatible
                break

        if compatible:
            curAdjMatrix = deepcopy(initialAdjMat)
            cur_edges_left_sub1=deepcopy(edges_left_sub1)
            for V1 in range(0,numVertices_sub1):
                curAdjMatrix[numVertices_sub1+V2][V1]=edgesCombo[edge][V1]
                curAdjMatrix[V1][numVertices_sub1+V2]=edgesCombo[edge][V1]
                cur_edges_left_sub1[V1]-=edgesCombo[edge][V1]

            generateGraphs_0920(V2+1,ListIndices,cur_edges_left_sub1,curAdjMatrix) # call the recursive function for the next vertex


# 09/17/2018 - generate graphs for each pair of g1 and g2, so that only edge combinations that are compatible are generated
def enumerateAll_0917():

    #global Edges_left_sub1
    #global Edges_left_sub2
    global first
    #global genGraphs_time
    
    basicAdjMat = [[0] * (numVertices_sub1+numVertices_sub2) for a in range(numVertices_sub1+numVertices_sub2)] # the starting matrix with every element 0
    Edges_left_sub1 = [4]*numVertices_sub1
    Edges_left_sub2 = [4]*numVertices_sub2

    for g2 in Sub2_Graphs: # for each sub 2 graph
        
        #if g2.graphID == "3_1" or g2.graphID == "3_2" or g2.graphID == "3_3" or g2.graphID == "3_4":
            continue
        
        #if g2.graphID == "4_1" or g2.graphID == "4_2" or g2.graphID == "4_3" or g2.graphID == "4_4" or g2.graphID == "4_5" or g2.graphID == "4_6" or g2.graphID == "4_7" or g2.graphID == "4_8" or g2.graphID == "4_9" or g2.graphID == "4_10" or g2.graphID == "4_11" or g2.graphID == "4_12" or g2.graphID == "4_13" or g2.graphID == "4_14" or g2.graphID == "4_16" or g2.graphID == "4_17" or g2.graphID == "4_18":
            continue

        #print "Sub 2 graph id: %s"%g2.graphID

        basicAdjMat[0:numVertices_sub2+numVertices_sub1] = [[0] * (numVertices_sub1+numVertices_sub2) for a in range(numVertices_sub1+numVertices_sub2)]
        Edges_left_sub2[0:numVertices_sub2] = [4]*numVertices_sub2
        for i in range(0,numVertices_sub2): # copying the sub 2 graph adjMatrix and number of edges left
            for j in range(0,numVertices_sub2):
                if i != j:
                    basicAdjMat[i+numVertices_sub1][j+numVertices_sub1] = g2.adjMatrix[i][j]
                    Edges_left_sub2[i]-=g2.adjMatrix[i][j]
        
        
        # 09/20 - Adding the listIndices part of the enumerateAll() function here, but slightly changed
        ListIndices=[]
        for i in range(0,numVertices_sub2):
            
            tempIndex=[]
            for e in range(0,len(edgesCombo)):
                
                deg = 0
                for j in range(0,numVertices_sub1):
                    deg+=edgesCombo[e][j]
                
                if deg <= Edges_left_sub2[i]:
                    tempIndex.append(e)
                else:
                    break
            ListIndices.append(tempIndex)

        #for i in range(0,numVertices_sub2):
        #    for j in range(0,len(ListIndices[i])):
        #        print ListIndices[i][j],
        #    print


        for g1 in Sub1_Graphs: # for each sub 1 graph
            
            #if g1.graphID != "3_2":
            #    continue
            
            initialAdjMat = deepcopy(basicAdjMat)
            first=True
            Edges_left_sub1[0:numVertices_sub1] = [4]*numVertices_sub1
            for i in range(0,numVertices_sub1): # copying the sub 1 graph adjMatrix and number of edges left
                for j in range(0,numVertices_sub1):
                    if i != j:
                        initialAdjMat[i][j] = g1.adjMatrix[i][j]
                        Edges_left_sub1[i]-=g1.adjMatrix[i][j]
            
            #start_time_2=time.clock()
            generateGraphs_0920(0,ListIndices,Edges_left_sub1,initialAdjMat) # add sub 2 graph g2 to the sub 1 graph g1
            #end_time_2=time.clock()
            #genGraphs_time+=end_time_2-start_time_2



# The main function
#start_time = time.clock()

numVertices = int(sys.argv[1])
numVertices_sub1 = int(sys.argv[2]) # 07/06/2018 - the number of vertices of two subgraphs that will be combined
numVertices_sub2 = int(sys.argv[3])

UniqueGraphs = dict() # 09/14/2018 - implementation of the search using a python dictionary
total_new_graphs=0

# 07/06/2018 - changes to work with sub 2 graphs with any number of vertices
Sub1_Graphs = [] # Dual graphs with vertex number sub1
Sub2_Graphs = [] # Dual graphs with vertex number sub2
edgesCombo = [] # all allowed edges combination by which the one of the new vertices of sub 2 graphs can be added to sub 1 graphs
first=True
#genGraphs_time=0
#rules_time=0
#determine_time=0
#dict_time=0
#Edges_left_sub1 = [4]*numVertices_sub1
#Edges_left_sub2 = [4]*numVertices_sub2

adjMatFile = "V%dAdjDG_%d_%d"%(numVertices,numVertices_sub1,numVertices_sub2) # output files for newly enumerated graphs
eigenFile = "%dEigen_%d_%d"%(numVertices,numVertices_sub1,numVertices_sub2)

# 07/06/2018 - to read both subgraphs that will be used
readSubDualGraphs(numVertices_sub1,Sub1_Graphs) # read the eigen values and adjacency matrices for dual graphs for Sub1
readSubDualGraphs(numVertices_sub2,Sub2_Graphs) # read the eigen values and adjacency matrices for dual graphs for Sub2

#start_time_1=time.clock()
# 07/06/2018 - changes to work with adding not just one vertex, but graphs of sub 2 vertex numbers to graphs of sub 1 vertex numbers
#print "Generating all allowed edge combinations by which the new graphs with be added to starting graphs read before"
curEdgeCombo=[]
for i in range(0,numVertices_sub1):
     curEdgeCombo.append(0)
edgesCombo.append(curEdgeCombo) # generate all allowed edge combinations for 0 egdes
for i in range(0,numVertices_sub1):
     curEdgeCombo[i] = 0
genEdgesCombo(1,0,curEdgeCombo) # generate all allowed edge combinations for 1 egdes
for i in range(0,numVertices_sub1):
    curEdgeCombo[i] = 0
genEdgesCombo(2,0,curEdgeCombo) # generate all allowed edge combinations for 2 egdes
for i in range(0,numVertices_sub1):
    curEdgeCombo[i] = 0
genEdgesCombo(3,0,curEdgeCombo) # generate all allowed edge combinations for 3 egdes
for i in range(0,numVertices_sub1):
    curEdgeCombo[i] = 0
genEdgesCombo(4,0,curEdgeCombo) # generate all allowed edge combinations for 4 egdes
print "Total number of possible edge combinations generated for one vertex: %d"%len(edgesCombo)
#print the edge combinations for testing
#for i in range(0,len(edgesCombo)):
#    for j in range(0,numVertices_sub1):
#        print edgesCombo[i][j],
#    print
#end_time_1=time.clock()
#print "The total time taken in seconds to generate edges for one single vertex:",
#print (end_time_1 - start_time_1)

# Enumerate all graphs
enumerateAll_0917()

#end_time = time.clock()

#print "The total time taken in seconds to generate graphs (not the checking part for sub 2 graph): ",
#print genGraphs_time

#print "Time spent in follow rules function: ",
#print rules_time

#print "Time spent in determine unique function: ",
#print determine_time

#print "Time spent in checking for isomorphism: ",
#print dict_time

#print "The total time taken in seconds is:",
#print (end_time - start_time)
