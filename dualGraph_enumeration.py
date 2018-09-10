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

# 09/07/2018 - binary tree node to store dual graphs for better search time
class BinaryTreeNode:

    def __init__(self,values=None):
        self.key = values
        self.graphs = []
        self.leftNode = None
        self.rightNode = None
    
    def addGraph(self,adjMatrix,graphID,adjMatFile,eigenFile):
        g = DualGraph(len(adjMatrix),adjMatrix,graphID,self.key)
        g.printEigen(eigenFile)
        printMat(g.adjMatrix,adjMatFile)
        self.graphs.append(g)
    
    def insert(self,eigenvalues,adjMatrix,total_graphs,adjMatFile,eigenFile):

        if self.key == None: # first graph to be created
            total_graphs=1
            self.key = eigenvalues
            self.addGraph(adjMatrix,total_graphs,adjMatFile,eigenFile)
            return total_graphs

        else: # the tree is not empty

            curNode = self
            while(curNode != None):

                if(eigenvalues < curNode.key): # go to the left subtree as eigenvalues are smaller
                    if curNode.leftNode != None:
                        curNode = curNode.leftNode
                    else:
                        total_graphs+=1
                        curNode.leftNode = BinaryTreeNode(eigenvalues)
                        curNode.leftNode.addGraph(adjMatrix,total_graphs,adjMatFile,eigenFile)
                        return total_graphs
            
                elif(eigenvalues > curNode.key): #go to the right subtree as eigenvalues are greater
                    if curNode.rightNode != None:
                        curNode = curNode.rightNode
                    else:
                        total_graphs+=1
                        curNode.rightNode = BinaryTreeNode(eigenvalues)
                        curNode.rightNode.addGraph(adjMatrix,total_graphs,adjMatFile,eigenFile)
                        return total_graphs

                else: # decimalArray is the same as this node's eigenvalues, so now check for isomorphism
                    searchGraphs = curNode.graphs
                    flag=False

                    for g in searchGraphs: # for each graph in this node, check if the current graph is isomorphic or not
                        flag=checkIsomorphism(g.adjMatrix,adjMatrix)
                        if flag: # if this is already present, just return
                            return total_graphs

                    # code will come here if the current graph is not present, so needs to be added
                    total_graphs+=1
                    curNode.addGraph(adjMatrix,total_graphs,adjMatFile,eigenFile)
                    return total_graphs


# function to read the eigen values and adjacency matrices for dual graphs - 05/18/2018
# changes to make it read any vertex number dual graphs that will be used as subgraphs - 07/06/2018
def readSubDualGraphs(vertex,Graphs):

    #print "Reading starting graphs with %d vertices"%(vertex)
    
    prevAdjFile = "V%dAdjDG_map_sort"%(vertex) # output files
    prevEigenFile = "%dEigen_map_sort"%(vertex)

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
    
    global total_new_graphs
    # calculate the laplacian and eigen values
    decimalArray = calcEigenValues(adjMatrix)

    # determine if this graph in already there or not
    #id = "NA"
    #for g in range(0,len(UniqueGraphs)):
    #    id = UniqueGraphs[g].match(decimalArray,adjMatrix) #06/11/2018 - adding the second argument of adjMatrix to determine unique graphs
    #    if id != "NA":
            #print "Graph already found %d"%id
    #        return
    #if id == "NA":
        #print "New graph being added with id %d"%(len(UniqueGraphs)+1)
    #    UniqueGraphs.append(DualGraph(len(adjMatrix),adjMatrix,len(UniqueGraphs)+1,decimalArray))
    #    printMat(adjMatrix,adjMatFile)
    #    UniqueGraphs[len(UniqueGraphs)-1].printEigen(eigenFile)

    # 09/07/2018 - implementation of a binary tree for insertion and search
    total_new_graphs=UniqueGraphs.insert(decimalArray,adjMatrix,total_new_graphs,adjMatFile,eigenFile)


# recursive function to enumerate dual graphs. Each new vertex is connected to each previous vertex one at a time with 1,2, and 3 edges (again one at a time)
# This misses graphs where the new vertex can be connected to more than one previous vertices - 05/16/2018
def enumerate(adjMatrix,curVertex):
    
    if(curVertex == numVertices): #recursion ending condition
        #print "Full graph connected"
        flag = False
        flag = followRules(adjMatrix) # add self loops and see if this matrix follows the rules of dual graphs or not
        if flag:
            determineUnique(adjMatrix)
        return

    #print "Adding vertex %d to the graph"%curVertex
    for prev in range(0,curVertex): #range over all previous vertices
        #print "Connecting it to Vertex %d"%prev
        for edges in range(1,4): #range over all possible number of edges between two vertices
            #print "With %d edges"%edges
            curAdjMatrix = deepcopy(adjMatrix) #using the original adjacency matrix
            curAdjMatrix[prev][curVertex] = edges
            curAdjMatrix[curVertex][prev] = edges
            degree_cur=0
            degree_prev=0
            for i in range(0,numVertices):
                degree_prev += curAdjMatrix[prev][i]
                degree_cur += curAdjMatrix[curVertex][i]
            if degree_cur > 4 or degree_prev > 4:
                next
                #print "Graph not accoring to rules. Not going forward"
            else:
                #print "Going forward with enumeration"
                enumerate(curAdjMatrix,curVertex+1)


# function to enumerate graphs with numVertices by adding a new vertex/or graph with sub 2 number of vertices to all graphs with sub 1 number of vertices by using edge combinations in edgesCombo - 07/06/2018
# for each dual graph starting point, this checks the compatibility of each edge combination, and generates unique graphs that follows the rules
# 07/06/2018 - the three arguments are: the adjacency matrix (without self loops) of the sub 2 graph that will be added to all sub 1 graphs, list of total number of edges to be added to sub 1 graph's vertices, list of indices in the edgeCombo array for edges between each vertex of sub 2 and sub 1 graphs (each corresponding to the same entry in the list of total number of edges)
def enumerateAllSub1(toAddSubMat,numEdgestoAdd,edgesComboIndices):

    for g in Sub1_Graphs:

	#if g.graphID != "13137":
	#	continue
        
        #print "Adding to starting graph %s"%g.graphID
        
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


# 07/08/2018 - to select edge combinations that are compatible with a given sub 2 graph (from which this function is called)
# 07/30/2018 - adding the last argument for completely symmetric graphs (2_1,2_2,2_3,3_5, and 4_19)
def selectEdgeCombos(V,edgeConn,prevIndex,prevNumEdges,ListIndices,ListNumEdges,index):

    if V == numVertices_sub2: # all vertices have been covered for this edge combinations
        # add the edge combination to the list of edge indices total number of edges to list of total number of edges
        ListIndices.append(prevIndex)
        ListNumEdges.append(prevNumEdges)
        return

    #for e in range(0,len(edgesCombo)):
    for e in range(index,len(edgesCombo)):

        deg = 0
        for i in range(0,numVertices_sub1):
            deg+=edgesCombo[e][i]

        if deg + edgeConn[V] <= 4: # this edge is compatible with this sub 2 vertex V

            curNumEdges = deepcopy(prevNumEdges)
            for i in range(0,numVertices_sub1):
                curNumEdges[i]+=edgesCombo[e][i] # adding the number of edges that will have to be added to sub 1 graph vertices using this edge combo

            curIndex = deepcopy(prevIndex)
            curIndex[V]=e # adding the index for the compatible edge combination for the sub 2 vertex V

            selectEdgeCombos(V+1,edgeConn,curIndex,curNumEdges,ListIndices,ListNumEdges,0)
            #selectEdgeCombos(V+1,edgeConn,curIndex,curNumEdges,ListIndices,ListNumEdges,e) # use this for completely symmetric graphs

        #else:
            #print "Edge combination %d not compatible with vertex %d"%(e+1,V+1)


# 07/06/2018 - function to combine edge combinations for all vertices in sub 2 vertices graphs, and then calling enumerateAllSub1 for each sub 2 vertex graph
def enumerateAll():

    for g in Sub2_Graphs:
        
        #if g.graphID != sys.argv[4]: # to work only for one specific sub2 graph - used for dividing large runs
        #    continue
        #if g.graphID == "1" or g.graphID == "2" or g.graphID == "4" or g.graphID == "8": #specifically for the ones that are not non-separable blocks for 3 vertices
        #    continue
        
        #print "Analzying edge combinations to be added for graph %s"%g.graphID

        basicAdjMat = [] # adjacency matrix for graph g without self loops
        edgeConn = []
        for i in range(0,numVertices_sub2):
            deg=0
            tempArray=[]
            for j in range(0,numVertices_sub2):
                if i != j: # not counting self loops
                    tempArray.append(g.adjMatrix[i][j])
                    deg+=g.adjMatrix[i][j]
                else:
                    tempArray.append(0)
            edgeConn.append(deg) # storing the connections already present in the graph
            basicAdjMat.append(tempArray)

        ListNumEdges=[]
        ListIndices=[]
        tempIndex=[]
        tempNumEdges=[]
        
        for i in range(0,numVertices_sub2):
            tempIndex.append(0)
        for i in range(0,numVertices_sub1):
            tempNumEdges.append(0)

        #selectEdgeCombos(0,edgeConn,tempIndex,tempNumEdges,ListIndices,ListNumEdges) # making a list of edge combinations that are compatible with this sub 2 graph
        selectEdgeCombos(0,edgeConn,tempIndex,tempNumEdges,ListIndices,ListNumEdges,0)
        
        del ListIndices[0] # because the first is always not adding any edges to by vertices, 0 always
        del ListNumEdges[0]
	
        #print "Number of edge combinations for all vertices compatible with this graph: %d"%(len(ListNumEdges))
        #for i in range(0,len(ListIndices)):
        #    for j in range(0,numVertices_sub2):
        #        print ListIndices[i][j],
        #    print

        enumerateAllSub1(basicAdjMat,ListNumEdges,ListIndices) # generate new graphs by adding this sub 2 graph to all sub 1 graphs using the edge combinations selected above


# The main function
start_time = time.clock()

numVertices = int(sys.argv[1])
numVertices_sub1 = int(sys.argv[2]) # 07/06/2018 - the number of vertices of two subgraphs that will be combined
numVertices_sub2 = int(sys.argv[3])



#UniqueGraphs = [] # New unique dual graphs that will be enumerated
UniqueGraphs = BinaryTreeNode() # 09/07/2018 - binary search tree
total_new_graphs=0

# 07/06/2018 - changes to work with sub 2 graphs with any number of vertices
Sub1_Graphs = [] # Dual graphs with vertex number sub1
Sub2_Graphs = [] # Dual graphs with vertex number sub2
edgesCombo = [] # all allowed edges combination by which the one of the new vertices of sub 2 graphs can be added to sub 1 graphs

number_file = int(sys.argv[4])
adjMatFile = "time_test/V%dAdjDG_%d_%d_%d"%(numVertices,numVertices_sub1,numVertices_sub2,number_file) # output files for newly enumerated graphs
eigenFile = "time_test/%dEigen_%d_%d_%d"%(numVertices,numVertices_sub1,numVertices_sub2,number_file)

# following lines were used to generate the graphs with 2 vertices (as there is only one previous vertex to connect to)
#adjMatrix_init = []
#for i in range(0,numVertices): # initializing the adjacency matrix
#    tempArray = []
#    for j in range(0,numVertices):
#        tempArray.append(0.0000)
#    adjMatrix_init.append(tempArray)
#enumerate(adjMatrix_init,1)

# 07/06/2018 - to read both subgraphs that will be used
readSubDualGraphs(numVertices_sub1,Sub1_Graphs) # read the eigen values and adjacency matrices for dual graphs for Sub1
readSubDualGraphs(numVertices_sub2,Sub2_Graphs) # read the eigen values and adjacency matrices for dual graphs for Sub2

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
#print "Total number of possible edge combinations generated for one vertex: %d"%len(edgesCombo)
#print the edge combinations for testing
#for i in range(0,len(edgesCombo)):
#    for j in range(0,numVertices_sub1):
#        print edgesCombo[i][j],
#    print

# Enumerate all graphs
enumerateAll()

end_time = time.clock()

print "The total time taken in seconds is:",
print (end_time - start_time)
