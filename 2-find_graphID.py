#Cigdem 2017
#Reads submatrices and determines subgraph IDs (using the same fnc from dualgraphs.py)

import sys
import os
from numpy import *
import numpy.linalg as LA
from decimal import *
from copy import deepcopy
from itertools import permutations

#keys = []
#values = []
#graphID = []

#04/20/2018 - S.J. to make the assigning of graph IDs faster
class DualGraph:

	def __init__(self,a,b):
		self.graphID = a
		self.values = b

	def match(self,b):
		if self.values == b:
			return self.graphID
		else:
			return "NA"


def loadEigenvalues(num_vertices):
	decimalPlace = Decimal("0.00000001")
	f = open("%dEigen" %(num_vertices))
	line = f.readline()
	#keys.append(line[1:-1])
	key = line[1:-1] # 04/20/2018 - S.J.
	while(len(line) > 0):
		tArray = []
		while(len(line) > 0 and line[0] != '>'):
			#tArray.append(line[:-1])
                	tArray.append(Decimal(str(line[:-1])).quantize(decimalPlace))
			line = f.readline()
		#values.append(tArray) # 4/20/2018 - S.J.
		if num_vertices == 2:
			Vertex2Graphs.append(DualGraph(key,tArray))
		elif num_vertices == 3:
			Vertex3Graphs.append(DualGraph(key,tArray))
		elif num_vertices == 4:
			Vertex4Graphs.append(DualGraph(key,tArray))
		elif num_vertices == 5:
			Vertex5Graphs.append(DualGraph(key,tArray))
		elif num_vertices == 6:
			Vertex6Graphs.append(DualGraph(key,tArray))
		elif num_vertices == 7:
			Vertex7Graphs.append(DualGraph(key,tArray))
		elif num_vertices == 8:
			Vertex8Graphs.append(DualGraph(key,tArray))
		elif num_vertices == 9:
			Vertex9Graphs.append(DualGraph(key,tArray))
		if(len(line) > 0):
			#keys.append(line[1:-1])	
			key = line[1:-1] # 04/20/2018 - S.J.
		line = f.readline()	
	f.close()


#04/20/2018 - S.J.
Vertex2Graphs = []
Vertex3Graphs = []
Vertex4Graphs = []
Vertex5Graphs = []
Vertex6Graphs = []
Vertex7Graphs = []
Vertex8Graphs = []
Vertex9Graphs = []


#04/20/2018 - S.J.
loadEigenvalues(2)
loadEigenvalues(3)
loadEigenvalues(4)
loadEigenvalues(5)
loadEigenvalues(6)
loadEigenvalues(7)
loadEigenvalues(8)
loadEigenvalues(9)

myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/list_bp","r").readlines()]
outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/all_subgraph_IDs_test.txt","w")
#outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/all_block_IDs.txt","w")
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/all_structures.txt","r").readlines()]
#outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/all_subgraph_IDs.txt","w")
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/all_structures.txt","r").readlines()]
#outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/all_subgraph_IDs.txt","w")


for file1 in myfiles:
    #print file1
    #matrix_files=os.listdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s/"%file1)
    #matrix_files=os.listdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/submatrices/%s/"%file1)
    #matrix_files=os.listdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/submatrices/%s/"%file1)
    matrix_files=os.listdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices_subgraphs/%s/"%file1)
    count=0 #keep track for number of subgraphs
    total_count=len(matrix_files) # S.J. 10/26/2017 - to get the correct order of matrix reading
    #print total_count
    #for mymatrix in matrix_files:
    for x in range(1,total_count+1): # S.J. 10/26/2017 - to get the correct order of matrix reading
	mymatrix="matrix%d.txt"%(x)
        #print mymatrix
	#adjMatrix=loadtxt("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s/%s"%(file1, mymatrix),dtype='i')
        #adjMatrix=loadtxt("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/submatrices/%s/%s"%(file1, mymatrix),dtype='i')
        #adjMatrix=loadtxt("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/submatrices/%s/%s"%(file1, mymatrix),dtype='i')
	adjMatrix=loadtxt("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices_subgraphs/%s/%s"%(file1, mymatrix),dtype='i')
	count+=1
	#We need to calculate the laplacian matrix from adj matrix. L=D-A
        laplacian = []
        i=0
        for row in adjMatrix:
            #print row
            laprow=-row #take the negative
            laprow[i]+=sum(row) #calc the diagonal elements
            laplacian.append(laprow)
            i+=1


        #print adjMatrix
        #print laplacian

        N=len(adjMatrix)
        if N==1:
            outfile.write("%s\t%d\t%s\n"%(file1, count, "1_1"))
            print file1, count, "1_1"
        elif N>9:
            outfile.write("%s\t%d\t%d\n"%(file1, count, N))
            print file1, count, N
        else:
            #keys=[] 04/20/2018 - S.J.
            #values=[]
            #loadEigenvalues(N)
    
            eigen = sort(LA.eigvals(laplacian))
	    #print eigen
            decimalArray = []
            #decimalPlace = Decimal("0.0001")
	    decimalPlace = Decimal("0.00000001")
            for i in eigen:
		if isinstance(i, complex): # 04/20/2018 - S.J.
			print "I have a complex eigen value"
                	decimalArray.append(Decimal(str(i.real)).quantize(decimalPlace))
		else:
                	decimalArray.append(Decimal(str(i)).quantize(decimalPlace))
	    # 04/20/2018 - S.J.
	    loc = "NA"
	    if N == 2:
		for g in range(0,len(Vertex2Graphs)):
			loc = Vertex2Graphs[g].match(decimalArray)
			if loc != "NA":
				break
	    elif N == 3:
		for g in range(0,len(Vertex3Graphs)):
			loc = Vertex3Graphs[g].match(decimalArray)
			if loc != "NA":
				break
	    elif N == 4:
		for g in range(0,len(Vertex4Graphs)):
			loc = Vertex4Graphs[g].match(decimalArray)
			if loc != "NA":
				break
	    elif N == 5:
		for g in range(0,len(Vertex5Graphs)):
			loc = Vertex5Graphs[g].match(decimalArray)
			if loc != "NA":
				break
	    elif N == 6:
		for g in range(0,len(Vertex6Graphs)):
			loc = Vertex6Graphs[g].match(decimalArray)
			if loc != "NA":
				break
	    elif N == 7:
		for g in range(0,len(Vertex7Graphs)):
			loc = Vertex7Graphs[g].match(decimalArray)
			if loc != "NA":
				break
	    elif N == 8:
		for g in range(0,len(Vertex8Graphs)):
			loc = Vertex8Graphs[g].match(decimalArray)
			if loc != "NA":
				break
	    elif N == 9:
		for g in range(0,len(Vertex9Graphs)):
			loc = Vertex9Graphs[g].match(decimalArray)
			if loc != "NA":
				break

	    # 04/20/2019 - commented out by S.J.
            #loc = -1
            #for i in range(0,len(values)):
            #    tArray = []
            #    for j in range(0,len(values[i])):
            #        tArray.append(Decimal(str(values[i][j])).quantize(decimalPlace))
            #    if decimalArray == tArray:
            #        loc = i

            #address negative 0 output - 04/20/2018 - next two for loops commented out by S.J. 
            #for i in range(0,len(decimalArray)):
            #    if str(decimalArray[i])[0] == "-":
            #        decimalArray[i] = Decimal(str(decimalArray[i])[1:]).quantize(decimalPlace)
            #evNum = 1
            #for i in decimalArray:
                #print "Eigenvalue %d: " %(evNum) + str(i)
            #    evNum+= 1
            #print "%s" %(keys[loc])
	    #if loc == -1: # added S.J. 11/09/2017 to make sure we don't assign graph IDs my mistake even if we don't have them in the library
	    if loc == "NA": # 04/20/2018 - modified by S.J. - added S.J. 11/09/2017 to make sure we don't assign graph IDs my mistake even if we don't have them in the library
                outfile.write("%s\t%d\t%d\n"%(file1, count, N)) #just write the number of vertices for this as well
                print file1, count, N
            else: # we have a graph ID that matches the eigen values in the library
                # 04/20/2019 - commented by S.J.
		#graphID=(keys[loc])
                #print file1, count, graphID #print filename (rnastrand ID), subgraph number, and graph ID
                #outfile.write("%s\t%d\t%s\n"%(file1, count, graphID)) #print filename (rnastrand ID), subgraph number, and graph ID
                print file1, count, loc #print filename (rnastrand ID), subgraph number, and graph ID
                outfile.write("%s\t%d\t%s\n"%(file1, count, loc)) #print filename (rnastrand ID), subgraph number, and graph ID
outfile.close()
