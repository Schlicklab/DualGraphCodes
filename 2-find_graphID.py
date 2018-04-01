#Cigdem 2017
#Reads submatrices and determines subgraph IDs (using the same fnc from dualgraphs.py)

import sys
import os
from numpy import *
import numpy.linalg as LA
from decimal import *
from copy import deepcopy
from itertools import permutations

keys = []
values = []
graphID = []
def loadEigenvalues(num_vertices):
	f = open("%dEigen" %(num_vertices))
	line = f.readline()
	keys.append(line[1:-1])
	while(len(line) > 0):
		tArray = []
		while(len(line) > 0 and line[0] != '>'):
			tArray.append(line[:-1])
			line = f.readline()
		values.append(tArray)
		if(len(line) > 0):
			keys.append(line[1:-1])	
		line = f.readline()	
	f.close()



#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/list_bp","r").readlines()]
#outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/all_subgraph_IDs.txt","w")
myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/all_structures.txt","r").readlines()]
outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/all_subgraph_IDs.txt","w")
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/all_structures.txt","r").readlines()]
#outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/all_subgraph_IDs.txt","w")


for file1 in myfiles:
    #print file1
    #matrix_files=os.listdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s/"%file1)
    matrix_files=os.listdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/submatrices/%s/"%file1)
    #matrix_files=os.listdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/submatrices/%s/"%file1)
    count=0 #keep track for number of subgraphs
    total_count=len(matrix_files) # S.J. 10/26/2017 - to get the correct order of matrix reading
    #print total_count
    #for mymatrix in matrix_files:
    for x in range(1,total_count+1): # S.J. 10/26/2017 - to get the correct order of matrix reading
	mymatrix="matrix%d.txt"%(x)
        #print mymatrix
	#adjMatrix=loadtxt("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s/%s"%(file1, mymatrix),dtype='i')
        adjMatrix=loadtxt("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/submatrices/%s/%s"%(file1, mymatrix),dtype='i')
        #adjMatrix=loadtxt("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/submatrices/%s/%s"%(file1, mymatrix),dtype='i')
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
            keys=[]
            values=[]
            loadEigenvalues(N)
    
            eigen = sort(LA.eigvals(laplacian))
            decimalArray = []
            decimalPlace = Decimal("0.0001")
            for i in eigen:
                decimalArray.append(Decimal(str(i)).quantize(decimalPlace))
            loc = -1
            for i in range(0,len(values)):
                tArray = []
                for j in range(0,len(values[i])):
                    tArray.append(Decimal(str(values[i][j])).quantize(decimalPlace))
                if decimalArray == tArray:
                    loc = i
            #address negative 0 output
            for i in range(0,len(decimalArray)):
                if str(decimalArray[i])[0] == "-":
                    decimalArray[i] = Decimal(str(decimalArray[i])[1:]).quantize(decimalPlace)
            evNum = 1
            for i in decimalArray:
                #print "Eigenvalue %d: " %(evNum) + str(i)
                evNum+= 1
            #print "%s" %(keys[loc])
	    if loc == -1: # added S.J. 11/09/2017 to make sure we don't assign graph IDs my mistake even if we don't have them in the library
                outfile.write("%s\t%d\t%d\n"%(file1, count, N)) #just write the number of vertices for this as well
                print file1, count, N
            else: # we have a graph ID that matches the eigen values in the library
                graphID=(keys[loc])
                print file1, count, graphID #print filename (rnastrand ID), subgraph number, and graph ID
                outfile.write("%s\t%d\t%s\n"%(file1, count, graphID)) #print filename (rnastrand ID), subgraph number, and graph ID
outfile.close()
