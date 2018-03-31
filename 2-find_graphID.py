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



myfiles=[x.strip() for x in open("C:/Users/cigdems/Desktop/NYU_WORK/Dual_Partitioning/rnastrand_results/all_structures.txt","r").readlines()]
outfile=open("C:/Users/cigdems/Desktop/NYU_WORK/Dual_Partitioning/rnastrand_results/all_subgraph_IDs.txt","w")                
for file1 in myfiles:
    #print file1
    matrix_files=os.listdir("C:/Users/cigdems/Desktop/NYU_WORK/Dual_Partitioning/rnastrand_results/submatrices/%s/"%file1)
    count=0 #keep track for number of subgraphs
    for mymatrix in matrix_files:
        adjMatrix=loadtxt("C:/Users/cigdems/Desktop/NYU_WORK/Dual_Partitioning/rnastrand_results/submatrices/%s/%s"%(file1, mymatrix),dtype='i')
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
            graphID=(keys[loc])
            print file1, count, graphID #print filename (rnastrand ID), subgraph number, and graph ID
            outfile.write("%s\t%d\t%s\n"%(file1, count, graphID)) #print filename (rnastrand ID), subgraph number, and graph ID
outfile.close()
