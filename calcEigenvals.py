#!/usr/local/bin/python
import sys
from numpy import *
import numpy.linalg as LA
from decimal import *
from copy import deepcopy
from itertools import permutations

adjMatrix = []
degMatrix = []
tempArray = []
# to re-calculate the iegn values of the dual graph libraries with higher precision. Number of vertices will have to be changed
graph_num=0
size=0
degree=0
file = open("V7AdjDG",'r')
line = file.readline()
for line in file:
	if line != '\n':
		size +=  1
		tempArray = []
		degree = 0
		for x in line.split():
			degree += float(x)
			tempArray.append(float(x))
			#if graph_num == 144:
			#	print x
		adjMatrix.append(tempArray)
		tempArray = []
		for i in range (1,size):
			tempArray.append(0.0000)
		tempArray.append(degree)
		for i in range (size,7):
			tempArray.append(0.0000)
		degMatrix.append(tempArray)

		if size == 7:
			graph_num += 1
			print ">7_%d" %(graph_num)
			#print adjMatrix
			#print degMatrix
			laplacian = array(degMatrix) - array(adjMatrix)
			#print laplacian
			eigen = sort(LA.eigvals(laplacian))
			try:	
				decimalArray = []
                		decimalPlace = Decimal("0.00000001")
               			for i in eigen:
                        		decimalArray.append(Decimal(str(i.real)).quantize(decimalPlace))
				for i in range(0,len(decimalArray)):
                        		if str(decimalArray[i])[0] == "-":
                                		decimalArray[i] = Decimal(str(decimalArray[i])[1:]).quantize(decimalPlace)
					print '{0:.8f}'.format(decimalArray[i])
			except:
				print "Weird eigen values"
				print eigen
			adjMatrix = []
			degMatrix = []
			laplacian = []
			size = 0
file.close()
