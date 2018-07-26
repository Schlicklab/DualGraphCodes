#! /usr/local/bin/python
import sys
from numpy import *
import numpy.linalg as LA
from decimal import *
from copy import deepcopy
from itertools import permutations


keys = []
values = []
graphID = []
def loadEigenvalues(num_vertices):
	f = open("%dEigen" %(num_vertices-1))
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



class Base:
	index = None  #nucleotide Index
	indexBP = None  #base paired to 
	nt = None #NT value
	active = None
	helixNumber = 0 #what helix is this a part of?
	def initialize(self,indexv,ntv,indexBPv): 
		self.index = int(indexv)
		self.indexBP = int(indexBPv)
		self.nt = str(ntv)
		self.active = True

class Loop:
	start = None
	end = None
	def __init__(self):
		pass

class Helix:
	start = None
	end = None
	flag = "" #ARE YOU A LOOP?
	Loop = None
	connected = None
	edges = 0
	def __init__(self):
		self.connected = []

# S.J. 02/11/2018 - to keep track of edges
class Edge:
	Helix1 = None
	Helix2 = None
	start = None
	end = None
	def __init__(self):
		pass

class RNAInfo:
	Bases = None
	Loops = None
	Helices = None
	Edges = None # S.J. 02/11/2018 to keep track of edges
	numVert = None
	adjMatrix = []
	degMatrix = []
	laplacian = None
	def __init__(self):
		self.Bases = [0]
		self.Loops = [0]
		self.Helices = [0]
		self.Edges = [0] # S.J. 02/11/2018 to keep track of edges
		self.numVert = 0
	def makeMatrices(self):
		for i in range(1,len(self.Helices)):
			tArray = []
			for j in range(1,len(self.Helices)):
				tArray.append(0)
			self.adjMatrix.append(tArray)
		for i in range(1,len(self.Helices)):
			tArray = []
			for j in range(1,len(self.Helices)):
				tArray.append(0)
			self.degMatrix.append(tArray)

	def addBase(self,baseA):
		self.Bases.append(baseA)
	def printOut(self,whichBase=1000):
		if whichBase == 1000:
			for i in range(1,len(self.Bases)):
				print "%d\t%d\t%s\t%d" %(self.Bases[i].index,self.Bases[i].indexBP,self.Bases[i].nt,self.Bases[i].helixNumber)
		for i in range(1,len(self.Helices)):
			print "for helix %d: start=%d, end=%d, flag=%s" %(i,self.Helices[i].start,self.Helices[i].end,self.Helices[i].flag)
		for i in range(1,len(self.Loops)):
			print "for loop %d: start=%d, end=%d" %(i,self.Loops[i].start,self.Loops[i].end)
	def printConnections(self):
		for i in range(1,len(self.Helices)):
			print "helix %d is connected to: %s and has %d edges." %(i,str(self.Helices[i].connected),self.Helices[i].edges)
	def printAdj(self):
		print "Adjacency Matrix:"
		for i in self.adjMatrix:
			print i
	def printDeg(self):
		print "Degree Matrix:"
		for i in self.degMatrix:
			print i
	def printLpl(self):
		print "Laplacian Matrix:"
        	for i in self.laplacian:
			print i
	def printHelices(self):
        	for i in range(1,len(self.Helices)):
                #print "Vertex %d: start_pos=%d, end_pos=%d, flag=%s" %(i,self.Helices[i].start,self.Helices[i].end,self.Helices[i].flag)
                #print "Vertex %d: start_pos: (%d, %d), end_pos: (%d, %d)" %(i,self.Helices[i].start,self.Bases[self.Helices[i].start].indexBP,self.Helices[i].end,self.Bases[self.Helices[i].end].indexBP)
                	#print "Vertex %d: first strand: (%d, %d), second strand: (%d, %d)" %(i,self.Helices[i].start,self.Helices[i].end,self.Bases[self.Helices[i].end].indexBP,self.Bases[self.Helices[i].start].indexBP)
                	print "Vertex %d: first strand: %d %d second strand: %d %d" %(i,self.Helices[i].start,self.Helices[i].end,self.Bases[self.Helices[i].end].indexBP,self.Bases[self.Helices[i].start].indexBP)
    
	# S.J. 02/11/2018 - to print edges information
	def printEdges(self):
		for i in range(1,len(self.Edges)):
			print "Edge: helix 1: %d helix 2: %d strand: %d %d" %(self.Edges[i].Helix1,self.Edges[i].Helix2,self.Edges[i].start,self.Edges[i].end)

	def printOrder(self,vOrder):
		order = []
		prevHelix = 0
		for i in range(1,len(self.Bases)):
			currHelix=self.Bases[i].helixNumber
			if currHelix != 0 and currHelix != prevHelix:
				prevHelix = currHelix
				if currHelix != 0:
					order.append(vOrder[currHelix-1])
		print "5'-" + str(order) + "-3'"

	def clear(self):
		Bases = None
		Loops = None
		Helices = None
		numVert = None
		adjMatrix = []
		degMatrix = []
		laplacian = None
		

### makeMatrices ####
#####################
def makeMatrices(RNA):
	for i in range(1,len(RNA.Helices)):
		tArray = []
		for j in range(1,len(RNA.Helices)):
			tArray.append(0)
		RNA.adjMatrix.append(tArray)
	for i in range(1,len(RNA.Helices)):
		tArray = []
		for j in range(1,len(RNA.Helices)):
			tArray.append(0)
		RNA.degMatrix.append(tArray)
		
		
#Translate information from the CT file into an RNA class
def getCTInfo(arg):
	f = open(arg)
	RNA = RNAInfo()
	line = f.readline()
	while(line.split()[0] != '1'):
		line = f.readline()
	while(len(line.split()) > 1):
		oneBase = Base()
		oneBase.initialize(line.split()[0],line.split()[1],line.split()[4])
		RNA.addBase(oneBase)
		line = f.readline()
	f.close()	
	return RNA
##Translate information from BPSEQ file into an RNA class
def getBPSEQInfo(arg):
	f = open(arg)
	RNA = RNAInfo()
	line = f.readline()
	while(line.split()[0] != '1'):
		line = f.readline()
	while(len(line.split()) > 1):
		oneBase = Base()
		oneBase.initialize(line.split()[0],line.split()[1],line.split()[2])
		RNA.addBase(oneBase)
		line = f.readline()
	f.close()	
	return RNA

##Translate information from an adjacency matrix into an RNA class - S.J. 07/05/2018
def getAdjMatInfo(arg):
        f = open(arg)
        RNA = RNAInfo()
        lines = f.readlines()
        for i in range(0,len(lines)):
		RNA.Helices.append(Helix())
                line = lines[i]
                tempArray = []
                degree = 0
                for x in line.split():
                        degree += float(x)
                        tempArray.append(float(x))
                RNA.adjMatrix.append(tempArray)
                tempArray = []
                for j in range (1,i+1):
                        tempArray.append(0.0000)
                tempArray.append(degree)
                for j in range (i+1,len(lines)):
                        tempArray.append(0.0000)
                RNA.degMatrix.append(tempArray)
	f.close()
        return RNA

## Translate information from a dot bracket notation into an RNA class - S.J. 07/25/2018
# first line of a dot bracket notation should contain the sequence, and the second line should contain the dot bracket notation
def getDotBracketInfo(arg):
    f = open(arg)
    RNA = RNAInfo()
    
    line = f.readline().strip()
    sequence = [c for c in line] # stores the base identity of the residue
    line = f.readline().strip()
    dotb = [c for c in line]
    f.close()
    
    stack_bp=[]
    base_pair=[] # to store the residue number of base pair
    
    for i in range(0,len(sequence)):
        base_pair.append(0) # initializing base pairs with all 0's
    
    for i in range(0,len(sequence)): # looping over number of bases based on sequence, as the dot bracket line can contain extra things (from Vienna format)
        if dotb[i] == "(": # opening bracket, add the res number to stack
            stack_bp.append(i)
        elif dotb[i] == ")": # closing bracket, then pop the res number from the stack, and assign base_pairs to both residues
            i_bp = stack_bp.pop()
            base_pair[i]=i_bp+1 #as numbers start from 1
            base_pair[i_bp]=i+1

    #add bases to RNA
    for i in range(0,len(sequence)):
        oneBase = Base()
        oneBase.initialize(i+1,sequence[i],base_pair[i])
        RNA.addBase(oneBase)
    
    return RNA


#Determine whether or not there are pseudoknots in the structure
def pseudoKnots(RNA):
	for i in range(1,len(RNA.Bases)-1):
		if RNA.Bases[i].indexBP > 0:
			for j in range(i+1,len(RNA.Bases)):
				if RNA.Bases[j].indexBP > 0:
					if (j < RNA.Bases[i].indexBP and RNA.Bases[i].indexBP < RNA.Bases[j].indexBP):
						return True
	return False	

### countHelices ####
#####################
#This method counts the number of helices and loops
def countHelices(RNA):
	nHelix = 1
	i = 1
	#find the first.
	while (RNA.Bases[i].indexBP==0):
		i += 1

	RNA.Bases[i].helixNumber = nHelix
	RNA.Bases[i].active = False
	RNA.Bases[RNA.Bases[i].indexBP].helixNumber = nHelix
	RNA.Bases[RNA.Bases[i].indexBP].active = False
	RNA.Helices.append(Helix())
	RNA.Helices[nHelix].start = i;
	RNA.Helices[nHelix].end = i;
	i+=1
	for j in range(i,len(RNA.Bases)):
		if(RNA.Bases[j].indexBP>0 and RNA.Bases[j].active == True):
            		
			if RNA.Bases[j].indexBP+1 != RNA.Bases[j-1].indexBP:
				nHelix += 1
				RNA.Helices.append(Helix())
				RNA.Helices[nHelix].start = j;
				RNA.Helices[nHelix].end = j;
			RNA.Bases[j].helixNumber = nHelix
			RNA.Bases[j].active = False
			RNA.Bases[RNA.Bases[j].indexBP].helixNumber = nHelix
			RNA.Bases[RNA.Bases[j].indexBP].active = False
			RNA.Helices[nHelix].end = j;
		else:
			if RNA.Bases[j].indexBP==0:
				RNA.Bases[j].helixNumber = 0

	for i in range(1,len(RNA.Helices)):
		helixEnd = RNA.Helices[i].end
		if clearPath(RNA,helixEnd,RNA.Bases[helixEnd].indexBP):
			loop = Loop()
			loop.start = RNA.Helices[i].start
			loop.end = RNA.Bases[RNA.Helices[i].start].indexBP
			RNA.Loops.append(loop)
			RNA.Helices[i].flag = 'L'
			RNA.Helices[i].Loop = loop



### changeHelices ####
#####################
#Combines helices if they are only separated by one unpaired NT
def changeHelices(RNA):	
	changes = []
	for i in range(1,len(RNA.Helices)-1):  
		#never do this to loops
		if RNA.Helices[i].flag == 'L' and RNA.Helices[i+1].flag == 'L':
			pass
		else:
			helix2fiveStart = RNA.Helices[i+1].start
			helix2fiveEnd = RNA.Helices[i+1].end
			helix2threeEnd = RNA.Bases[RNA.Helices[i+1].start].indexBP
			helix2threeStart = RNA.Bases[RNA.Helices[i+1].end].indexBP
			helix1fiveEnd = RNA.Helices[i].end
			helix1fiveStart = RNA.Helices[i].start
			helix1threeStart = RNA.Bases[RNA.Helices[i].end].indexBP
			helix1threeEnd = RNA.Bases[RNA.Helices[i].start].indexBP
			# if statement added by S.J. 01/30/2018 to not combine helices that don't have base pairs in the correct order
			if helix2threeEnd > helix1threeStart: # 3' of helix 1 starts before 3' of helix 2 ends, therefore cannot be combined
				continue
			Total5P = abs(helix2fiveStart - helix1fiveEnd)-1
			Total3P = abs(helix1threeStart - helix2threeEnd)-1
			if ((abs(Total5P + Total3P) < 2) or (abs(Total5P) == 1 and abs(Total3P) == 1)):
				changes.append(i)
	for i in changes: #change bases
		j = 1
		##Base Change
		while(RNA.Bases[j].helixNumber <=i):
			j += 1
		for k in range(j,len(RNA.Bases)):
			if RNA.Bases[k].helixNumber != 0 and RNA.Bases[k].helixNumber>i:
				RNA.Bases[k].helixNumber -= 1
		
		RNA.Helices[i].end = RNA.Helices[i+1].end
		if RNA.Helices[i+1].flag == 'L':
			RNA.Helices[i].flag = 'L'
			RNA.Helices[i].Loop = RNA.Helices[i+1].Loop 
			RNA.Helices[i].Loop.start = RNA.Helices[i].start 
			RNA.Helices[i].Loop.end = RNA.Helices[i].end 
		del RNA.Helices[i+1]
		for m in range(0,len(changes)):
			if changes[m] > i:
				changes[m] -= 1

	
	singleHelices = []
	for i in range(1,len(RNA.Helices)):
		if RNA.Helices[i].start == RNA.Helices[i].end:
			singleHelices.append(i)
			fivePrime = RNA.Helices[i].start
			threePrime  = RNA.Bases[fivePrime].indexBP
			print "Helix %d is a single base-pair helix with 5' = %d and 3' = %d!" %(i,fivePrime,threePrime)
	for i in singleHelices:
		fivePrime = RNA.Helices[i].start
		threePrime  = RNA.Bases[fivePrime].indexBP
		RNA.Bases[fivePrime].indexBP = 0
		RNA.Bases[fivePrime].helixNumber = 0
		RNA.Bases[threePrime].indexBP = 0
		RNA.Bases[threePrime].helixNumber = 0
		
		j = 1
		while(j<len(RNA.Bases) and RNA.Bases[j].helixNumber <=i):
			j+=1
			for k in range(j,len(RNA.Bases)):
				if RNA.Bases[k].helixNumber != 0 and RNA.Bases[k].helixNumber>i:
					RNA.Bases[k].helixNumber -= 1
		del RNA.Helices[i]
		for m in range(0,len(singleHelices)):
			if singleHelices[m] > i:
				singleHelices[m] -= 1
	
	#redo loops if you removed single helices
	if len(singleHelices) > 0:
		for i in range(1,len(RNA.Helices)):
			helixEnd = RNA.Helices[i].end
			if clearPath(RNA,helixEnd,RNA.Bases[helixEnd].indexBP):
				loop = Loop()
				loop.start = RNA.Helices[i].start
				loop.end = RNA.Bases[RNA.Helices[i].start].indexBP
				RNA.Loops.append(loop)
				RNA.Helices[i].flag = 'L'
				RNA.Helices[i].Loop = loop

	

### connectHelices ####
#####################
#method to count the number of connections 
def connectHelices(RNA):	
#first connect loops to themselves
	nEdges=0; # S.J. 02/11/2018 - count edges
	for i in range(1,len(RNA.Helices)):
		if RNA.Helices[i].flag == 'L':
			RNA.Helices[i].edges += 2
			RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
			nEdges += 1
			RNA.Edges[nEdges].Helix1 = i
			RNA.Edges[nEdges].Helix2 = i
			RNA.Edges[nEdges].start = RNA.Helices[i].end
			RNA.Edges[nEdges].end = RNA.Bases[RNA.Helices[i].end].indexBP
			RNA.adjMatrix[i-1][i-1] = 2
	for i in range(1,len(RNA.Helices)-1):
		for j in range(i+1,len(RNA.Helices)):
			helix2fiveStart = RNA.Helices[j].start
			helix2fiveEnd = RNA.Helices[j].end
			helix2threeEnd = RNA.Bases[RNA.Helices[j].start].indexBP
			helix2threeStart = RNA.Bases[RNA.Helices[j].end].indexBP
			helix1fiveEnd = RNA.Helices[i].end
			helix1fiveStart = RNA.Helices[i].start
			helix1threeStart = RNA.Bases[RNA.Helices[i].end].indexBP
			helix1threeEnd = RNA.Bases[RNA.Helices[i].start].indexBP
            

    

			helix2 = [helix2fiveStart, helix2fiveEnd, helix2threeEnd, helix2threeStart]
			helix1 = [helix1fiveStart, helix1fiveEnd, helix1threeEnd, helix1threeStart]
                      
			

			if (clearPath(RNA,helix1fiveEnd,helix2fiveStart) or (helix2fiveStart - helix1fiveEnd)==1):
				increment(RNA,i,j)
				RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
				nEdges += 1
				RNA.Edges[nEdges].Helix1 = i
				RNA.Edges[nEdges].Helix2 = j
				RNA.Edges[nEdges].start = helix1fiveEnd
				RNA.Edges[nEdges].end = helix2fiveStart

			if (clearPath(RNA,helix2threeEnd,helix1threeStart) or (helix1threeStart - helix2threeEnd)==1):
				increment(RNA,i,j)
				RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
				nEdges += 1
				RNA.Edges[nEdges].Helix1 = i
				RNA.Edges[nEdges].Helix2 = j
				RNA.Edges[nEdges].start = helix2threeEnd
				RNA.Edges[nEdges].end = helix1threeStart

			if (clearPath(RNA,helix1fiveEnd,helix2threeStart) or (helix2threeStart - helix1fiveEnd)==1):
				increment(RNA,i,j)
				RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
				nEdges += 1
				RNA.Edges[nEdges].Helix1 = i
				RNA.Edges[nEdges].Helix2 = j
				RNA.Edges[nEdges].start = helix1fiveEnd
				RNA.Edges[nEdges].end = helix2threeStart

			if (clearPath(RNA,helix2threeEnd,helix1fiveStart) or (helix1fiveStart - helix2threeEnd)==1):
				increment(RNA,i,j)
				RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
				nEdges += 1
				RNA.Edges[nEdges].Helix1 = i
				RNA.Edges[nEdges].Helix2 = j
				RNA.Edges[nEdges].start = helix2threeEnd
				RNA.Edges[nEdges].end = helix1fiveStart

			if (clearPath(RNA,helix1threeEnd,helix2fiveStart) or (helix2fiveStart -helix1threeEnd==1)):
				increment(RNA,i,j)
				RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
				nEdges += 1
				RNA.Edges[nEdges].Helix1 = i
				RNA.Edges[nEdges].Helix2 = j
				RNA.Edges[nEdges].start = helix1threeEnd
				RNA.Edges[nEdges].end = helix2fiveStart

			if pseudoKnots(RNA):
				if (clearPath(RNA,helix2fiveEnd,helix1threeStart)):
					increment(RNA,i,j)
					RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
					nEdges += 1
					RNA.Edges[nEdges].Helix1 = i
					RNA.Edges[nEdges].Helix2 = j
					RNA.Edges[nEdges].start = helix2fiveEnd
					RNA.Edges[nEdges].end = helix1threeStart

				if (clearPath(RNA,helix1threeEnd,helix2threeStart)):
					increment(RNA,i,j)
					RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
					nEdges += 1
					RNA.Edges[nEdges].Helix1 = i
					RNA.Edges[nEdges].Helix2 = j
					RNA.Edges[nEdges].start = helix1threeEnd
					RNA.Edges[nEdges].end = helix2threeStart

                #Added by CSB:
				if (clearPath(RNA,helix1threeStart,helix2fiveStart)):
   					increment(RNA,i,j)
					RNA.Edges.append(Edge()) # S.J. 02/11/2018 - store edge information
					nEdges += 1
					RNA.Edges[nEdges].Helix1 = i
					RNA.Edges[nEdges].Helix2 = j
					RNA.Edges[nEdges].start = helix1threeStart
					RNA.Edges[nEdges].end = helix2fiveStart


    
    
	for m in range(1,len(RNA.Helices)):
		RNA.degMatrix[m-1][m-1] = RNA.Helices[m].edges

def correctHNumbers(RNA):
	for i in range(1,len(RNA.Helices)):
		for j in range(RNA.Helices[i].start,RNA.Helices[i].end+1):
			#if RNA.Bases[j].helixNumber == 0:
				RNA.Bases[j].helixNumber = i
		for l in range(RNA.Bases[RNA.Helices[i].end].indexBP,RNA.Bases[RNA.Helices[i].start].indexBP+1):		
			#if RNA.Bases[l].helixNumber == 0:
				RNA.Bases[l].helixNumber = i

#def calcEigen(RNA,arg):
def calcEigen(RNA): # S.J. 07/05/2018 - removing the last argument as that is not needed
	if len(RNA.Helices)==2:
		print "1_1" 
	elif len(RNA.Helices)>10:
		print "TMV,%d" %(len(RNA.Helices)-1)
	else:
		loadEigenvalues(len(RNA.Helices))
		RNA.laplacian = array(RNA.degMatrix) - array(RNA.adjMatrix)
		RNA.printLpl()
		eigen = sort(LA.eigvals(RNA.laplacian))
		decimalArray = []
		#decimalPlace = Decimal("0.0001")
		decimalPlace = Decimal("0.00000001")
		for i in eigen:
			decimalArray.append(Decimal(str(i)).quantize(decimalPlace))
		loc = -1
		for i in range(0,len(values)):
			tArray = []
			for j in range(0,len(values[i])):
				tArray.append(Decimal(str(values[i][j])).quantize(decimalPlace))
			if decimalArray == tArray:
				loc = i
                		#print decimalArray #added by S.J. 11/09/2017 for testing
                		#print tArray
                		#print loc
                		#print keys[loc]
		#address negative 0 output
		for i in range(0,len(decimalArray)):
			if str(decimalArray[i])[0] == "-":	
				decimalArray[i] = Decimal(str(decimalArray[i])[1:]).quantize(decimalPlace)
		evNum = 1
		for i in decimalArray:
			print "Eigenvalue %d: " %(evNum) + str(i)
			evNum+= 1
		#print loc  #added by S.J. 11/09/2017 for testing
		#print len(keys)
		#print keys[-2]
		#print keys[0]
		#print keys[1]
		#print keys[2]
		#print keys[38594]
		#print keys[38595]
		#the following if statement added by S.J. 11/09/2017 to make sure that is we don't find the eigen values then we accidently don't assign a graph ID, because python does not give errors for negative indices
		if loc == -1: # if matching graph is not found
			print "TMV,%d" %(len(RNA.Helices)-1)
			return 0
		print "Graph ID: %s" %(keys[loc])
		#print "<a href=http://www.biomath.nyu.edu/rag/dual_topology.php?topo=%s>Graph ID: %s</a>" %(keys[loc], keys[loc])
		graphID.append(keys[loc])
		return 1 # added by S.J. 11/09/2017 to return 1 if successful in assigning graph ID

def label(RNA):
		#Code to LABEL, please address any PROBLEMS!###ALSO FIX GRAPHID
		ID = graphID[0]
		vertexOrder = None
		vNum = int(ID.split('_')[0])
		ID = int(ID.split('_')[1])
		g = open("V%dAdjDG" %(vNum),'r')
		a = []
		#print ((ID-1)*(vNum+2)+2)
		for j in range(0,((ID-1)*(vNum+2)+2)-1):
			g.readline()
		for k in range(0,vNum):
			tempA = []
			for l in g.readline().split():
				tempA.append(int(float(l)))
			a.append(tempA)
		g.close()
		b = RNA.adjMatrix
		#print a
		#print b
		c = deepcopy(a)
		num = []
		for i in range(0,len(a)):
			num.append(i)
		for i in list(permutations(num)):
			listI = list(i) 
			for j in range(0,len(c)):
				jI = listI[j]
				for k in range(0,len(c)):
					kI= listI[k]
					c[j][k] = a[jI][kI]
			if c==b:
				for i in range(0,len(listI)):
					listI[i]+=1
				vertexOrder = listI
				break
		#print str(vertexOrder)		
		if vertexOrder != None:			
			RNA.printOrder(vertexOrder)
		else:
			print "Graph isomorphism."		
		

def main():			
	#for arg in sys.argv[1:]: # S.J. 07/05/2018 - eliminating the for loop, so this will run on one file at a time
	# S.J. 07/04/2018 - adding to read in adjacency matrices
	adjMatTrue = False
	if sys.argv[1] == "-adj_mat":
		adjMatTrue = True
                RNA = getAdjMatInfo(sys.argv[2])
    	elif sys.argv[1] == "-dotb": # S.J. 07/25/2018 - flag to read RNA info in dotbracket notation
        	RNA = getDotBracketInfo(sys.argv[2])
	elif sys.argv[1][-2:] == "ct": 
		RNA = getCTInfo(sys.argv[1])
	else:
		RNA = getBPSEQInfo(sys.argv[1])

	if not adjMatTrue:
		countHelices(RNA) 
		changeHelices(RNA)
		RNA.makeMatrices()
		connectHelices(RNA)
	print "Number of Vertices: " + str(len(RNA.Helices)-1)
	RNA.printAdj()
        RNA.printDeg()
	
	if not adjMatTrue:
        	RNA.printHelices()
       		RNA.printEdges() # S.J. 02/11/2018

        success=calcEigen(RNA) # catching the return value S.J. 11/09/2017 - removing the last argument as that is not need 07/05/2018
        
	if not adjMatTrue:
		correctHNumbers(RNA)

        if len(RNA.adjMatrix)==1 or len(RNA.adjMatrix)>9:
        	print "No matching graph exists because vertex number is either 1 or greater than 10."
	elif success == 0: # no graph ID was assigned as eigen values not in the library S.J. 11/09/2017
		print "No matching graph exists (even if the vertex number is between 2 and 9)."
        else:
		if not adjMatTrue:
        		label(RNA)

#check if there are only zeroes between nt (start) and nt (end)
def clearPath(RNA,start,end):
	if end<start:
		return False
	for i in range(start+1,end):
		"false for %d and %d" %(start,end)
		if RNA.Bases[i].indexBP !=  0:
			return False
	return True

def increment(RNA,i,j):
	RNA.Helices[i].edges += 1
	RNA.Helices[j].edges += 1
	RNA.adjMatrix[i-1][j-1] += 1
	RNA.adjMatrix[j-1][i-1] += 1

if __name__ == "__main__":
    main()

