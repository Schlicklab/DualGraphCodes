#!/usr/local/bin/python

# 07/10/2018 - Swati Jain - (S.J.) - script to generate adjacency matrices for subgraphs created after partitioning and assign dual graph IDs to them
# combination of previous scripts 1-read_all_outputs.py and 2-find_graphID.py that were written by Cigdem Sevim Bayrak (CSB) and modified by S.J.
# takes the input as the partitioning output (should be in the directpry named output in the path given (without the last slash)), needs two directories called edges and submatrices to be already present in the path.

import re
import numpy
import random
import os
from ClassesFunctions import *

def readDualGraphs():

    Graphs=[]
    DualGraphsLib.append(Graphs) # for vertex == 1, no graphs are there
    for i in range(2,10): # will read dual graphs from 2-9 vertices (10 used as the range function stops before the last number)

        Graphs=[]
        file_eigen = "%dEigen"%i
        file_adjMat = "V%dAdjDG"%i

        loadEigenvalues(Graphs,i,file_eigen) # load eigen values for dual graphs for vertex number i
        loadAdjMatrices(Graphs,i,file_adjMat) # load adjacency matrices for dual graphs for vertex number i

        DualGraphsLib.append(Graphs)


# flag to specify if the partitioning only contains the blocks or all possible subgraphs
file_list = sys.argv[2]
path = sys.argv[3]
sub_string = "" # the string used to identify when the partitioning output starts a new block/subgraph
submatrix_dir = "" # directories based on blocks or subgraphs
edges_dir = ""
part_out_dir = ""

if sys.argv[1] == "-blocks":
    sub_string = "New Block"
    sub_outfile=open("%s/all_block_IDs.txt"%path,"w")
    submatrix_dir = "submatrices"
    edges_dir = "edges"
    part_out_dir = "output_correct"
elif sys.argv[1] == "-subgraphs":
    sub_string = "New Subgraph"
    sub_outfile=open("%s/all_subgraph_IDs.txt"%path,"w")
    submatrix_dir = "submatrices_subgraphs"
    edges_dir = "edges_subgraphs"
    part_out_dir = "output_subgraphs"

# read the dual graphs from the existing dual graph libraries (dervied from enumeration)
DualGraphsLib=[]
readDualGraphs()

# read the filenames from the list of files provided (without the .pdb or .bpseq extension)
myfiles=[x.strip() for x in open(file_list,"r").readlines()]

for filename in myfiles:
    #print filename
    
    if not os.path.isdir("%s/%s/%s"%(path,submatrix_dir,filename)):
        os.mkdir("%s/%s/%s"%(path,submatrix_dir,filename))
    if not os.path.isdir("%s/%s/%s"%(path,edges_dir,filename)):
        os.mkdir("%s/%s/%s"%(path,edges_dir,filename))

    #read dual_graph_partitioning output
    inputfile=open("%s/%s/%s_output.txt"%(path,part_out_dir,filename),"r")
    mylines=inputfile.readlines()
    inputfile.close()

    count = 0
    for i in range(len(mylines)):
        line=mylines[i]
        if sub_string in line:

            count+=1 #number of subgraphs
            edges=[x.strip() for x in mylines[i+2].strip().split('-')] #read the next line. Get the edges (pairs) by splitting from "-"
            numbers= [int(x) for x in re.findall(r"[\w']+",mylines[i+2])] #read only the numbers from a line like (11,3) - (10,11) - (4,10) - (3,4) -
            vertices=list(set(numbers)) #get the unique numbers. these are the vertices of the corresponding subgraph. i.e. 3,4,10,11
            #print edges
            #print vertices
            
            matrix = []
            for a in range(0,len(vertices)):
                tempArray = []
                for b in range(0,len(vertices)):
                    tempArray.append(0)
                matrix.append(tempArray)
            
            #matrix = numpy.zeros((len(vertices),len(vertices)),int) #adjacency matrix of the subgraph
            for j in range(len(edges)-1): #last entry is empty as there is a "-" at ythe end of each line. Or simply use del edges[-1]
                edge=edges[j]
                indices=[int(x) for x in re.findall(r"[\w']+",edge)] #read each pair (edge). i.e.  (11,3)
                
                #print indices[0], indices[1]
                m=vertices.index(indices[0]) #determine the order (index) of the first vertex of the edge. For (11,3) it is 11 and the index is 3 according to vertices
                n=vertices.index(indices[1]) ##determine the order (index) of the first vertex of the edge. It is 3 and the index is 0
                
                matrix[m][n]+=1 #increase the number of connections in the adjacency matrix. matrix[3][0] will be increased 1
                matrix[n][m]+=1 #since the matrix is symmetric, increase matrix[0][3] 1.
                
            #print matrix
            #print len(matrix)
            outfile=open("%s/%s/%s/matrix%d.txt"%(path,submatrix_dir,filename,count),"w")
            for ii in range(len(matrix)):
                for jj in range(len(matrix[0])):
                    outfile.write("%d\t"%matrix[ii][jj])
                outfile.write("\n")
            outfile.close()
            outfile=open("%s/%s/%s/edges%d.txt"%(path,edges_dir,filename,count),"w")
            outfile.write(mylines[i+2]) #Write the edges
            outfile.close()

            # calculate eigen values for the current subgraph matrix
            N=len(vertices)
            if N==1:
                sub_outfile.write("%s\t%d\t%s\n"%(filename, count, "1_1"))
                print filename, count, "1_1"
            elif N>9:
                sub_outfile.write("%s\t%d\t%d\n"%(filename, count, N))
                print filename, count, N
            else:
                eigen = calcEigenValues(matrix) # calculate the eigen values for the subgraph matrix
                #printEigenValues(eigen)

                loc = "NA"
                for g in DualGraphsLib[N-1]: # search for matching dual graph in the list of dual graphs with the correct vertex number
                    loc = g.match(eigen,matrix)
                    if loc != "NA": # match found
                        break

                if loc == "NA": # 04/20/2018 - modified by S.J. - added S.J. 11/09/2017 to make sure we don't assign graph IDs my mistake even if we don't have them in the library
                    sub_outfile.write("%s\t%d\t%d\n"%(filename, count, N)) #just write the number of vertices for this as well
                    print filename, count, N
                else:
                    print filename, count, loc #print filename (rnastrand ID), subgraph number, and graph ID
                    sub_outfile.write("%s\t%d\t%s\n"%(filename, count, loc)) #print filename (rnastrand ID), subgraph number, and graph ID

sub_outfile.close()
