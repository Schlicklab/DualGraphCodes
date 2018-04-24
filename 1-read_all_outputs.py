##Cigdem S Bayrak
##December, 2016

#Reads output file of the dual graph partitioning code (python & c++ from run_code.py)
#Creates directory for each structure
#Writes each submatrix (adj matrix of each subgraph) under submatrices directory
#Outputs edges.txt and subgraphs.png into subgraph/name_of_the_structure/ directory

import re
import numpy
import random
from igraph import *
import os

#read dual_graph_partitioning output

myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/list_bp_full","r").readlines()]
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/all_structures.txt","r").readlines()]
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/all_structures.txt","r").readlines()]

for filename in myfiles:
    print filename
    if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s"%filename):
        os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s"%filename)
    #if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/submatrices/%s"%filename):
    #    os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/submatrices/%s"%filename)
    #if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/submatrices/%s"%filename):
        #os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/submatrices/%s"%filename)
    #if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices_subgraphs/%s"%filename):
    #    os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices_subgraphs/%s"%filename)
    
    if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/edges/%s"%filename):
        os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/edges/%s"%filename)
    #if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/edges/%s"%filename):
    #    os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/edges/%s"%filename)
    #if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/edges/%s"%filename):
        #os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/edges/%s"%filename)
    #if not os.path.isdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/edges_subgraphs/%s"%filename):
    #    os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/edges_subgraphs/%s"%filename)

    inputfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/output/%s_output.txt"%filename,"r")
    #inputfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/output/%s_output.txt"%filename,"r")
    #inputfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/%s_output.txt"%filename,"r")
    #inputfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/output/%s_Subgraphs.txt"%filename,"r")
    mylines=inputfile.readlines()
    inputfile.close()

    count=0
    for i in range(len(mylines)):
        line=mylines[i]
        if "New Block" in line:
        #if "New Subgraph" in line:
            count+=1 #number of subgraphs
            edges=[x.strip() for x in mylines[i+2].strip().split('-')] #read the next line. Get the edges (pairs) by splitting from "-"
            numbers= [int(x) for x in re.findall(r"[\w']+",mylines[i+2])] #read only the numbers from a line like (11,3) - (10,11) - (4,10) - (3,4) -
            vertices=list(set(numbers)) #get the unique numbers. these are the vertices of the corresponding subgraph. i.e. 3,4,10,11
        
            #print edges
            #print vertices
            matrix = numpy.zeros((len(vertices),len(vertices)),int) #adjacency matrix of the subgraph
            for j in range(len(edges)-1): #last entry is empty as there is a "-" at ythe end of each line. Or simply use del edges[-1]
                edge=edges[j]
                indices=[int(x) for x in re.findall(r"[\w']+",edge)] #read each pair (edge). i.e.  (11,3)
            
                #print indices[0], indices[1]
                m=vertices.index(indices[0]) #determine the order (index) of the first vertex of the edge. For (11,3) it is 11 and the index is 3 according to vertices
                n=vertices.index(indices[1]) ##determine the order (index) of the first vertex of the edge. It is 3 and the index is 0
            
                matrix[m][n]+=1 #increase the number of connections in the adjacency matrix. matrix[3][0] will be increased 1
                matrix[n][m]+=1 #since the matrix is symmetric, increase matrix[0][3] 1.
            #print matrix
            outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s/matrix%d.txt"%(filename,count),"w")
            #outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/submatrices/%s/matrix%d.txt"%(filename,count),"w")
            #outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/submatrices/%s/matrix%d.txt"%(filename,count),"w")
            #outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices_subgraphs/%s/matrix%d.txt"%(filename,count),"w")
            for ii in range(len(matrix)):
                for jj in range(len(matrix[0])):
                    outfile.write("%d\t"%matrix[ii][jj])
                outfile.write("\n")
            outfile.close()
            outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/edges/%s/edges%d.txt"%(filename,count),"w")
            #outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/edges/%s/edges%d.txt"%(filename,count),"w")
            #outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/edges/%s/edges%d.txt"%(filename,count),"w")
            #outfile=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/edges_subgraphs/%s/edges%d.txt"%(filename,count),"w")
            outfile.write(mylines[i+2]) #Write the edges
            outfile.close()
        
#        random.seed(1)
#        g=Graph() #generate the dual graph
#        g=Graph.Adjacency(matrix.tolist(), mode=ADJ_UNDIRECTED)
#        
#        mylayout=g.layout_kamada_kawai()
#        bbox=BoundingBox(3200,3200) #we determine the bbox and margin so that there will be white space around the plot. Otherwise it doesnt show the plot completely.
#        mylabels=vertices #put vertex numbers as labels
#        
#        figure=plot(g, vertex_color="red", vertex_frame_color="red",  edge_color="black", edge_width=10, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(500,500,500,500), vertex_size=80) #vertex_label=mylabels,
#        figure.save("/Users/cs4367/Desktop/Dual_Partitioning/rRNA_ribovision/subgraphs/DM_LSU/subgraph%d.png"%count)











