#! /usr/local/bin/python
import sys
import re #for regular expression
import os

#code to extract the atomic fragments of subgraphs (as divided the the dual graph partitioning code) and make a database with atomic fragments corresponding to dual graph IDs

class Helix:
    #number = None
    #start_5 = None
    #end_5 = None
    #start_3 = None
    #end_3 = None
    #Residues = None

    def __init__(self,a,b,c,d,e):
        self.number = a
        self.start_5 = b
        self.end_5 = c
        self.start_3 = d
        self.end_3 = e
        self.Residues = [x for x in range(self.start_5,self.end_5+1)] + [x for x in range(self.start_3,self.end_3+1)] # residue numbers in the helix

    def printHelix(self):
        print "Helix: %d %d %d %d %d" %(self.number,self.start_5,self.end_5,self.start_3,self.end_3)
        for i in self.Residues:
            print i,
        print

class Edge:
    #helix1 = None
    #helix2 = None
    #start = None
    #end = None
    #Residues = None

    def __init__(self,a,b,c,d):
        self.helix1 = a
        self.helix2 = b
        self.start = c
        self.end = d
        self.Residues = [x for x in range(self.start,self.end+1)] # residue numbers in the edges

    def printEdge(self):
        print "Edge: %d %d %d %d" %(self.helix1,self.helix2,self.start,self.end)
        for i in self.Residues:
            print i,
        print

class Subgraph:
    #Vertices = None
    #Residues = None
    #graph_id = None

    def __init__(self,a,b,c):
        self.Vertices = a
        self.Residues = list(set(b)) # get unique residue numbers
        self.Residues.sort()
        self.graph_id = c

    def printSubgraph(self,file):
        file.write("Subgraph ID: " + str(self.graph_id) + "\t")
        #print "Vertices: ",
        for v in self.Vertices:
            file.write(str(v-1) + " ") # just for pritning the graph vertices start from 0, as in case of the tree graphs as well
        file.write("\n")
        #for r in self.Residues:
            #file.write(str(r) + " ")
        #file.write("\n")

class Base:

    def __init__(self):
        self.pdb_lines = []

    def printBase(self,file):
        for i in range(0,len(self.pdb_lines)):
            file.write(self.pdb_lines[i])

def main():

    # read the basepairs file

    basepairs_file = "/Users/sj78/Documents/labwork/RNAFilesfromPDB_Aug31/FinalBPSEQs_LC/output/%s_basepairs.txt" %(sys.argv[1])
    output_file = "/Users/sj78/Documents/labwork/RNAFilesfromPDB_Aug31/FinalBPSEQs_LC/output_subgraphs/%s_output.txt" %(sys.argv[1])
    subgraph_file = "/Users/sj78/Documents/labwork/RNAFilesfromPDB_Aug31/FinalBPSEQs_LC/all_subgraph_IDs.txt"
    pdb_file = "/home/sj78/labwork/AllRNADataset_August2018/RNAPDBs/%s.pdb" %(sys.argv[1])
    #bpseq_file = "/home/sj78/labwork/NonRedundantList_Oct2017/BPSEQs/%s.bpseq" %(sys.argv[1])

    Edges = []
    Helices = []
    Subgraphs = []
    Bases = []

    file = open(basepairs_file,"r") # reading the basepairs file
    bplines = [x.strip() for x in file.readlines()]
    file.close()
    for line in bplines:
        words = re.split('\s+',line)
        if words[0] == "Vertex":
            Helices.append(Helix(int(words[1].split(":")[0]),int(words[4]),int(words[5]),int(words[8]),int(words[9])))
        if words[0] == "Edge:":
            Edges.append(Edge(int(words[3]),int(words[6]),int(words[8]),int(words[9])))

    #for i in range(0,len(Helices)):
    #   Helices[i].printHelix()
    #for i in range(0,len(Edges)):
    #   Edges[i].printEdge()

    # read the output of the dual graph partitioning code file and subgraphs id file

    file = open(output_file,"r") # subgraph vertices file
    lines = file.readlines()
    file.close()

    file = open(subgraph_file,"r") # subgraph ids file
    subgraph_lines = [x.strip() for x in file.readlines()]
    file.close()

    file_index = 0
    for i in range(0,len(subgraph_lines)):
        cols = re.split('[\s]+',subgraph_lines[i])
        if cols[0] == sys.argv[1]: # the subgraph ids start from here in this file
            file_index = i
            break

    for i in range(0,len(lines)):
        
        if "New Subgraph" in lines[i]: # contains a new block
            sub_edges=[x.strip() for x in lines[i+2].strip().split('-')] # list of edges in the subgraph (i+2 as the second line after the New Block contains the edges)
            del sub_edges[-1]
            numbers = [int(x)+1 for x in re.findall('[\w]+',lines[i+2])] # re = regular expression and \w means all words that contain a-zA-Z0-9
            vertices=list(set(numbers)) # list of helices in the subgraph
            vertices.sort()
            cols = re.split('[\s]+',subgraph_lines[file_index]) #subgraph id
            file_index += 1
            
            subgraph_residues = list()
            for v in vertices:
                for h in range(0,len(Helices)): # search every helix in the list of Helices
                    if v == Helices[h].number:
                        subgraph_residues += Helices[h].Residues # add the list of residues from the helix to the list of residues of the subgraph

            for e in sub_edges:
                helices = [int(x)+1 for x in re.findall('[\w]+',e)]
                for f in range(0,len(Edges)):
                    if (helices[0] == Edges[f].helix1 and helices[1] == Edges[f].helix2) or (helices[0] == Edges[f].helix2 and helices[1] == Edges[f].helix1):
                        subgraph_residues += Edges[f].Residues # add the list of residues from the edge to the list of residues of the subgraph

            # creating the subgraph instance
            Subgraphs.append(Subgraph(vertices,subgraph_residues,cols[2]))


    # read the pdb file

    file = open(pdb_file,"r")
    pdbContent = file.readlines() # reading the pdb file
    file.close()

    prevRes = ""
    curRes = ""
    resNum = -1
    count = 0
    for i in range(0,len(pdbContent)):
        count += 1
        curRes = pdbContent[i][21:26] # read the chain id and the residue number of this line
        if curRes != prevRes: # new residue
            Bases.append(Base())
            resNum += 1
        Bases[resNum].pdb_lines.append(pdbContent[i]) # add the pdb line to the corresponding base
        prevRes = curRes

    #write the subgraph information into the file in the database
    write_subgraph = "/home/sj78/labwork/RAG-3Dual_Jan2019/Subgraphs/%s-subgraphs" %(sys.argv[1])
    file = open(write_subgraph,"w")
    for i in range(0,len(Subgraphs)):
        check_id = re.split('_',Subgraphs[i].graph_id)
        if int(check_id[0]) > 9: # if the subgraph has more than 9 vertices, do not print
            continue
        Subgraphs[i].printSubgraph(file)
    file.close()

    #read the bpseq file
    #file = open(bpseq_file,"r")
    #bpseqContent = file.readlines() # reading the pdb file
    #file.close()

    #writing the atomic fragments corresponding to the subgraphs in the database
    for i in range(0,len(Subgraphs)):
        
        check_id = re.split('_',Subgraphs[i].graph_id)
        if int(check_id[0]) > 9: # if the subgraph has more than 9 vertices
            continue

        if not os.path.isdir("/home/sj78/labwork/RAG-3Dual_Jan2019/Results/%s"%Subgraphs[i].graph_id): # create the subgraph id directpry if necessary
            os.mkdir("/home/sj78/labwork/RAG-3Dual_Jan2019/Results/%s"%Subgraphs[i].graph_id)
            os.mkdir("/home/sj78/labwork/RAG-3Dual_Jan2019/Results/%s/AllAtom"%Subgraphs[i].graph_id)
            #os.mkdir("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/Dual_RAG-3D/Results/%s/SecStruct"%Subgraphs[i].graph_id)

        #write the atomic fragment file
        atomic_frag_file = "/home/sj78/labwork/RAG-3Dual_Jan2019/Results/%s/AllAtom/%s-%s-%d-AA-frgt.pdb" %(Subgraphs[i].graph_id,sys.argv[1],Subgraphs[i].graph_id,i+1)
        file = open(atomic_frag_file,"w")
        for res in Subgraphs[i].Residues:
            Bases[res-1].printBase(file)
        file.close()

        #write the sec struct fragment file
        #bpseq_frag_file = "/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/Dual_RAG-3D/Results/%s/SecStruct/%s-%s-%d-AA-frgt.bpseq" %(Subgraphs[i].graph_id,sys.argv[1],Subgraphs[i].graph_id,i+1)
        #file = open(bpseq_frag_file,"w")
        #for res in Subgraphs[i].Residues:
        #   file.write(bpseqContent[res]) #bpseq base lines start from 1
        #file.close()

if __name__ == "__main__":
    main()
