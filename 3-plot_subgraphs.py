from igraph import *
from ClassesFunctions import *

#name="PDB_00746"
#ID="4_4"

#Read the matrix
#graph=sys.argv[1]
#pdb=sys.argv[2]
#number=sys.argv[3]
#g=Graph()
#g=Graph.Read_Adjacency("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/adj_matrices/%s_matrix.txt"%pdb, sep=" ", mode=ADJ_UNDIRECTED)
#g=Graph.Read_Adjacency("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/submatrices/%s/matrix%s.txt"%(pdb,number), sep="\t", mode=ADJ_UNDIRECTED)
#g=Graph.Read_Adjacency("/Users/cs4367/Desktop/Dual_Partitioning/rnastrand_results/adj_matrices/%s_matrix.txt"%pdb, sep=" ", mode=ADJ_UNDIRECTED)
#g=Graph.Read_Adjacency("/Users/sj78/Documents/labwork/RNAFilesfromPDB_Aug31/FinalBPSEQs/adj_matrices/4V8M_BE_matrix.txt", sep=" ", mode=ADJ_UNDIRECTED)
#g=Graph.Read_Adjacency("/Users/sj78/Documents/labwork/Papers_Writeups/DualGraphEnumeration/elsarticle-template/Figs_material/Incompatible_3_1.txt", sep="\t", mode=ADJ_UNDIRECTED)

#mylayout=g.layout_circle()
#bbox=BoundingBox(400,400)

#mylabels=[]
#for i in range(1,7):
#    mylabels.append(str(i))

#Red figures
#figure=plot(g, vertex_color="red", vertex_frame_color="black", vertex_size=40, edge_color="red", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))

#Black figures
#figure=plot(g, vertex_color="black", vertex_frame_color="black", vertex_size=50, edge_color="black", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(200,200,200,200))


#figure.save("./all_db/subgraphs/subgraph%s_%s.png"%(ID,name))
#figure.save("/home/sj78/labwork/Papers_Writeups/DualGraph_partitioning/elsarticle-template/Figs_material/%s.png"%graph)
#figure.save("/home/sj78/labwork/Papers_Writeups/DualGraph_partitioning/elsarticle-template/Figs_material/%s_sub.png"%graph)
#figure.save("/Users/sj78/Documents/labwork/Papers_Writeups/DualGraphEnumeration/elsarticle-template/Figs_material/Incompatible_6.png")


Graphs = []
num_vertices = 9
eigen_file = "9Eigen_map_sort"
adj_file = "V9AdjDG_map_sort"
loadEigenvalues(Graphs,num_vertices,eigen_file)
loadAdjMatrices(Graphs,num_vertices,adj_file)

for g1 in Graphs:

    g=Graph()
    g.add_vertices(num_vertices)
    #g = Graph.Adjacency(g1.adjMatrix,mode=ADJ_UNDIRECTED)
    #g=Graph.Read_Adjacency("/Users/sj78/Documents/labwork/DualGraphCodes/test_file", sep=" ", mode=ADJ_UNDIRECTED)
    #mylayout=g.layout_fruchterman_reingold()
    mylayout=g.layout_circle()
    bbox=BoundingBox(400,400)

    for i in range(0,num_vertices):
    
        for j in range(i,num_vertices):
    
            for k in range(0,int(g1.adjMatrix[i][j])):

                g.add_edge(i,j)

    print g.get_adjacency()

    if g1.graphID == "9_18405" or g1.graphID == "9_19203" or g1.graphID == "9_20559" or g1.graphID == "9_20569" or g1.graphID == "9_20790" or g1.graphID == "9_21508" or g1.graphID == "9_3051" or g1.graphID == "9_35189" or g1.graphID == "9_35458" or g1.graphID == "9_3854" or g1.graphID == "9_38596" or g1.graphID == "9_38597" or g1.graphID == "9_38598" or g1.graphID == "9_38599" or g1.graphID == "9_4495" or g1.graphID == "9_49214" or g1.graphID == "9_86359":
        figure=plot(g, "/Users/sj78/Documents/labwork/Papers_Writeups/DualGraphEnumeration/elsarticle-template/Figs_material/AllGraphImages/%s.png"%g1.graphID,vertex_color="red", vertex_frame_color="black", vertex_size=40, edge_color="red", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
    else:
    	figure=plot(g, "/Users/sj78/Documents/labwork/Papers_Writeups/DualGraphEnumeration/elsarticle-template/Figs_material/AllGraphImages/%s.png"%g1.graphID, vertex_color="black", vertex_frame_color="black", vertex_size=40, edge_color="black", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
    #figure.save("/Users/sj78/Documents/labwork/Papers_Writeups/DualGraphEnumeration/elsarticle-template/Figs_material/AllGraphImages/%s.png"%g1.graphID)

