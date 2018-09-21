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
#g=Graph.Read_Adjacency("/Users/sj78/Documents/GrantsChapters/FiguresForRNAGraphsChapter/Figs_material/1s72_dualmat.txt", sep=" ", mode=ADJ_UNDIRECTED)
#g=Graph.Read_Adjacency("/Users/sj78/Documents/labwork/DualGraphCodes/test_file", sep=" ", mode=ADJ_UNDIRECTED)

#mylayout=g.layout_circle()
#bbox=BoundingBox(800,800)

#Red figures
#figure=plot(g, vertex_color="red", vertex_frame_color="black", vertex_size=50, edge_color="red", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(200,200,200,200))

#Black figures
#figure=plot(g, vertex_color="black", vertex_frame_color="black", vertex_size=50, edge_color="black", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(200,200,200,200))


#figure.save("./all_db/subgraphs/subgraph%s_%s.png"%(ID,name))
#figure.save("/home/sj78/labwork/Papers_Writeups/DualGraph_partitioning/elsarticle-template/Figs_material/%s.png"%graph)
#figure.save("/home/sj78/labwork/Papers_Writeups/DualGraph_partitioning/elsarticle-template/Figs_material/%s_sub.png"%graph)
#figure.save("/Users/sj78/Documents/GrantsChapters/FiguresForRNAGraphsChapter/Figs_material/1s72_2Ddual.png")


Graphs = []
num_vertices = 4
eigen_file = "4Eigen_map_sort"
adj_file = "V4AdjDG_map_sort"

loadEigenvalues(Graphs,num_vertices,eigen_file)
loadAdjMatrices(Graphs,num_vertices,adj_file)


for g1 in Graphs:

    if g1.graphID != "4_13":
	continue

    g=Graph()
    g.add_vertices(num_vertices)
    #g = Graph.Adjacency(g1.adjMatrix,mode=ADJ_UNDIRECTED)
    #g=Graph.Read_Adjacency("/Users/sj78/Documents/labwork/DualGraphCodes/test_file", sep=" ", mode=ADJ_UNDIRECTED)
    #mylayout=g.layout_fruchterman_reingold()
    mylayout=g.layout_circle()
    bbox=BoundingBox(800,800)

    for i in range(0,num_vertices):
    
        for j in range(i,num_vertices):
    
            for k in range(0,int(g1.adjMatrix[i][j])):

                g.add_edge(i,j)

    print g.get_adjacency()

    figure=plot(g, vertex_color="red", vertex_frame_color="black", vertex_size=50, edge_color="red", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(200,200,200,200))
    figure.save("/Users/sj78/Documents/labwork/Papers_Writeups/DualGraphEnumeration/Figs_material/%s.png"%g1.graphID)
