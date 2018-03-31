from igraph import *

name="PDB_00746"
ID="4_4"

#Read the matrix
g=Graph()
#g=Graph.Read_Adjacency("/Users/cs4367/Desktop/Dual_Partitioning/rnastrand_results/adj_matrices/%s_matrix.txt"%name, sep=" ", mode=ADJ_UNDIRECTED)
g=Graph.Read_Adjacency("C:/Users/cigdems/Desktop/matrix.txt", sep=" ", mode=ADJ_UNDIRECTED)
mylayout=g.layout_circle()
bbox=BoundingBox(800,800)

#Red figures
figure=plot(g, vertex_color="red", vertex_frame_color="red", vertex_size=50, edge_color="red", edge_width=2, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(200,200,200,200))

#Black figures
#figure=plot(g, vertex_color="black", vertex_frame_color="black", vertex_size=50, edge_color="black", edge_width=2, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(200,200,200,200))


#figure.save("./all_db/subgraphs/subgraph%s_%s.png"%(ID,name))
figure.save("C:/Users/cigdems/Desktop/hairpin_dual.png")
