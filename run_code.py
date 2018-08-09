import os,sys

path="/home/sj78/labwork/NonRedundantList_Oct2017/BPSEQs/"

myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/list_bp_ext","r").readlines()]
for name in myfiles:
    print name
    os.system("python dualGraphs.py %s/%s > /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/output/%s_basepairs.txt"%(path,name,name.split(".")[0]))
    #n=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/adj_matrices/n.txt","r").readline()
    #os.system("./dualgraph.out -input /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/adj_matrices/%s_matrix.txt -len %s -output /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/output/%s_output.txt -all /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/output/%s_Subgraphs.txt"%(name.split(".")[0],n,name.split(".")[0],name.split(".")[0]))
