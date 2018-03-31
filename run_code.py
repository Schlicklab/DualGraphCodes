#Need to change file paths in dualGraphs.py

import os,sys
#name=sys.argv[1]
path="/Users/cs4367/Desktop/Dual_Partitioning/ct_files_from_ribovision/with_pknot/"
#names=os.listdir(path)

#names=["DM_LSU.bpseq", "EC_LSU.bpseq",	"HS_LSU.bpseq",	"PF_LSU.bpseq", "SC_LSU.bpseq"]
names=["DM_LSU.bpseq"]
#os.system("python dualGraphs.py /Users/cs4367/Documents/BPSEQ_2016/%s.bpseq"%name)
for name in names:
    print name
    os.system("python dualGraphs.py %s/%s > /Users/cs4367/Desktop/Dual_Partitioning/rRNA_ribovision/output/%s_basepairs.txt"%(path,name,name.split(".")[0]))
    n=open("/Users/cs4367/Desktop/Dual_Partitioning/rRNA_ribovision/adj_matrices/n.txt","r").readline()
    os.system("./dualgraph.out -input /Users/cs4367/Desktop/Dual_Partitioning/rRNA_ribovision/adj_matrices/%s_matrix.txt -len %s -output /Users/cs4367/Desktop/Dual_Partitioning/rRNA_ribovision/output/%s_output.txt"%(name.split(".")[0],n,name.split(".")[0]))
