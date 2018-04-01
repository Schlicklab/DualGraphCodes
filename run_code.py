#Need to change file paths in dualGraphs.py
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/all_structures.txt","r").readlines()]
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/all_structures.txt","r").readlines()]

#for filename in myfiles:

import os,sys
#name=sys.argv[1]
#path="/home/sj78/labwork/NonRedundantList_Oct2017/BPSEQs/"
path="/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/ct_files_from_ribovision/with_pknot/"
#path="/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/"
#names=os.listdir(path)

#names=["DM_SSU.bpseq", "EC_SSU.bpseq",	"HS_SSU.bpseq",	"PF_SSU.bpseq", "SC_SSU.bpseq", "DM_LSU.bpseq", "EC_LSU.bpseq",  "HS_LSU.bpseq", "PF_LSU.bpseq", "SC_LSU.bpseq"]
#names=["DM_LSU.bpseq"]
#names=["PDB_00573.ct","PDB_00573_1mol.ct"]
#os.system("python dualGraphs.py /Users/cs4367/Documents/BPSEQ_2016/%s.bpseq"%name)
#myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/list_bp_ext","r").readlines()]
myfiles=[x.strip() for x in open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/all_structures.txt","r").readlines()]
#for name in names:
for name in myfiles:
    print name
    #os.system("python dualGraphs.py %s/%s > /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/output/%s_basepairs.txt"%(path,name,name.split(".")[0]))
    os.system("python dualGraphs.py %s/%s > /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/output/%s_basepairs.txt"%(path,name,name.split(".")[0]))
    #os.system("python dualGraphs.py %s%s > %s%s_basepairs.txt"%(path,name,path,name.split(".")[0]))
    #n=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/adj_matrices/n.txt","r").readline()
    n=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/adj_matrices/n.txt","r").readline()
    #n=open("/home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/CODES/Test/n.txt","r").readline()
    #os.system("./dualgraph.out -input /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/adj_matrices/%s_matrix.txt -len %s -output /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/NonRed2017_results/output/%s_output.txt"%(name.split(".")[0],n,name.split(".")[0]))
    os.system("./dualgraph.out -input /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/adj_matrices/%s_matrix.txt -len %s -output /home/sj78/labwork/DualGraphs_Cigdem/Dual_Partitioning/rRNA_ribovision/output/%s_output.txt"%(name.split(".")[0],n,name.split(".")[0]))
    #os.system("./dualgraph.out -input %s%s_matrix.txt -len %s -output %s%s_output.txt"%(path,name.split(".")[0],n,path,name.split(".")[0]))
