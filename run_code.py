import os,sys

pdbfiles=[]
vertices=[]
#n=open("/home/sj78/labwork/LWilliams_work/Test_xrna/adj_matrices/n.txt","r").readline()
#n="38"
#os.system("./dualgraph.out -input /home/sj78/labwork/LWilliams_work/Test_xrna/adj_matrices/%s_matrix.txt -len %s -output /home/sj78/labwork/LWilliams_work/Test_xrna/output/%s_output.txt"%(name.split(".")[0],n,name.split(".")[0]))
#path="/home/sj78/labwork/AllRNADataset_August2018"
path="/Users/sj78/Documents/labwork/RNAFilesfromPDB_Aug31/FinalBPSEQs"

mylines=[x.strip() for x in open(sys.argv[1],"r").readlines()]
for line in mylines:

    cols = [x for x in line.split()]
    numbers = [int(x) for x in cols[1].split('_')]
    if numbers[0] != 1:
        pdbfiles.append(cols[0])
        vertices.append(numbers[0])

for i in range(0,len(pdbfiles)):
    print pdbfiles[i],
    print vertices[i]
    os.system("./dualgraph.out -input %s/adj_matrices/%s_matrix.txt -len %d -output %s/output/%s_output.txt -all %s/output_subgraphs/%s_output.txt"%(path,pdbfiles[i],vertices[i],path,pdbfiles[i],path,pdbfiles[i]))


#myfiles=[x.strip() for x in open("/home/sj78/labwork/AllRNADataset_August2018/temp_list","r").readlines()]
#myfiles=[x.strip() for x in open("/Users/sj78/Documents/labwork/RNAFilesfromPDB_Aug31/FinalBPSEQs/doAgain_list","r").readlines()]
#for name in myfiles:
#    print name
#    #os.system("python dualGraphs.py %s/%s.bpseq > /home/sj78/labwork/AllRNADataset_August2018/output/%s_basepairs.txt"%(path,name,name))
#    os.system("python dualGraphs.py %s/%s.bpseq > %s/output/%s_basepairs.txt"%(path,name,path,name))
