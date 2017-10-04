#!/usr/bin/python
from optparse import OptionParser
import math
def infer_anc():
    # parse command
    parser = OptionParser()
    parser.add_option("-p", "--par", dest="par", metavar="FILE", help="input parameter file path")
    (options, args) = parser.parse_args()
    input_par = options.par

    # input paramenter and files
    input_file = open(input_par, "r")
    input = input_file.readlines()
    for line in input:
        if line.split(":")[0]=="geno":
            geno_file = open(line.split()[1],"r")
        if line.split(":")[0]=="snp":
            snp_file = open(line.split()[1],"r")
        if line.split(":")[0]=="ind":
            ind_file = open(line.split()[1],"r")
        if line.split(":")[0]=="snpwt":
            snpwt_file = open(line.split()[1],"r")
        if line.split(":")[0]=="predpcoutput":
            output_file = open(line.split()[1],"w")
	
    # par_file = open("./inferancestry.info","r")
    # wt_file = open("./snpwt."+pop,"r")
    geno = geno_file.readlines()
    snp  = snp_file.readlines()
    ind  = ind_file.readlines()
    weight = snpwt_file.readlines()
    
    # assign number of PC to predict
    npc = len(weight[0].rstrip().split())
    npop = npc+1
    shrinkage = [float(x) for x in weight[0].rstrip().split()]
    coef_pop = [float(x) for x in weight[4].rstrip().split()]
    wt_number = float(len(weight)-5)

    # for strand check
    TRANS = { "T": "A", "A": "T", "G": "C", "C": "G" }        
    
    # create snpwt dictionary
    snpwt = {}
    for row in range(5,len(weight)):
        snpwt_line = weight[row].rstrip().split(' ')
        snpwt[snpwt_line[0]] = snpwt_line[1:(npc+4)] 

    # calculate predicted PCs
    for col in range(len(geno[0].rstrip())):
        print >>output_file, ind[col].split()[0], ind[col].split()[2], 
        pc_adj_list = []
        for pc in range(npc):
            WX = 0.0
            n = 0.0
            for row in range(len(geno)):
                snp_input = snp[row].rstrip().split()
                snpwt_tmp = snpwt.get(snp_input[0], "NA")
                if snpwt_tmp != "NA":
                    X = float(geno[row].rstrip()[col])
                    if X == 9.0: continue
                    elif (snpwt_tmp[0] == snp_input[4] and snpwt_tmp[1] == snp_input[5]) or (snpwt_tmp[0] == TRANS[snp_input[4]] and snpwt_tmp[1] == TRANS[snp_input[5]]):
                        P = float(snpwt_tmp[2])
                        W = float(snpwt_tmp[(pc+3)])
                        WX = WX + W*((X-P*2.0)/(math.sqrt((P)*(1.0-(P)))))
                        n += 1.0
                    elif (snpwt_tmp[0] == snp_input[5] and snpwt_tmp[1] == snp_input[4]) or (snpwt_tmp[0] == TRANS[snp_input[5]] and snpwt_tmp[1] == TRANS[snp_input[4]]): 
                        X = 2.0-X
                        P = float(snpwt_tmp[2])
                        W = float(snpwt_tmp[(pc+3)])
                        WX = WX + W*((X-P*2.0)/(math.sqrt((P)*(1.0-(P)))))
                        n += 1.0
                    else: continue
            pred_pc_adj = WX*(wt_number/n)/shrinkage[pc]
            pc_adj_list.append(pred_pc_adj)

        percent_adj_list = [0.0]*npop
        for i in range(npop):
            for pc in range(npc):
                percent_adj_list[i] += coef_pop[(pc+i*npop)]*pc_adj_list[pc]
            percent_adj_list[i] += coef_pop[(npc+i*npop)]
        # percent_adj_list.append(1.0-sum(percent_adj_list))
        
        for i in range(len(percent_adj_list)):
            if percent_adj_list[i] < 0.0:
                percent_adj_list[i] = 0.0
            if percent_adj_list[i] > 100.0:
                percent_adj_list[i] = 100.0
        if sum(percent_adj_list) > 1.0:
            sum_denom = sum(percent_adj_list)
            for i in range(len(percent_adj_list)):
                percent_adj_list[i] = percent_adj_list[i]/sum_denom
        
        result_list = pc_adj_list+percent_adj_list 
        print >>output_file, int(n), 
        for i in result_list:
            print >>output_file, i,
        print >>output_file, " "
            
    geno_file.close()
    snp_file.close()
    ind_file.close()
    snpwt_file.close()
    output_file.close()
        
if __name__ == "__main__":
    infer_anc()
