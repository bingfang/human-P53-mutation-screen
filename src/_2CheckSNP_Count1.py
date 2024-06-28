#!/usr/local/bin/python3.6
from Bio.Seq import Seq
import glob
def main():
    
    list1 = ['AGGCTGCTCCCCCC','AGGCTGCCCCCCCC'] # exon4_CCCCC
    list2 = ['ctacaagcagtcac'.upper(),'CTACAAGTAGTCAC'] # exon5_Q165ter
    list3 = ['TGTACCACCATCC','TGTACCCCCATCCA'] # exon7_5end
    list4 = ['gaaccggaggccc'.upper(),'GAACCAGAGGCC'] # exon7_CCCCCR248Q    
    list5 = ['AGTACTCCCCTG','AGTACTCCCCTG'] # exon5_5end
    list6 = ['TGTACCACCATCC','TGTACCCCCATCCA'] # exon7_5end

    dic={"ex4_CCCC":list1, "ex5_Q165ter":list2, "ex7_5end":list3, "ex7_Q248Q":list4, "ex5_5end":list5, "ex7_5end":list6}
    for i in dic:
        print(dic[i])
        outputfilename = str(i) + "_snp_check1.txt" 
        for name in glob.glob('./1*2*.txt'):
            inputfilename = str(name)
            with open(inputfilename, 'r') as f:
                data_in = f.read().rstrip().split('\n')
            WT = Seq(dic[i][0])
            WT_RC = WT.reverse_complement()
            mut= Seq(dic[i][1])
            mut_RC = mut.reverse_complement() 
            wt_F, wt_R, mut_F, mut_R,F_rate,R_rate = check_SNP(data_in,str(WT),str(WT_RC),str(mut),str(mut_RC))
            print(F_rate,R_rate )
            with open(outputfilename, "a") as f:
                f.write(inputfilename[2:]+'\n')
                f.write("{}{}{}{}\n".format("WT_Forward\t",str(WT),"\t",wt_F))
                f.write("{}{}{}{}\n".format("WT_Reward\t",str(WT_RC),"\t",wt_R))
                f.write("{}{}{}{}{}{:.2%}\n".format("mut_Forward\t",str(mut),"\t",mut_F,"\t",F_rate))
                f.write("{}{}{}{}{}{:.2%}\n".format("mut_Reward\t",str(mut_RC),"\t",mut_R,"\t",R_rate))
        

    
def check_SNP(data_in,str1,str2,str3,str4):
    wt_F = 0
    mut_F = 0
    wt_R = 0
    mut_R = 0
    for line in data_in:
        if str1 in line:
            wt_F += 1
        if str2 in line:
            wt_R += 1
        if str3 in line:
            mut_F += 1
        if str4 in line:
            mut_R += 1
    F_rate= float(mut_F)/((float(wt_F)+float(mut_F))+1)
    R_rate= float(mut_R)/((float(wt_R)+float(mut_R))+1)     
    return wt_F, wt_R, mut_F, mut_R, F_rate,R_rate

main()

