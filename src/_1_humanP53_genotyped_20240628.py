#!/usr/local/bin/python3.6

from Bio.Seq import Seq
import glob
from operator import itemgetter

# check format of input files
# check if P53_Miseq_screen.txt is pre-exist
# check the filter for mutation rate

# exon sequences with splicing sites.
exon2_F = "agcagccagactgccttccgggtcactgccatggaggagccgcagtcagatcctagcgtcgagccccctctgagtcaggaaacattttcagacctatggaaactgt".upper()
exon2_R = Seq(exon2_F).reverse_complement()
exon2_R = str(exon2_R)

exon3_F = "agacttcctgaaaacaacgttctggt".upper()
exon3_R = Seq(exon3_F).reverse_complement()
exon3_R = str(exon3_R)

exon4_F = "agtcccccttgccgtcccaagcaatggatgatttgatgctgtccccggacgatattgaacaatggttcactgaagacccaggtccagatgaagctcccagaatgccagaggctgctccccccgtggcccctgcaccagcagctcctacaccggcggcccctgcaccagccccctcctggcccctgtcatcttctgtcccttcccagaaaacctaccagggcagctacggtttccgtctgggcttcttgcattctgggacagccaagtctgtgacttgcacggt".upper()
exon4_R = Seq(exon4_F).reverse_complement()                          
exon4_R = str(exon4_R)

exon5_F ="agtactcccctgccctcaacaagatgttttgccaactggccaagacctgccctgtgcagctgtgggttgattccacacccccgcccggcacccgcgtccgcgccatggccatctacaagcagtcacagcacatgacggaggttgtgaggcgctgcccccaccatgagcgctgctcagatagcgatggt".upper()
exon5_R = Seq(exon5_F).reverse_complement()
exon5_R = str(exon5_R)

exon6_F = "aggtctggcccctcctcagcatcttatccgagtggaaggaaatttgcgtgtggagtatttggatgacagaaacacttttcgacatagtgtggtggtgccctatgagccgcctgaggt".upper()
exon6_R = Seq(exon6_F).reverse_complement()
exon6_R = str(exon6_R)

exon7_F = "aggttggctctgactgtaccaccatccactacaactacatgtgtaacagttcctgcatgggcggcatgaaccggaggcccatcctcaccatcatcacactggaagactccaggt".upper()
exon7_R = Seq(exon7_F).reverse_complement()
exon7_R = str(exon7_R)

exon8_F = "agtggtaatctactgggacggaacagctttgaggtgcgtgtttgtgcctgtcctgggagagaccggcgcacagaggaagagaatctccgcaagaaaggggagcctcaccacgagctgcccccagggagcactaagcgaggt".upper()
exon8_R = Seq(exon8_F).reverse_complement()
exon8_R = str(exon8_R)

exon9_F =  "agcactgcccaacaacaccagctcctctccccagccaaagaagaaaccactggatggagaatatttcacccttcaggt".upper()
exon9_R = Seq(exon9_F).reverse_complement()
exon9_R = str(exon9_R)

exon10_F = "aggaccagaccagctttcaaaaagaaaattgttaaagagagcatgaaaatggttctatgactttgcctgatacagatgctacttgacttacgatggtgttacttcctgataaactcgtcgtaagttgaaaatattgt".upper()
exon10_R = Seq(exon10_F).reverse_complement()
exon10_R = str(exon10_R)

exon11_F = "agatccgtgggcgtgagcgcttcgagatgttccgagagctgaatgaggccttggaactcaaggatgcccaggctgggaaggagccaggggggagcagggctcactccaggt".upper()
exon11_R = Seq(exon11_F).reverse_complement()
exon11_R = str(exon11_R)

exon12_F = "agccacctgaagtccaaaaagggtcagtctacctcccgccataaaaaactcatgttcaagacagaagggcctgactcagactgac".upper()
exon12_R = Seq(exon12_F).reverse_complement()
exon12_R = str(exon12_R)


def main():
    
    for name in glob.glob('../data/1*2*.txt'):                            
        inputfilename = str(name)
        print(inputfilename)
        outputfilename1 = "../data/human_P53_Miseq_screen.txt"
        outputfilename2 = "../data/Top100_seq" + str(name[20:-4])+".txt"
        outputfilename3 = "../data/Check_seq_all_human.txt"
        
        # input one sample each time
        with open(inputfilename, 'r') as f:                         
            data_in = f.read().rstrip().split('\n')
        
        # at the beginning, counted_line have 3 columns [seq, count_read, Percentage in total reads]                                      
        counted_line = get_unique(data_in)
        
        # start a new dictionary.
        exon_reads = {"Exon2_3":0, "Exon4":0, "Exon5":0, "Exon6":0, "Exon7":0, "Exon8":0, "Exon9":0, "Exon10":0, "Exon11":0, "Exon12":0,"unknown":0} 
        
        # renewed counted_line have 4 columns [seq, count_read, Percentage in total reads,Exon]
        counted_line, exon_reads  = count_exon_reads(counted_line,exon_reads)
        print(len(counted_line))
        
        # renewed counted_line have 5 columns [seq, count_read, Percentage in total reads, Exon, Percentage in each exon]
        counted_line = filter_reads(counted_line, exon_reads)
        print(len(counted_line))
        
        # output a summary file all samples including top 100 reads for each sample
        with open(outputfilename1, 'a') as f:
            f.write("{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}\n".format("Sample ID","\t","Order","\t","Sequence","\t","Reads","\t","Percentage in total reads","\t","Exon","\t","Percentage in each exon","\t","Genotype"))
            for i in range(0,100):
                f.write("{}{}{}{}{}{}{:d}{}{:.2%}{}{}{}{:.2%}{}{}\n".format(inputfilename[2:],"\t",str(i+1),"\t",counted_line[i][0],"\t",counted_line[i][1],"\t",counted_line[i][2],"\t",counted_line[i][3],"\t",counted_line[i][4],"\t",counted_line[i][5]))  
        # output top 50 reads for each sample in fasta format
        with open(outputfilename2, 'w') as f:
            for i in range(0,100):
                f.write(">"+inputfilename[20:-4]+"_"+str(i+1)+"\n") 
                f.write(counted_line[i][0]+"\n")
        # output filtered reads (>1%) for each sample in fasta format. The actual mutation rate is about 4% considering the total reads including read1,read2, and many bad reads.
        with open(outputfilename3, 'a') as f:
            for i in range(0,300):
                if  counted_line[i][4] > 0.01 and counted_line[i][5] == "Check":
                    f.write(">"+inputfilename[20:-4]+"_"+str(i+1)+"\n")
                    f.write(counted_line[i][0]+"\n")	

# abstract read1 and read2 and put into seq_read12 list.
# create a set with only unique reads.
# count each unique read. 
# calculate percentage of each unique read in total reads             
def get_unique(data_in):
    seq_read12 =[] 
    for i in range(len(data_in)):  
        barcoderead1read2=data_in[i].split('\t')
        seq_read12.append(barcoderead1read2[1])
        seq_read12.append(barcoderead1read2[2])
    total_read = len(seq_read12)
    unique_seq = set(seq_read12)                            ### create a set with only unique read.              
    counted_line = []
    for seq in unique_seq: 
        if len(seq) >100:                                   ### remove primer dimers
            count_read = seq_read12.count(seq)              ### count each unique read.
            percent = float(count_read)/float(total_read)
            counted_line.append([seq, count_read, percent]) ### counted_line is a list of list.
    return counted_line
    
# mark exon number for each Read
# calculate total reads for each exon
# sort reads by count_read
       
def count_exon_reads(counted_line,exon_reads):       
    for line in counted_line:
        if "GGGTTGGAAGTGTCTCATGCTG" in line[0] or "AGCAGTCAGAGGACCAGGTC" in line[0]:
            line.append("Exon2_3")
            exon_reads["Exon2_3"] += float(line[1])
        elif "GACTGCTCTTTTCACCCATCTAC" in line[0] or "CGGCCAGGCATTGAAGTC" in line[0]:
            line.append("Exon4")
            exon_reads["Exon4"] += float(line[1])              
        elif "CTGACTTTCAACTCTGTCTCCTTCC" in line[0] or "agcaatcagtgaggaatcagagg".upper() in line[0]:
            line.append("Exon5")
            exon_reads["Exon5"] += float(line[1])
        elif "cctctgattcctcactgattgct".upper() in line[0] or "AGCCCTGTCGTCTCTCC" in line[0]:
            line.append("Exon6")
            exon_reads["Exon6"] += float(line[1])
        elif "CCTCATCTTGGGCCTGTGTTATC" in line[0] or "TGATGAGAGGTGGATGGGTA" in line[0]:
            line.append("Exon7")
            exon_reads["Exon7"] += float(line[1])
        elif "CTTAGGCTCCAGAAAGGACAAGG" in line[0] or "TCTCCTCCACCGCTTCTTG" in line[0]:
            line.append("Exon8")
            exon_reads["Exon8"] += float(line[1])
        elif "CCTCAGATTCACTTTTATCACCTTTCCT" in line[0] or "TTAGTTAGCTACAACCAGGAGCCA" in line[0]:
            line.append("Exon9")   
            exon_reads["Exon9"] += float(line[1])    
        elif "GCTAACTAACTTCAGAACACCAACTTATACC" in line[0] or "AGCAGGCTAGGCTAAGCTATG" in line[0]:
            line.append("Exon10")
            exon_reads["Exon10"] += float(line[1])
        elif "TGAACCATCTTTTAACTCAGGTACTGTG" in line[0] or "TGAAGGCAGGATGAGAATGGAATCCTATG" in line[0]:
            line.append("Exon11")
            exon_reads["Exon11"] += float(line[1])
        elif "CACTCATGTGATGTCATCTCTCCTCC" in line[0] or "CTTCTGACGCACACCTATTGC" in line[0]:
            line.append("Exon12")
            exon_reads["Exon12"] += float(line[1])
        else:
            line.append("unknown")
            exon_reads["unknown"] += float(line[1])
    print("unique reads:", len(counted_line))
    print(exon_reads)    
    return counted_line, exon_reads

# calculate percentage of each read in total reads of each exon 
def filter_reads(counted_line, exon_reads):
    for line in counted_line:
        if line[3] == "Exon2_3":
            percent_in_exon = float(line[1])/exon_reads["Exon2_3"]
            line.append(percent_in_exon)
            if exon2_F in line[0] or exon2_R in line[0]:
                line.append("WT_exon2")
            elif exon3_F in line[0] or exon3_R in line[0]:
                line.append("WT_exon3")
            else:
                line.append("Check")
        elif line[3] == "Exon4":
            percent_in_exon = float(line[1])/exon_reads["Exon4"]
            line.append(percent_in_exon)
            if  exon4_F in line[0] or exon4_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")       
        elif line[3] == "Exon5":
            percent_in_exon = float(line[1])/exon_reads["Exon5"]
            line.append(percent_in_exon)
            if exon5_F in line[0] or exon5_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")                
        elif line[3] == "Exon6":
            percent_in_exon = float(line[1])/exon_reads["Exon6"]
            line.append(percent_in_exon)
            if exon6_F in line[0] or exon6_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon7":
            percent_in_exon = float(line[1])/exon_reads["Exon7"]
            line.append(percent_in_exon)
            if exon7_F in line[0] or exon7_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon8":
            percent_in_exon = float(line[1])/exon_reads["Exon8"]
            line.append(percent_in_exon)
            if exon8_F in line[0] or exon8_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon9":
            percent_in_exon = float(line[1])/exon_reads["Exon9"]
            line.append(percent_in_exon)
            if exon9_F in line[0] or exon9_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon10":
            percent_in_exon = float(line[1])/exon_reads["Exon10"]
            line.append(percent_in_exon)
            if exon10_F in line[0] or exon10_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon11":
            percent_in_exon = float(line[1])/exon_reads["Exon11"]
            line.append(percent_in_exon)
            if exon11_F in line[0] or exon11_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")
        elif line[3] == "Exon12":
            percent_in_exon = float(line[1])/exon_reads["Exon12"]
            line.append(percent_in_exon)
            if exon12_F in line[0] or exon12_R in line[0]:
                line.append("WT")
            else:
                line.append("Check")     
        elif line[3] == "unknown":
            percent_in_exon = float(line[1])/exon_reads["unknown"]
            line.append(percent_in_exon)
            line.append("check")
            print(line)
        else:
            print("exon is unknown")
        print(len(line))
    ### sort unique sequence by count_read
    counted_line.sort(key=itemgetter(1), reverse=True)   # cannot assign counted_line.sort(key=itemgetter(1), reverse=True) to a new variable    
    return counted_line
    
main()    
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	