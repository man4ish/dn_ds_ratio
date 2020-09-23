#!/usr/bin/python
import json
import csv
import re
import sys
import os
import pandas as pd

def get_codon():
    ''' readinfg codon json'''
    with open('json_data/codon.txt') as cf:
       codon_data = json.load(cf)
    return codon_data

def getmutation_table():
    '''reading mutation codon'''
    with open('json_data/mutation_codon.txt') as mcf:
       mutation_codon_data = json.load(mcf)
    return mutation_codon_data

def get_all_possible_path():
    '''reading all possible path'''
    with open('json_data/all_possible_path.txt') as appf:
         all_possible_path = json.load(appf)
    return all_possible_path

def get_triplets(seq, gene_id):
    codon_list = []
    mutation_codon_data = getmutation_table()
    for i in range(int(len(seq)/3)):
        codon = []
        start = 3*i
        end = 3*i+3
        codon.append(gene_id)
        triplet =seq[start:end]
        N = 0
        S = 0
        N1 = mutation_codon_data[triplet + "_N_1"]
        N2 = mutation_codon_data[triplet + "_N_2"]
        N3 = mutation_codon_data[triplet + "_N_3"]
        S1 = mutation_codon_data[triplet + "_S_1"]
        S2 = mutation_codon_data[triplet + "_S_2"]
        S3 = mutation_codon_data[triplet + "_S_3"]
        N = N + (N1 + N2 + N3)
        S = S + (S1 + S2 + S3)
        codon.append(seq[start:end])
        codon.append(start+1)
        codon.append(end)
        codon.append(str(3*i+1) + ", " + str(3*i+2) + ", " + str(3*i+3))
        codon.append(N)
        codon.append(S)
        codon_list.append(codon)

    return codon_list
    
def read_refseq(fasta_file):
    seq = ''
    with open(fasta_file) as fp:
        line = fp.readline()
        while line:
            if not line.startswith(">"):
               line = line.strip()
               seq = seq + line
            line = fp.readline()
    return(seq)

def gen_codonlist(seq, start, stop, gene_id):
    subseq = seq[start-1:stop]
    return get_triplets(subseq, gene_id)
    
def get_gff_file(gff_file):
    with open(gff_file) as fp:
        line = fp.readline()
        while line:
            if not line.startswith("#"):
               line = line.strip()
               rec = line.split("\t")
               cds_start = rec[3]
               cds_end   = rec[4]
               id = rec[8]
            gene_id = id.split("\"")[1]
            line = fp.readline()

def filter_ann(ann_field):
    if '>' in ann_field:
        return True
    else:
        return False

def read_vcf(vcf_file):
    varlist = []
    #['Chr01', '8', '.', 'C', 'G', '65.28', 'FS_filter;SOR_filter', '', 'GT:DP:AD', '0/1:300:100,200']
    with open(vcf_file) as fp:
        line = fp.readline()
        while line:
            if not line.startswith("#"):
               var = []
               coverage = {"A" : 0, "C": 0, "G" : 0, "T":0}
               line = line.strip()
               rec = line.split("\t")
               var.append(rec[0])
               var.append(rec[3])
               var.append(rec[4])
               var.append(rec[1])
               annotation = rec[7]
               print(annotation.split("|"))
               var_field = filter(filter_ann, annotation.split("|"))
               seq = read_refseq("sample.fa")
               pos_in_codon = 0
               for variation in var_field:
                   pos_in_codon = variation[2:(variation.find('>') - 1)]
                   print(pos_in_codon)
               var.append(pos_in_codon)  # get from snpeff results
               mod = int(pos_in_codon) % 3
               codon_start = (int(pos_in_codon) // 3) * 3 + 1
               print("codon_start=" + str(codon_start))
               print(seq)
               codon = seq[codon_start-1:codon_start + 2]
               var.append(codon)      #get from snpeff results


               var_class = annotation.split("|")[1]


               var.append(var_class)            #get from snpeff results
               format = (rec[8]).split(":")
               DP_index = format.index("DP")
               AD_index = format.index("AD")
               sample = (rec[9]).split(":")

               #print(annotation.split("|"))

               print("DP=" + sample[DP_index] + "\tAD=" +sample[AD_index])
               coverage[rec[3]] = (sample[AD_index]).split(",")[0]
               coverage[rec[4]] = (sample[AD_index]).split(",")[1]  #for single allele
               print(coverage) 
               var.append(coverage["A"])
               var.append(coverage["C"])
               var.append(coverage["G"])
               var.append(coverage["T"])
               varlist.append(var)

            line = fp.readline()
    return varlist      

get_gff_file("sample.gtf")
seq = read_refseq("sample.fa")

'''
data = read_vcf("snpeff_sample.ann.vcf")
print(data)

if os.path.exists("variant_info.tsv"):
  os.remove("variant_info.tsv")

with open('variant_info.tsv', 'a') as myfile:
    wr = csv.writer(myfile, delimiter='\t')
    for data_list in data:
        wr.writerow(data_list)
'''
print(gen_codonlist(seq, 1, 18, "gene_id1_cds1"))

#codonlist = get_triplets(seq)


'''
Calculating Nd and Sd for input vcf file
with open("input.vcf") as fvp:
    vcfline = fvp.readline()
    while vcfline:
       if not vcfline.startswith("#"):
          vcfline = vcfline.strip()
       fine = fvp.readline()
'''
