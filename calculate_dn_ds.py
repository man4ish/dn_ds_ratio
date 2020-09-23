#!/usr/bin/python
import json
import csv
import os

def get_codon():
    ''' readinfg codon json'''
    with open('json_data/codon.json') as cf:
       codon_data = json.load(cf)
    return codon_data

def getmutation_table():
    '''reading mutation codon'''
    with open('json_data/mutation_codon.json') as mcf:
       mutation_codon_data = json.load(mcf)
    return mutation_codon_data

def get_all_possible_path():
    '''reading all possible path'''
    with open('json_data/all_possible_path.json') as appf:
         all_possible_path = json.load(appf)
    return all_possible_path

def get_all_possible_combination(n):
    return n*(n-1)/2

def get_coverage_product(num_A, num_C, num_G, num_T):
    num_AC = num_A*num_C
    num_AG = num_A*num_G
    num_AT = num_A*num_T
    num_CG = num_C*num_G
    num_CT = num_C*num_T
    num_GT = num_G*num_T
    coverage_list = [num_AC, num_AG, num_AT, num_CG, num_CT, num_GT]
    return coverage_list

def possible_codon(key, alleles):
    '''generate a list of possible codon with given position and allele list'''

    codon = {}
    ref  = key.split("-")[2]
    for i in range(3):
        codon[i+1] = [ref[i]]

    for allele in alleles:
        for pos in allele:
            alt_allele_list = allele[pos]
            for alt_allele in alt_allele_list:
                codon[pos].append(alt_allele)

    possible_cdn_list = []
    for i in range(len(codon[1])):
        for j in range(len(codon[2])):
            for k in range(len(codon[3])):
                possible_cdn_list.append(codon[1][i] + codon[2][j] + codon[3][k])
    return possible_cdn_list


def get_all_possible_codon(codon_list):
    '''get all possible codon with given allele'''
    print(codon_list)
    cdn_dict = {}
    for cdn_lst in codon_list:
        key = str(cdn_lst[0]) + "-" + str(cdn_lst[6]) + "-" + cdn_lst[7]
        if not key in cdn_dict:

           cdn_dict[key] = []
           pos_dict = {}
           pos_in_codon =cdn_lst[5]
           pos_dict[pos_in_codon] = [cdn_lst[2]]
           cdn_dict[key].append(pos_dict)
           #allele_frq_dict = {}
           #allele_frq_dict[pos_in_codon] = cdn_lst[9]
           #cdn_dict[key].append(allele_frq_dict)
        else:
            pos_dict = {}
            pos_in_codon = cdn_lst[5]
            pos_dict[pos_in_codon] = [cdn_lst[2]]
            #allele_frq_dict = {}
            #allele_frq_dict[pos_in_codon] = cdn_lst[9]
            cdn_dict[key].append(pos_dict)
            #cdn_dict[key].append(allele_frq_dict)

    for key in cdn_dict:
        all_codon = possible_codon(key, cdn_dict[key])
        print(all_codon)

    return cdn_dict


def get_triplets(seq, gene_id):
    '''generate triplets from ref seq'''

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
    '''read fasta file'''

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
    '''generate codon list for gene'''

    subseq = seq[start-1:stop]
    return get_triplets(subseq, gene_id)
    
def read_gff_file(gff_file):
    '''read gtf file'''

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
    '''parse vcf file'''

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
               if(mod == 0):
                   mod = 3
               codon_start = (int(pos_in_codon) // 3) * 3 + 1
               print("codon_start=" + str(codon_start))
               print(seq)
               codon = seq[codon_start-1:codon_start + 2]
               var.append(mod)
               var.append(codon_start)
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
               var.append(coverage)
               #var.append(coverage["A"])
               #var.append(coverage["C"])
               #var.append(coverage["G"])
               #var.append(coverage["T"])
               varlist.append(var)

            line = fp.readline()
    return varlist      

gff_file  = read_gff_file("sample.gtf")
seq = read_refseq("sample.fa")


result_list = gen_codonlist(seq, 1, 18, "gene_id1_cds1")
if os.path.exists("codon_results.tsv"):
  os.remove("codon_results.tsv")
with open('codon_results.tsv', 'a') as cdr_file:
    cdr = csv.writer(cdr_file, delimiter='\t')
    for cd_list in result_list:
        cdr.writerow(cd_list)

data = read_vcf("snpeff_sample.ann.vcf")

get_all_possible_codon(data)
exit(data)
if os.path.exists("variant_info.tsv"):
  os.remove("variant_info.tsv")
with open('variant_info.tsv', 'a') as myfile:
    wr = csv.writer(myfile, delimiter='\t')
    for data_list in data:
        wr.writerow(data_list)



'''
Calculating Nd and Sd for input vcf file
with open("input.vcf") as fvp:
    vcfline = fvp.readline()
    while vcfline:
       if not vcfline.startswith("#"):
          vcfline = vcfline.strip()
       fine = fvp.readline()
'''
