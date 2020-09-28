import os
import json
class parse_dataUtils:

    def getmutation_table(self):
        '''reading mutation codon'''

        with open('json_data/mutation_codon.json') as mcf:
            mutation_codon_data = json.load(mcf)
        return mutation_codon_data

    def get_codon(self, seq, gene_id):
        codon_list = []
        mutation_codon_data = self.getmutation_table()

        for i in range(int(len(seq) / 3)):
            codon = []
            start = 3 * i
            end = 3 * i + 3
            codon.append(gene_id)
            triplet = seq[start:end]
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
            #codon.append((cds_start - 1) + start + 1)
            #codon.append((cds_start - 1) + end)
            #codon.append(str((cds_start - 1) + 3 * i + 1) + ", " + str((cds_start - 1) + 3 * i + 2) + ", " + str(
            #    (cds_start - 1) + 3 * i + 3))
            codon.append(N)
            codon.append(S)
            codon_list.append(codon)

        return codon_list

    def get_triplets(self, seq, chr_dict):
        '''generate triplets from ref seq'''

        for chr in chr_dict.keys():
            for genes in chr_dict[chr]:
                trascipt_seq = ''
                for gene_id, cds_list in genes.items():
                    for cds_coordinates in cds_list:
                        min = int(cds_coordinates[0])
                        max = int(cds_coordinates[1])
                        print(min, max)
                        trascipt_seq  = trascipt_seq  + seq[min-1:max]
                    print(gene_id, trascipt_seq)
                    print(self.get_codon(trascipt_seq, gene_id))




        '''
        codon_list = []
        mutation_codon_data = self.getmutation_table()

        for i in range(int(len(seq) / 3)):
            codon = []
            start = 3 * i
            end = 3 * i + 3
            codon.append(gene_id)
            triplet = seq[start:end]
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
            codon.append((cds_start -1) + start +1)
            codon.append((cds_start - 1) + end)
            codon.append(str((cds_start -1) + 3 * i + 1) + ", " + str((cds_start -1) + 3 * i + 2) + ", " + str((cds_start -1) +3 * i + 3))
            codon.append(N)
            codon.append(S)
            codon_list.append(codon)'''
        #return codon_list

    def read_gff_file(self, gff_file):
        '''read gtf file'''
        chr_dict = {}

        with open(gff_file) as fp:
            for line in fp:
                if not line.startswith("#"):
                    line = line.strip()
                    rec = line.split("\t")
                    chr = rec[0]

                    if(rec[2] == "gene"):
                        print(rec[8])
                        id = (rec[8]).split(";")[0]
                        gene_name = id.split("=")[1]
                        print(gene_name)
                        gene_dict = {}
                        cds_flag = 0
                        print(chr)

                    if(rec[2] == "CDS"):
                        cds_start = rec[3]
                        cds_end = rec[4]
                        cds_flag = 1
                        if(gene_name in gene_dict.keys()):
                            gene_dict[gene_name].append([cds_start, cds_end])
                        else:
                            gene_dict[gene_name] = [[cds_start, cds_end]]

                    if(cds_flag):
                        if(chr in chr_dict.keys()):
                            if(gene_dict not in chr_dict[chr]):
                                chr_dict[chr].append(gene_dict)
                        else:
                            chr_dict[chr] = [gene_dict]

        return chr_dict

    def filter_ann(self, ann_field):
        if '>' in ann_field:
            return True
        else:
            return False

    def read_vcf(self, vcf_file):
        '''parse vcf file'''

        varlist = []
        with open(vcf_file) as fp:
            line = fp.readline()
            while line:
                if not line.startswith("#"):
                    var = []
                    coverage = {"A": 0, "C": 0, "G": 0, "T": 0}
                    line = line.strip()
                    rec = line.split("\t")
                    var.append(rec[0])
                    var.append(rec[3])
                    var.append(rec[4])
                    var.append(rec[1])
                    annotation = rec[7]
                    var_field = filter(self.filter_ann, annotation.split("|"))
                    seq = self.read_refseq("sample.fa")
                    pos_in_codon = 0
                    for variation in var_field:
                        pos_in_codon = variation[2:(variation.find('>') - 1)]
                    var.append(pos_in_codon)  # get from snpeff results
                    mod = int(pos_in_codon) % 3
                    quotient = int(pos_in_codon) // 3
                    if (mod == 0):
                        quotient = quotient - 1
                        codon_start = quotient * 3 + 1
                    else:
                        codon_start = quotient * 3 + 1
                    codon = seq[codon_start - 1:codon_start + 2]
                    if (mod == 0):
                        mod = 3
                    var.append(mod)
                    var.append(codon_start)
                    var.append(codon)  # get from snpeff results

                    var_class = annotation.split("|")[1]

                    var.append(var_class)  # get from snpeff results
                    format = (rec[8]).split(":")
                    DP_index = format.index("DP")
                    AD_index = format.index("AD")
                    sample = (rec[9]).split(":")

                    coverage[rec[3]] = (sample[AD_index]).split(",")[0]
                    coverage[rec[4]] = (sample[AD_index]).split(",")[1]  # for single allele
                    var.append(coverage)
                    varlist.append(var)
                line = fp.readline()
        return varlist

    def read_refseq(self, fasta_file):
        '''read fasta file'''

        seq = ''
        with open(fasta_file) as fp:
            line = fp.readline()
            while line:
                if not line.startswith(">"):
                    line = line.strip()
                    seq = seq + line
                line = fp.readline()
        return (seq)

if __name__ == "__main__":
    pu = parse_dataUtils()
    dir = "/Users/manishkumar/Desktop/apps/SNPGenie/Analysis_with_ADinInfo/code/data/two_gene_one_transcript_multiple_cds_one_positive_one_negative"
    sequence =  pu.read_refseq(os.path.join(dir, "sample.fa"))
    vcf_data = pu.read_vcf(os.path.join(dir, "sample.ann.vcf"))
    gff_data = pu.read_gff_file(os.path.join(dir, "sample_file.gff"))

    pu.get_triplets(sequence, gff_data)

    #print(gff_data)

    #print(sequence)

    
