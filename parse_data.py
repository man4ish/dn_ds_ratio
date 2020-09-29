
import json
import csv
import os

class parse_dataUtils:

    def getmutation_table(self):
        '''reading mutation codon'''

        with open('json_data/mutation_codon.json') as mcf:
            mutation_codon_data = json.load(mcf)
        return mutation_codon_data

    def get_codon(self, chr, seq, gene_id, diffmap):
        codon_list = []
        mutation_codon_data = self.getmutation_table()

        for i in range(int(len(seq) / 3)):
            codon_num = i+1
            codon = []
            start = 3 * i
            end = 3 * i + 3
            codon.append(chr)
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
            cds_start = start +1
            cds_end =  end
            codon.append(diffmap[start])
            codon.append(diffmap[start + 2])
            nt_pos1 = diffmap[start]
            nt_pos2 = diffmap[start+1]
            nt_pos3 = diffmap[start + 2]
            codon.append(str(nt_pos1) + ", " + str(nt_pos2) + ", " + str(nt_pos3))
            codon.append(codon_num)
            codon.append(N)
            codon.append(S)
            codon_list.append(codon)

        return codon_list

    def get_triplets(self, seq, chr_dict):
        '''generate triplets from ref seq'''
        codon_list = []
        for chr in chr_dict.keys():
            posmap = {}

            gdiff = 0
            cds_pos = 1
            for genes in chr_dict[chr]:
                trascipt_seq = ''
                diffmap = {}
                cds_global_pos_map = {}
                for gene_id, cds_list in genes.items():

                    c = 1

                    #min = cds_list[0][0]
                    #print("min = " + min)

                    for cds_coordinates in cds_list:
                        if(c > 1):
                            #print(cds_coordinates)
                            #print(c)
                            #print(cds_list[c-1])
                            #print(cds_list[c-2])
                            diff  = int((cds_list[c-1])[0]) - int((cds_list[c-2])[1]) -1

                            #print(diff)
                            gdiff = gdiff + diff
                            #print(gdiff)
                            #print(">>>>>>>")
                        min = int(cds_coordinates[0])
                        max = int(cds_coordinates[1])

                        for i in range(min, max+1):
                            #posmap[i] = i
                            diffmap[cds_pos] = gdiff
                            #
                            cds_pos = cds_pos + 1
                        c = c + 1

                        trascipt_seq  = trascipt_seq  + seq[min-1:max]
                    #print(gene_id, trascipt_seq)
                    #print(posmap)
                    #print(diffmap)

                    fkey  = list(diffmap.keys())[0]

                    count = 0
                    for key in diffmap.keys():
                        #print(str(key-fkey ) + "\t" + str(diffmap[key] + fkey + count))
                        cds_global_pos_map[key-fkey ] = diffmap[key] + fkey + count
                        count = count + 1
                    #print(cds_global_pos_map)
                    codon_list.append(self.get_codon(chr, trascipt_seq, gene_id, cds_global_pos_map))
        return codon_list

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

    def filter_codon_change(self, ann_field):
        if 'p.' in ann_field:
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
                    print(line)
                    var = []
                    coverage = {"A": 0, "C": 0, "G": 0, "T": 0}
                    line = line.strip()
                    rec = line.split("\t")
                    var.append(rec[0])
                    var.append(rec[3])
                    var.append(rec[4])
                    #var.append(rec[1])
                    annotation = rec[7]
                    var_field = filter(self.filter_ann, annotation.split("|"))
                    mut_filed = filter(self.filter_codon_change, annotation.split("|"))


                    #exit(annotation.split("|"))
                    seq = self.read_refseq("sample.fa")
                    pos_in_codon = 0
                    field_num = 0
                    pos_field = 0
                    mutation_field = ""
                    for variation in var_field:
                        pos_in_codon = variation[2:(variation.find('>') - 1)]
                        mutation_field += variation
                        print(variation)
                    var.append(pos_in_codon)  # get from snpeff results
                    print(mutation_field)

                    codon_number = 0
                    for mutation in mut_filed:
                        codon_number = mutation[5:6]
                        #print(mutation)

                    var.append(codon_number)


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
    var_list = pu.read_vcf(os.path.join(dir, "sample.ann.vcf"))

    if os.path.exists("variant_info.tsv"):
        os.remove("variant_info.tsv")
    with open('variant_info.tsv', 'w') as variant_tmp_file:
        var_temp = csv.writer(variant_tmp_file, delimiter='\t')
        var_temp.writerow(["#chr", "ref", "alt", "pos", "codon number", "pos in codon", "codon start", "codon", "mutation type", "coverage"])
        for var_gene_list in var_list:
            var_temp.writerow(var_gene_list)

    gff_data = pu.read_gff_file(os.path.join(dir, "sample_file.gff"))

    codon_list = pu.get_triplets(sequence, gff_data)
    #print(codon_list)

    if os.path.exists("codon_results_temp.tsv"):
        os.remove("codon_results_temp.tsv")
    with open('codon_results_temp.tsv', 'w') as cdr_tmp_file:
        cdr_temp = csv.writer(cdr_tmp_file, delimiter='\t')
        cdr_temp.writerow(["#chr", "gene", "codon", "codon start", "codon end", "codon positions", "codon number", "N", "S"])
        for gene_codon_list in codon_list:
            for codon in gene_codon_list:
                cdr_temp.writerow(codon)


