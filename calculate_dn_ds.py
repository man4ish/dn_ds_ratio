#!/usr/bin/python
import json
import csv
import os

class calculate_dNdS_Utils:
    def get_codon(self):
        ''' readinfg codon json'''

        with open('json_data/codon.json') as cf:
            codon_data = json.load(cf)
        return codon_data

    def get_codon(self):
        ''' readinfg codon json'''

        with open('json_data/codon.json') as cf:
            codon_data = json.load(cf)
        return codon_data

    def getmutation_table(self):
        '''reading mutation codon'''

        with open('json_data/mutation_codon.json') as mcf:
            mutation_codon_data = json.load(mcf)
        return mutation_codon_data

    def get_all_possible_path(self):
        '''reading all possible path'''

        with open('json_data/all_possible_path.json') as appf:
            all_possible_path = json.load(appf)
        return all_possible_path

    def calculate_dn_ds_ratio(self, codon_list):
        '''
        :param codon_list: list of codons with Ndiffs, Sdiffs, Nsites and Ssites values
        :return: ratio of dn and ds
        '''

        total_Ndiffs = 0
        total_Sdiffs = 0
        total_Nsites = 0
        total_Ssites = 0

        for cdn in codon_list:
            total_Nsites = total_Nsites + cdn[0]
            total_Ssites = total_Ssites + cdn[1]
            total_Ndiffs = total_Ndiffs + cdn[2]
            total_Sdiffs = total_Sdiffs + cdn[3]

        pn = total_Ndiffs / total_Nsites
        ps = total_Sdiffs / total_Ssites

        print("Ndiiffs = " + str(total_Ndiffs))
        print("Sdiiffs = " + str(total_Sdiffs))
        print("Nsites = " + str(total_Nsites))
        print("Ssites = " + str(total_Ssites))
        print("pn = " + str(pn))
        print("ps = " + str(ps))

        dn_ds_ratio = pn / ps

        return dn_ds_ratio

    def get_all_possible_combination(self, n):
        ''' get all possible combination formula = nc2'''
        return n * (n - 1) / 2

    def get_Ndiffs(self, coverage_dict):
        '''
        :param coverage_dict: dictionary of coverage
        :return: Ndiffs or Sdiffs
        '''
        num_A = int(coverage_dict['A'])
        num_C = int(coverage_dict['C'])
        num_G = int(coverage_dict['G'])
        num_T = int(coverage_dict['T'])
        DP = num_A + num_C + num_G + num_T

        num_AC = num_A * num_C
        num_AG = num_A * num_G
        num_AT = num_A * num_T
        num_CG = num_C * num_G
        num_CT = num_C * num_T
        num_GT = num_G * num_T

        possible_combination = self.get_all_possible_combination(DP)
        # coverage_list = [num_AC, num_AG, num_AT, num_CG, num_CT, num_GT]
        total_coverage_product = num_AC + num_AG + num_AT + num_CG + num_CT + num_GT
        Ndiffs = total_coverage_product / possible_combination

        return Ndiffs

    def get_allele_freq(self, pos, alleles):
        ''' get allele frquency for given pos'''
        for allele in alleles:
            if pos in allele.keys():
                # print(allele)
                alt_frq = allele[pos][1]['alt']['freq']
                ref_frq = allele[pos][0]['ref']['freq']
                allele_prop = int(alt_frq) / (int(ref_frq) + int(alt_frq))
                ref_prop = int(ref_frq) / (int(ref_frq) + int(alt_frq))
                # print(ref_prop)
                allele_dict = {"ref_freq": ref_prop, "alt_freq": allele_prop}
                return allele_dict
            else:
                allele_dict = {"ref_freq": 1, "alt_freq": 0}

        return allele_dict

    def possible_codon(self, key, alleles):
        '''generate a list of possible codon with given position and allele list'''

        codon = {}

        ref = key.split("-")[2]
        for i in range(3):
            codon[i + 1] = [ref[i]]

        for allele in alleles:
            for pos in allele:
                codon[pos].append(allele[pos][1]['alt']['allele'])

        possible_cdn_list = []
        possible_cdn_freq_list = []

        nt_frq1 = self.get_allele_freq(1, alleles)
        nt_frq2 = self.get_allele_freq(2, alleles)
        nt_frq3 = self.get_allele_freq(3, alleles)

        nt_pos1 = False
        nt_pos2 = False
        nt_pos3 = False

        for i in range(len(codon[1])):
            for j in range(len(codon[2])):
                for k in range(len(codon[3])):
                    pos1 = 1
                    pos2 = 2
                    pos3 = 3
                    nt1 = codon[1][i]
                    nt2 = codon[2][j]
                    nt3 = codon[3][k]

                    for allele in alleles:
                        if (pos1 in allele.keys()):
                            ref1 = allele[pos1][0]['ref']['allele']

                            if (nt1 == ref1):
                                nt_prop1 = nt_frq1["ref_freq"]
                            else:
                                nt_prop1 = nt_frq1["alt_freq"]
                            nt_pos1 = True

                        if (pos2 in allele.keys()):
                            ref2 = allele[pos2][0]['ref']['allele']

                            if (nt2 == ref2):
                                # print(codon[2][j] + "\t" + ref2)
                                # print(nt_frq2)
                                nt_prop2 = nt_frq2["ref_freq"]
                                # print(nt_prop2)
                            else:
                                nt_prop2 = nt_frq2["alt_freq"]
                                # print(nt_prop2)
                            nt_pos2 = True

                        if (pos3 in allele.keys()):
                            ref3 = allele[pos3][0]['ref']['allele']

                            if (nt3 == ref3):
                                nt_prop3 = nt_frq3["ref_freq"]
                            else:
                                nt_prop3 = nt_frq3["alt_freq"]
                            nt_pos3 = True

                        if (not nt_pos1):
                            nt_prop1 = 1
                        if (not nt_pos2):
                            nt_prop2 = 1
                        if (not nt_pos3):
                            nt_prop3 = 1

                    possible_cdn_freq_list.append([nt_prop1, nt_prop2, nt_prop3])
                    possible_cdn_list.append(nt1 + nt2 + nt3)

        possible_codon_dict = {"codon": possible_cdn_list, "allele_prop": possible_cdn_freq_list}

        return possible_codon_dict

    def get_Sites_ref(self, triplet):
        '''Calculate Nsites_ref and Ssites_ref from mutation table'''

        N = 0
        S = 0
        mutation_codon_data = self.getmutation_table()
        N1 = mutation_codon_data[triplet + "_N_1"]
        N2 = mutation_codon_data[triplet + "_N_2"]
        N3 = mutation_codon_data[triplet + "_N_3"]
        S1 = mutation_codon_data[triplet + "_S_1"]
        S2 = mutation_codon_data[triplet + "_S_2"]
        S3 = mutation_codon_data[triplet + "_S_3"]
        N = N + (N1 + N2 + N3)
        S = S + (S1 + S2 + S3)

        return {"Nsites_ref": N, "Ssites_ref": S}

    def get_Sites(self, all_codon):
        '''
        :param all_codon:
        :return: Nsites and Ssites
        '''

        c = 0
        Nsites = 0
        Ssites = 0
        for triplet in all_codon["codon"]:
            print(triplet)
            nt_freqs = all_codon["allele_prop"][c]
            freq_product = nt_freqs[0] * nt_freqs[1] * nt_freqs[2]
            # print(freq_product)
            c = c + 1
            Sites = self.get_Sites_ref(triplet)
            Nsites_ref = Sites['Nsites_ref']
            Ssites_ref = Sites['Ssites_ref']
            Nsites = Nsites + freq_product * Nsites_ref
            Ssites = Ssites + freq_product * Ssites_ref

        #print(str(Nsites) + "\t" + str(Ssites))
        return {"Nsites": Nsites, "Ssites": Ssites}

    def get_all_possible_codon(self, codon_list):
        '''get all possible codon with given allele'''

        cdn_dict = {}
        for cdn_lst in codon_list:
            key = str(cdn_lst[0]) + "-" + str(cdn_lst[6]) + "-" + cdn_lst[7]
            if not key in cdn_dict:

                cdn_dict[key] = []
                pos_dict = {}
                pos_in_codon = cdn_lst[5]
                ref_base = cdn_lst[1]
                freq = cdn_lst[9]
                alt_base = cdn_lst[2]
                allele_dict = [{'ref': {'allele': ref_base, 'freq': freq[ref_base]}},
                               {'alt': {'allele': alt_base, 'freq': freq[alt_base]}}]

                pos_dict[pos_in_codon] = allele_dict
                cdn_dict[key].append(pos_dict)

            else:
                pos_dict = {}
                pos_in_codon = cdn_lst[5]
                ref_base = cdn_lst[1]
                freq = cdn_lst[9]
                alt_base = cdn_lst[2]
                allele_dict = [{'ref': {'allele': ref_base, 'freq': freq[ref_base]}},
                               {'alt': {'allele': alt_base, 'freq': freq[alt_base]}}]
                pos_dict[pos_in_codon] = allele_dict

                cdn_dict[key].append(pos_dict)

        for key in cdn_dict:
            all_codon = self.possible_codon(key, cdn_dict[key])
            print("****")
            self.get_Sites(all_codon)
            print(all_codon)
            print("****")

        return cdn_dict

    def get_triplets(self, seq, gene_id):
        '''generate triplets from ref seq'''

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
            codon.append(start + 1)
            codon.append(end)
            codon.append(str(3 * i + 1) + ", " + str(3 * i + 2) + ", " + str(3 * i + 3))
            codon.append(N)
            codon.append(S)
            codon_list.append(codon)

        return codon_list

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

    def gen_codonlist(self, seq, start, stop, gene_id):
        '''generate codon list for gene'''

        subseq = seq[start - 1:stop]
        return self.get_triplets(subseq, gene_id)

    def read_gff_file(self, gff_file):
        '''read gtf file'''

        with open(gff_file) as fp:
            line = fp.readline()
            while line:
                if not line.startswith("#"):
                    line = line.strip()
                    rec = line.split("\t")
                    cds_start = rec[3]
                    cds_end = rec[4]
                    id = rec[8]
                gene_id = id.split("\"")[1]
                line = fp.readline()

    def filter_ann(self, ann_field):
        if '>' in ann_field:
            return True
        else:
            return False

    def read_vcf(self, vcf_file):
        '''parse vcf file'''

        varlist = []
        # ['Chr01', '8', '.', 'C', 'G', '65.28', 'FS_filter;SOR_filter', '', 'GT:DP:AD', '0/1:300:100,200']
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
                    print(annotation.split("|"))
                    var_field = filter(self.filter_ann, annotation.split("|"))
                    seq = self.read_refseq("sample.fa")
                    pos_in_codon = 0
                    for variation in var_field:
                        pos_in_codon = variation[2:(variation.find('>') - 1)]
                        print(pos_in_codon)
                    var.append(pos_in_codon)  # get from snpeff results
                    mod = int(pos_in_codon) % 3
                    quotient = int(pos_in_codon) // 3
                    if (mod == 0):
                        quotient = quotient - 1
                        print("remainder=" + pos_in_codon + "\t" + str(quotient))
                        codon_start = quotient * 3 + 1
                    else:
                        codon_start = quotient * 3 + 1
                    print("codon_start=" + str(codon_start))
                    print(seq)
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

                    print("DP=" + sample[DP_index] + "\tAD=" + sample[AD_index])
                    coverage[rec[3]] = (sample[AD_index]).split(",")[0]
                    coverage[rec[4]] = (sample[AD_index]).split(",")[1]  # for single allele
                    print(coverage)
                    var.append(coverage)
                    varlist.append(var)
                line = fp.readline()
        return varlist

if __name__ == "__main__":
    cu = calculate_dNdS_Utils()
    gff_file = cu.read_gff_file("sample.gtf")
    seq = cu.read_refseq("sample.fa")

    data_dict = {}

    data = cu.read_vcf("snpeff_sample.ann.vcf")
    if os.path.exists("variant_info.tsv"):
        os.remove("variant_info.tsv")
    with open('variant_info.tsv', 'a') as myfile:
        wr = csv.writer(myfile, delimiter='\t')
        for data_list in data:
            key = data_list[0] + '-' + str(data_list[6]) + '-' + data_list[7]
            print(key)
            if (key in data_dict.keys()):
                data_dict[key].append([data_list[8], data_list[9]])
            else:
                data_dict[key] = [[data_list[8], data_list[9]]]
            wr.writerow(data_list)

    all_possible_codon = cu.get_all_possible_codon(data)
    #print(all_possible_codon)

    result_list = cu.gen_codonlist(seq, 1, 18, "Chr01")

    if os.path.exists("codon_results.tsv"):
        os.remove("codon_results.tsv")
    with open('codon_results.tsv', 'a') as cdr_file:
        cdr = csv.writer(cdr_file, delimiter='\t')
        pn_ps_lst = []
        for cd_list in result_list:
            key = cd_list[0] + '-' + str(cd_list[2]) + '-' + cd_list[1]

            if (key in all_possible_codon.keys()):
                all_codon = cu.possible_codon(key, all_possible_codon[key])
                print("****")
                # print(cdn_dict[key])
                Sites = cu.get_Sites(all_codon)
                Nsites = Sites['Nsites']
                Ssites = Sites['Ssites']
                Ndiffs = 0
                Sdiffs = 0
                ###### Calculate Ndiffs and Sdiffs #####
                if (key in data_dict.keys()):
                    variants = data_dict[key]
                    print(variants)

                    for var in variants:

                        mutation_type = var[0]
                        coverage = var[1]
                        if (mutation_type == "synonymous_variant"):
                            Sdiffs = cu.get_Ndiffs(coverage)
                        else:
                            Ndiffs = cu.get_Ndiffs(coverage)

                print("Sdiifs " + str(Sdiffs) + "\t" + str(Ndiffs))

                cd_list.append(Nsites)
                cd_list.append(Ssites)
                cd_list.append(Ndiffs)
                cd_list.append(Sdiffs)
                pn_ps_lst.append([Nsites, Ssites, Ndiffs, Sdiffs])
            else:
                Nsites = cd_list[5]
                Ssites = cd_list[6]
                cd_list.append(Nsites)
                cd_list.append(Ssites)
                cd_list.append(0)
                cd_list.append(0)
                pn_ps_lst.append([Nsites, Ssites, 0, 0])

            cdr.writerow(cd_list)

        dn_ds_ratio = cu.calculate_dn_ds_ratio(pn_ps_lst)
        print("dn/ds ratio = " + str(dn_ds_ratio))  # need to check with spreasheet

