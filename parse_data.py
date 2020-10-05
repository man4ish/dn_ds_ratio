import json
import csv
import os
import re
import json
class parse_dataUtils:

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



        with open("dnds_statistics.tsv", "w") as stat_file:
            for cdn in codon_list:
                total_Nsites = total_Nsites + float(cdn[0])
                total_Ssites = total_Ssites + float(cdn[1])
                total_Ndiffs = total_Ndiffs + float(cdn[2])
                total_Sdiffs = total_Sdiffs + float(cdn[3])

            pn = total_Ndiffs / total_Nsites
            ps = total_Sdiffs / total_Ssites
            #dn_ds_ratio = 0  # hardcoded for testing
            dn_ds_ratio = pn / ps

            stat_file.write("Ndiiffs:\t" + str(total_Ndiffs) + "\n")
            stat_file.write("Sdiiffs:\t" + str(total_Sdiffs)+ "\n")
            stat_file.write("Nsites:\t" + str(total_Nsites)+ "\n")
            stat_file.write("Ssites:\t" + str(total_Ssites)+ "\n")
            stat_file.write("pn:\t" + str(pn)+ "\n")
            stat_file.write("ps:\t" + str(ps)+ "\n")
            stat_file.write("dn/ds:\t" + str(dn_ds_ratio)+ "\n")

        return dn_ds_ratio

    def get_all_possible_combination(self, n):
        ''' get all possible combination formula = nc2'''
        return n * (n - 1) / 2

    def get_Ndiffs(self, coverage):
        '''
        :param coverage_dict: dictionary of coverage
        :return: Ndiffs or Sdiffs
        '''
        coverage = coverage.replace("'", '"')
        coverage_dict = json.loads(coverage)


        print(coverage)
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

        print(str(num_AC) + "\t"+  str(num_AG) + "\t" + str(num_AT) +"\t" + str(num_CG) + "\t"+  str(num_CT) +"\t"+ str(num_GT))
        possible_combination = self.get_all_possible_combination(DP)

        total_coverage_product = num_AC + num_AG + num_AT + num_CG + num_CT + num_GT

        print("****" + str(total_coverage_product))
        Ndiffs = total_coverage_product / possible_combination
        print("(((((" + str(Ndiffs))

        return Ndiffs

    def get_Sites_ref(self, triplet):
        '''Calculate Nsites_ref and Ssites_ref from mutation table'''

        N = 0
        S = 0
        mutation_codon_data = self.getmutation_table()

        N1 = 0
        N2 = 0
        N3 = 0
        S1 = 0
        S2 = 0
        S3 = 0

        if(triplet + "_N_1" in mutation_codon_data.keys()):
            N1 = mutation_codon_data[triplet + "_N_1"]
        if (triplet + "_N_2" in mutation_codon_data.keys()):
            N2 = mutation_codon_data[triplet + "_N_2"]
        if (triplet + "_N_3" in mutation_codon_data.keys()):
            N3 = mutation_codon_data[triplet + "_N_3"]
        if (triplet + "_S_1" in mutation_codon_data.keys()):
            S1 = mutation_codon_data[triplet + "_S_1"]
        if (triplet + "_S_2" in mutation_codon_data.keys()):
            S2 = mutation_codon_data[triplet + "_S_2"]
        if (triplet + "_S_3" in mutation_codon_data.keys()):
            S3 = mutation_codon_data[triplet + "_S_3"]

        N = N + (N1 + N2 + N3)
        S = S + (S1 + S2 + S3)

        return {"Nsites_ref": N, "Ssites_ref": S}

    def get_Sites(self, all_codon):
        '''
        :param all_codon:
        :return: Nsites and Ssites
        '''

        if(len(all_codon) == 0):
            return {"Nsites": 0, "Ssites": 0}

        c = 0
        Nsites = 0
        Ssites = 0
        for triplet in all_codon["codon"]:
            print(triplet)
            nt_freqs = all_codon["allele_prop"][c]

            freq_product = nt_freqs[0] * nt_freqs[1] * nt_freqs[2]
            c = c + 1
            Sites = self.get_Sites_ref(triplet)

            Nsites_ref = Sites['Nsites_ref']
            Ssites_ref = Sites['Ssites_ref']
            Nsites = Nsites + freq_product * Nsites_ref
            Ssites = Ssites + freq_product * Ssites_ref
            print(Nsites)
            print(Ssites)

        return {"Nsites": Nsites, "Ssites": Ssites}

    def get_allele_freq(self, pos, alleles):
        ''' get allele frquency for given pos'''
        for allele in alleles:
            if pos in allele.keys():
                alt_frq = allele[pos][1]['alt']['freq']
                ref_frq = allele[pos][0]['ref']['freq']
                allele_prop = int(alt_frq) / (int(ref_frq) + int(alt_frq))
                ref_prop = int(ref_frq) / (int(ref_frq) + int(alt_frq))
                allele_dict = {"ref_freq": ref_prop, "alt_freq": allele_prop}
                return allele_dict
            else:
                allele_dict = {"ref_freq": 1, "alt_freq": 0}

        return allele_dict

    def possible_codon(self, key, alleles):
        '''generate a list of possible codon with given position and allele list'''

        codon = {}

        possible_cdn_list = []
        ref = key.split("-")[2]
        print(key)

        possible_cdn_freq_list = []
        possible_cdn_list = []

        if not ref:
            possible_cdn_freq_list.append([0, 0, 0])
            possible_cdn_list.append('---')
        else:
            print("&&&" + ref + "&&&&")
            for i in range(3):
                codon[i + 1] = [ref[i]]

            print(codon)
            for allele in alleles:
                for pos in allele:

                    pos = str(pos)
                    # exit(allele[pos][1]['alt']['allele'])
                    print(allele[pos][1]['alt']['allele'])
                    if (pos in codon.keys()):
                        codon[pos].append(allele[pos][1]['alt']['allele'])

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
                                    nt_prop2 = nt_frq2["ref_freq"]
                                else:
                                    nt_prop2 = nt_frq2["alt_freq"]
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
                freq = freq.replace("'",'"')
                print(freq)
                freq_dict = json.loads(freq)

                allele_dict = [{'ref': {'allele': ref_base, 'freq': freq_dict[ref_base]}},
                               {'alt': {'allele': alt_base, 'freq': freq_dict[alt_base]}}]

                pos_dict[pos_in_codon] = allele_dict
                cdn_dict[key].append(pos_dict)

            else:
                pos_dict = {}
                pos_in_codon = cdn_lst[5]
                ref_base = cdn_lst[1]
                freq = cdn_lst[9]
                freq = freq.replace("'", '"')
                print(freq)
                freq_dict = json.loads(freq)
                alt_base = cdn_lst[2]
                allele_dict = [{'ref': {'allele': ref_base, 'freq': freq_dict[ref_base]}},
                               {'alt': {'allele': alt_base, 'freq': freq_dict[alt_base]}}]
                pos_dict[pos_in_codon] = allele_dict

                cdn_dict[key].append(pos_dict)

        for key in cdn_dict:
            all_codon = self.possible_codon(key, cdn_dict[key])
            self.get_Sites(all_codon)

        return cdn_dict

    def getmutation_table(self):
        '''reading mutation codon'''

        with open('json_data/mutation_codon.json') as mcf:
            mutation_codon_data = json.load(mcf)
        return mutation_codon_data

    def reverse_complement(self, cds_seq):
        '''

        :param cds_seq:
        :return:
        '''
        'reverse'
        print(cds_seq)
        mapping = cds_seq.maketrans('ATCG',"TAGC")
        return cds_seq.translate(mapping)[::-1]


    def split_and_check(self, annot):
        ANNOT = {
            "synonymous_variant": 1,
            "missense_variant": 1,
            "stop_gained": 1,
            "start_lost": 1,
            "start_gained": 1,
            "stop_lost": 1
        }

        all_annot = annot.split("&")

        i = 0

        for annot1 in all_annot:
            if annot1 not in ANNOT:
                continue
            else:
                i = i + 1

        if (i > 0):
            return i
        else:
            return None

    def parse_annotation(self, info_string):
        """
        parses annotation string into a structure
        [gene_id, transcript_id,
        :param ann_string:
        :return:
        """
        #
        # Fields are delimited by ;
        # annotation field starts with ANN=
        # Each allele-effect in annotation field is separated by ,
        annotation_info = list()

        info = info_string.split(";")

        ann_string = None
        for j in info:
            if j.startswith("ANN="):
                ann_string = j

        if ann_string is None:
            return None
        # TODO: Add more variant effects that
        # TODO: may not be affecting the protein coding region

        allele_effect_list = ann_string.split(",")
        for a in allele_effect_list:
            eff = a.split("|")
            allele = eff[0].replace("ANN=", "")
            annot = eff[1]

            if (self.split_and_check(annot) is None):
                continue

            gene_id = eff[3]
            transcript_id = eff[6]
            base = eff[9]
            prot = eff[10]
            annotation_info.append([allele, annot, gene_id, transcript_id, base, prot])

        if annotation_info:
            return annotation_info
        else:
            return None

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
            #TODO: check if key exist otherwise make value = 0


            N1 = mutation_codon_data[triplet + "_N_1"]
            N2 = mutation_codon_data[triplet + "_N_2"]
            N3 = mutation_codon_data[triplet + "_N_3"]
            S1 = mutation_codon_data[triplet + "_S_1"]
            S2 = mutation_codon_data[triplet + "_S_2"]
            S3 = mutation_codon_data[triplet + "_S_3"]
            N = N + (N1 + N2 + N3)
            S = S + (S1 + S2 + S3)

            codon.append(seq[start:end])
            #cds_start = start +1
            #cds_end =  end
            codon.append(diffmap[start])
            codon.append(diffmap[start + 2])

            nt_pos1 = diffmap[start]
            nt_pos2 = diffmap[start+1]
            nt_pos3 = diffmap[start + 2]

            codon.append(str(nt_pos1) + "," + str(nt_pos2) + "," + str(nt_pos3))
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
                for gene_id, cds_dict in genes.items():


                    c = 1

                    #min = cds_list[0][0]
                    #print("min = " + min)
                    cds_list = cds_dict['cds_coordinates']
                    orientation = cds_dict['orient']
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
                    if(orientation == '-'):
                        trascipt_seq = self.reverse_complement(trascipt_seq)
                    #TODO: if orientation is '-' , reverse transcribe the sequence(trascipt_seq)

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
                    print(line)
                    rec = line.split("\t")
                    chr = rec[0]

                    if(rec[2] == "gene"):
                        print(rec[8])
                        id = (rec[8]).split(";")[0]
                        gene_name = id.split("=")[1]
                        orientation = rec[6]
                        print(orientation)
                        print(gene_name)
                        gene_dict = {}
                        cds_flag = 0
                        print(chr)

                    if(rec[2] == "CDS"):
                        cds_start = rec[3]
                        cds_end = rec[4]
                        cds_flag = 1
                        if(gene_name in gene_dict.keys()):
                            #gene_dict[gene_name] = {'orient': orientation, 'cds_coordinates': [[cds_start, cds_end]]}
                            gene_dict[gene_name]['cds_coordinates'].append([cds_start, cds_end])
                        else:
                            gene_dict[gene_name] = {'orient':orientation, 'cds_coordinates': [[cds_start, cds_end]]}
                            #gene_dict[gene_name] = [[cds_start, cds_end]]

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

    def read_vcf(self, vcf_file, seq):
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
                    pos = int(rec[1])
                    var.append(rec[1])

                    annotation = rec[7]
                    print("*******")
                    annot_list = self.parse_annotation(annotation)
                    print(annot_list)
                    cds_list = []
                    cds_pos_list = []
                    cds_num_list = []
                    mutation_type_list = []
                    #allle_dict ={}
                    codon_start_pos_list = []
                    postion_in_codon_list = []
                    if(annot_list):
                        for annot in annot_list:

                            m = re.match( r"\w.(\d+)\w>\w", annot[4])

                            cds_pos = m.groups()[0]
                            pos_in_codon = int(cds_pos) % 3
                            if(pos_in_codon == 0):
                                pos_in_codon = 3
                            codon_start_pos = pos - pos_in_codon + 1
                            codon_start_pos_list.append(str(codon_start_pos))
                            postion_in_codon_list.append(str(pos_in_codon))

                            print(pos_in_codon)
                            cds_pos_list.append(cds_pos)
                            codon_num = (annot[5])[5:6]
                            cds_num_list.append(codon_num)
                            cds_list.append(annot[2])
                            mutation_type_list.append(annot[1])
                        var.append(",".join(cds_num_list))
                        var.append(','.join(postion_in_codon_list))
                        var.append(','.join(codon_start_pos_list))
                    else:
                        var.append('-')
                        var.append('-')
                        var.append('-')
                        mutation_type_list.append('-')





                    print("*******")
                    var_field = filter(self.filter_ann, annotation.split("|"))
                    mut_filed = filter(self.filter_codon_change, annotation.split("|"))


                    #exit(annotation.split("|"))
                    #seq = self.read_refseq("sample.fa")
                    pos_in_codon = 0
                    field_num = 0
                    pos_field = 0
                    mutation_field = ""
                    for variation in var_field:
                        pos_in_codon = variation[2:(variation.find('>') - 1)]
                        mutation_field += variation
                        print(variation)

                    #var.append(pos_in_codon)  # get from snpeff results

                    print(mutation_field)

                    '''
                    codon_number = 0
                    for mutation in mut_filed:
                        codon_number = mutation[5:6]
                        #print(mutation)

                    var.append(codon_number)
                    '''

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
                    #var.append(mod)
                    #var.append(','.join(codon_start_pos_list))
                    var.append(codon)  # get from snpeff results

                    var_class = annotation.split("|")[1]
                    var.append(','.join(mutation_type_list))
                    #var.append(var_class)  # get from snpeff results
                    format = (rec[8]).split(":")
                    DP_index = format.index("DP")
                    AD_index = format.index("AD")
                    sample = (rec[9]).split(":")

                    coverage[str(rec[3])] = int((sample[AD_index]).split(",")[0])
                    coverage[str(rec[4])] = int((sample[AD_index]).split(",")[1])  # for single allele
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



    def merge_files(self, codon_file, variant_file):
        pos_dict = {}
        merged_list = []
        with open(codon_file) as cdfile:
            for line in cdfile:
                if(not line.startswith("#")):
                    line = line.rstrip()
                    cd_rec = line.split()
                    cds_positions = cd_rec[5]
                    cds_pos_llist = cds_positions.split(",")
                    for pos in cds_pos_llist:
                        pos_dict[pos] = line
        print(pos_dict)

        with open("corrected_variant_info.tsv", 'w') as mfile:
            with open(variant_file) as varfile:
                for line in varfile:
                    if (not line.startswith("#")):
                        line = line.rstrip()
                        var_rec = line.split("\t")
                        var_pos = var_rec[3]
                        if(var_pos in pos_dict.keys()):
                            codon = (pos_dict[var_pos]).split("\t")[2]
                            merged_line = line + "\t" + pos_dict[var_pos]
                            var_rec[7] = codon
                            mfile.write('\t'.join(var_rec) + "\n")
                            #merged_list.append(merged_line.split("\t"))
                            #print(merged_line)
                        merged_list.append(var_rec)
                    #print(var_pos)
        return  merged_list

    def generate_statistics(self, varinat_info_file, codon_info_file, all_possible_codon):
        '''
        :param varinat_info_file:
        :param codon_info_file:
        :param all_possible_codon:
        :return:
        '''

        '''read data from variant_info.tsv file'''
        data_dict = {}
        with open(varinat_info_file) as variantfle:
            for var_line in variantfle:
                var_line = var_line.rstrip()
                data_list = var_line.split("\t")
                #for data_list in data:
                key = data_list[0] + '-' + str(data_list[6]) + '-' + data_list[7]
                if (key in data_dict.keys()):
                    data_dict[key].append([data_list[8], data_list[9]])
                else:
                    data_dict[key] = [[data_list[8], data_list[9]]]

        '''read codon table and generate statistics'''
        with open(codon_info_file, 'r') as cdr_file:
            pn_ps_lst = []


            for line in cdr_file:
                if (line.startswith("#")):
                    continue
                line = line.rstrip()
                cd_list = line.split("\t")
                key = cd_list[0] + '-' + cd_list[3] + '-' + cd_list[2]
                print(cd_list)
                if (key in all_possible_codon.keys()):
                    print("****" + key + "****")
                    print(all_possible_codon[key])

                    all_codon = self.possible_codon(key, all_possible_codon[key])
                    #print(all_codon)
                    Sites = self.get_Sites(all_codon)
                    Nsites = Sites['Nsites']
                    Ssites = Sites['Ssites']

                    print(Nsites)
                    print(Ssites)
                    Ndiffs = 0
                    Sdiffs = 0
                    ###### Calculate Ndiffs and Sdiffs #####
                    if (key in data_dict.keys()):
                        variants = data_dict[key]
                        for var in variants:

                            mutation_type = var[0]
                            coverage = var[1]
                            if (mutation_type.find('synonymous_variant') != -1):
                                print("coverage =" + coverage)
                                Sdiffs = self.get_Ndiffs(coverage)
                            else:
                                Ndiffs = self.get_Ndiffs(coverage)

                    pn_ps_lst.append([Nsites, Ssites, Ndiffs, Sdiffs])
                else:
                    Nsites = cd_list[7]
                    Ssites = cd_list[8]
                    pn_ps_lst.append([Nsites, Ssites, 0, 0])

            #exit(pn_ps_lst)
            dn_ds_ratio = self.calculate_dn_ds_ratio(pn_ps_lst)
            print("dn/ds ratio = " + str(dn_ds_ratio))  # need to check with spreasheet

if __name__ == "__main__":
    pu = parse_dataUtils()
    #rev_comp = pu.reverse_complement("ATCGTGTT")


    dir = "/Users/manishkumar/Desktop/apps/SNPGenie/Analysis_with_ADinInfo/code/data/new_dataset"
    sequence =  pu.read_refseq(os.path.join(dir, "sample.fa"))
    var_list = pu.read_vcf(os.path.join(dir, "sample.ann.vcf"), sequence)

    print(var_list)

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

    merged_list = pu.merge_files("codon_results_temp.tsv", "variant_info.tsv")
    all_possible_codon = pu.get_all_possible_codon(merged_list)  # generating all possible codon
    pu.generate_statistics("corrected_variant_info.tsv", "codon_results_temp.tsv", all_possible_codon)
    #exit(all_possible_codon)

    #print(merged_list)

    #possible_codons = pu.get_all_possible_codon(merged_list)



