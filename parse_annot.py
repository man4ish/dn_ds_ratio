

#ann_string = "ANN=C|start_lost|HIGH|cds00001|GENE_cds00001|transcript|mRNA00001|protein_coding|1/4|c.2T>C|p.Met1?|2/19|2/19|1/5||WARNING_TRANSCRIPT_INCOMPLETE,C|upstream_gene_variant|MODIFIER|cds00005|GENE_cds00005|transcript|mRNA00002|protein_coding||c.-26T>C|||||26|WARNING_TRANSCRIPT_NO_START_CODON,C|upstream_gene_variant|MODIFIER|mrna00001|null|transcript|mrna00002|protein_coding||c.-25T>C|||||26|WARNING_TRANSCRIPT_NO_START_CODON,C|intragenic_variant|MODIFIER|gene00001|gene00001|gene_variant|gene00001|||n.2T>C||||||,C|non_coding_transcript_variant|MODIFIER|mrna00001|null|transcript|mrna00001|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON;LOF=(cds00001|GENE_cds00001|1|1.00)"
#ann_string = "ANN=T|splice_region_variant&synonymous_variant|LOW|cds00001|GENE_cds00001|transcript|mRNA00001|protein_coding|1/4|c.4C>T|p.Leu2Leu|4/19|4/19|2/5||WARNING_TRANSCRIPT_INCOMPLETE,T|upstream_gene_variant|MODIFIER|cds00005|GENE_cds00005|transcript|mRNA00002|protein_coding||c.-24C>T|||||24|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|mrna00001|null|transcript|mrna00002|protein_coding||c.-23C>T|||||24|WARNING_TRANSCRIPT_NO_START_CODON,T|intragenic_variant|MODIFIER|gene00001|gene00001|gene_variant|gene00001|||n.4C>T||||||,T|non_coding_transcript_variant|MODIFIER|mrna00001|null|transcript|mrna00001|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON"
#ann_string = "ANN=C|splice_acceptor_variant&splice_region_variant&intron_variant|HIGH|cds00001|GENE_cds00001|transcript|mRNA00001|protein_coding|3/3|c.18-1G>C||||||WARNING_TRANSCRIPT_INCOMPLETE,T|splice_acceptor_variant&splice_region_variant&intron_variant|HIGH|cds00001|GENE_cds00001|transcript|mRNA00001|protein_coding|3/3|c.18-1G>T||||||WARNING_TRANSCRIPT_INCOMPLETE,C|upstream_gene_variant|MODIFIER|cds00005|GENE_cds00005|transcript|mRNA00002|protein_coding||c.-3G>C|||||3|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|cds00005|GENE_cds00005|transcript|mRNA00002|protein_coding||c.-3G>T|||||3|WARNING_TRANSCRIPT_NO_START_CODON,C|upstream_gene_variant|MODIFIER|mrna00001|null|transcript|mrna00002|protein_coding||c.-2G>C|||||3|WARNING_TRANSCRIPT_NO_START_CODON,T|upstream_gene_variant|MODIFIER|mrna00001|null|transcript|mrna00002|protein_coding||c.-2G>T|||||3|WARNING_TRANSCRIPT_NO_START_CODON,C|intragenic_variant|MODIFIER|gene00001|gene00001|gene_variant|gene00001|||n.25G>C||||||,T|intragenic_variant|MODIFIER|gene00001|gene00001|gene_variant|gene00001|||n.25G>T||||||,C|non_coding_transcript_variant|MODIFIER|mrna00001|null|transcript|mrna00001|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON,T|non_coding_transcript_variant|MODIFIER|mrna00001|null|transcript|mrna00001|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON;LOF=(cds00001|GENE_cds00001|1|1.00)"
ann_string = "ANN=A|missense_variant&splice_region_variant|MODERATE|cds00005|GENE_cds00005|transcript|mRNA00002|protein_coding|1/4|c.3C>A|p.His1Gln|3/18|3/18|1/5||WARNING_TRANSCRIPT_NO_START_CODON,G|missense_variant&splice_region_variant|MODERATE|cds00005|GENE_cds00005|transcript|mRNA00002|protein_coding|1/4|c.3C>G|p.His1Gln|3/18|3/18|1/5||WARNING_TRANSCRIPT_NO_START_CODON,A|downstream_gene_variant|MODIFIER|cds00001|GENE_cds00001|transcript|mRNA00001|protein_coding||c.*3C>A|||||3|WARNING_TRANSCRIPT_INCOMPLETE,G|downstream_gene_variant|MODIFIER|cds00001|GENE_cds00001|transcript|mRNA00001|protein_coding||c.*3C>G|||||3|WARNING_TRANSCRIPT_INCOMPLETE,A|downstream_gene_variant|MODIFIER|mrna00001|null|transcript|mrna00001|protein_coding||c.*3C>A|||||3|WARNING_TRANSCRIPT_NO_START_CODON,G|downstream_gene_variant|MODIFIER|mrna00001|null|transcript|mrna00001|protein_coding||c.*3C>G|||||3|WARNING_TRANSCRIPT_NO_START_CODON,A|intragenic_variant|MODIFIER|gene00002|gene00002|gene_variant|gene00002|||n.30C>A||||||,G|intragenic_variant|MODIFIER|gene00002|gene00002|gene_variant|gene00002|||n.30C>G||||||,A|non_coding_transcript_variant|MODIFIER|mrna00001|null|transcript|mrna00002|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON,G|non_coding_transcript_variant|MODIFIER|mrna00001|null|transcript|mrna00002|protein_coding||||||||WARNING_TRANSCRIPT_NO_START_CODON"


def split_and_check(annot):
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
    		i = i +1
    
    if (i > 0):
       return i
    else:
       return None
   

def parse_annotation(info_string):
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
            allele = eff[0].replace("ANN=","")
            annot = eff[1]

            if (split_and_check(annot) is None):
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


ann = parse_annotation(ann_string)
print (ann)


#print ("\n\n\nPrinting split stuff\n\n\n")

#x = ann_string.split (",")
#for j in x:
#	print (j)






