import math
SAMPLES = ["JBM10","JBNZ41","JBC07"] #"JBM10","JBC07","JBNZ41","DS","933-7"

rule all:
    input:
        #expand("../results/extracted_genes/{sample}_matches.m8.trim2", sample=SAMPLES)
        expand("../results/extracted_genes/{sample}_merged.csv", sample=SAMPLES)

rule diamond:
    input:
        genes = "/media/tintri/results/extracted_genes/{sample}_genes.fa",
        database = "/media/tintri/code/taxmapper_supplement/databases/uniprot/uniprot_db.dmnd"
    output:
        "../results/extracted_genes/{sample}_matches.m8"
    threads: 16
    run:
        shell("echo '/home/sm/test/software/diamond/diamond blastx -p {threads} -d {input.database} -o {output} -q {input.genes}'")
        shell("/home/sm/test/software/diamond/diamond blastx -p {threads} -d {input.database} -o {output} -q {input.genes}")


rule trim_input_summary:
    input:
        "../results/extracted_genes/{sample}_matches.m8"
    output:
        temp("../results/extracted_genes/{sample}_matches.m8.trim")
    run:
        for entry in {input}:
            with open(entry[0],"r") as infile:
                read_dict = {}
                for line in infile:
                    line_list = line.strip().split("\t")
                    # read_name       ids   identity        aln_length     (e-value)[10]    bit-score[-1]
                    # gene_2|GeneMark.hmm|975_nt|dna.fa_5677	sp|B5XZG7|SYQ_KLEP3	46.7	321	145	7	1	963	257	551	3.7e-78	293.1
                    #try:
                    #if not line.startswith("read_name"):
                    if float(line_list[10]) < (1e-5):
                        if float(line_list[2]) > 30 and math.log10(float(line_list[10])+(1e-100)) < -1: # Identity > 30% and log10(e-value) < -1 # +(1e-100) because log of value 0 not calculable
                            if not line_list[0] in read_dict.keys():
                                read_dict[line_list[0]] = line_list
                            else:
                                #if (float(read_dict[line_list[0]][10])  <= (float(line_list[10]))):
                                if (float(read_dict[line_list[0]][10]) / math.log(int(read_dict[line_list[0]][3]))) <= (float(line_list[10]) / math.log(int(line_list[3]))): # check if e-value of next hit is higher
                                    continue
                                else:
                                    read_dict[line_list[0]] = line_list
                    #except ValueError:
                    #    pass
            with open(entry[0]+".trim","w") as outfile:
                for ekey in read_dict.keys():
                    outfile.write("\t".join(read_dict[ekey])+"\n")

rule merge:
    input:
        "../results/extracted_genes/{sample}_matches.m8.trim"
    output:
        "../results/extracted_genes/{sample}_matches.m8.trim2"
    run:
        # assign ko id
        with open("/media/tintri/code/taxmapper_supplement/databases/uniprot/uniprot2ko.txt","r") as ko_info:
            dict_ko = {}
            for line in ko_info:
                line = line.strip().split("\t")
                dict_ko[line[0]] = line[1]
        # assign gene name
        with open("/media/tintri/code/taxmapper_supplement/databases/uniprot/uniprot2gene.txt","r") as gene_info:
            dict_gene = {}
            for line in gene_info:
                line = line.strip().split("\t")
                dict_gene[line[0]] = line[1]
        # merge files
        for entry in {input}:
            with open(entry[0], "r") as infile:
                with open(entry[0] + "2", "w") as outfile:
                    read_dict = {}
                    for line in infile:
                        line_list = line.strip().split("\t")
                        ref = line_list[1].split("|")[1]
                        try:
                            new_ref = "|".join((line_list[1],dict_ko[ref],dict_gene[ref]))
                            entry = "\t".join((line_list[0],new_ref,"\t".join(line_list[2:])))
                            outfile.write(entry+"\n")
                        except KeyError:
                            try:
                                new_ref = "|".join((line_list[1], dict_ko[ref]))
                                entry = "\t".join((line_list[0], new_ref, "\t".join(line_list[2:])))
                                outfile.write(entry + "\n")
                            except KeyError:
                                try:
                                    new_ref = "|".join((line_list[1],"", dict_gene[ref]))
                                    entry = "\t".join((line_list[0], new_ref, "\t".join(line_list[2:])))
                                    outfile.write(entry + "\n")
                                except KeyError:
                                    print("KeyError on: "+str(ref))

rule merge_uniprot_KEGG:
    input:
        uniprot = "../results/extracted_genes/{sample}_matches.m8.trim2",
        kegg = "../results/extracted_genes/{sample}_complete.csv.trim"
    output:
        "../results/extracted_genes/{sample}_merged.csv"
    run:
        for entry in {input.uniprot}:
            matching_file = entry.split("/")[-1].split("_")[0]
            with open(entry, "r") as infile:
                read_dict = {}
                for line in infile:
                    line_list = line.strip().split("\t")
                    ko_id = line_list[1].split("|")[3]
                    e_value_uni = line_list[-2]
                    if len(ko_id) > 2:
                        read_dict[line_list[0]] = (ko_id,e_value_uni)

            with open("../results/extracted_genes/" + matching_file + "_complete.csv.trim", "r") as kegg_file:
                with open("../results/extracted_genes/" + matching_file + "_merged.csv", "w") as outfile:
                    i = 0
                    for line in kegg_file:
                        gene_info = line.strip().split("\t")
                        gene = gene_info[0]
                        #try:
                        if gene in read_dict.keys():
                            try:
                                kegg = gene_info[13]
                                if kegg.split(":")[1] != read_dict[gene][0]:
                                    i += 1
                                    if float(gene_info[10]) > float(read_dict[gene][1]):
                                        line = line.strip() + "\tko:" + read_dict[gene][0] + "\n"
                                        print("replaced:" + str(kegg) + "\twith\t" + str(read_dict[gene]) + " " + str(i))
                                    else:
                                        print("unequal:"+str(kegg)+"\t"+str(read_dict[gene])+" "+str(i))
                            except IndexError:
                                line = line.strip()+"\tko:"+read_dict[gene][0]+"\n"
                        outfile.write(line)

