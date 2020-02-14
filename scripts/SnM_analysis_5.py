#from ../metagenome_filter_module13a.py import SAMPLES3
#configfile: "config.yaml"
#bam=expand("sorted_reads/{sample}.bam", sample=config["samples"])
import math
import re
from collections import defaultdict
#subworkflow gene_annotation:
#    snakefile: "Snakefile_pipe_with_metagenome_klinik_server.py"
#SAMPLES = ["1006","199hm", "AR4D6",  "Chromulina",  "Epipyxis",  "FU18K-A",  "JBC27",  "JBMS11",  "LO226KS",  "LO234KE",  "WA18K-M","DS","933-7"]  #list_of_spezies_names # "1006","199hm",  "AR4D6",  "Chromulina",  "Epipyxis",  "FU18K-A",  "JBC27",  "JBMS11",  "LO226KS",  "LO234KE",  "WA18K-M"
non_axenic_SAMPLES = ["1006","199hm", "AR4D6",  "Chromulina",  "Epipyxis",  "FU18K-A",  "JBC27",  "JBMS11",  "LO226KS",  "LO234KE",  "WA18K-M"] #list_of_spezies_names contaminated with bacteria
axenic_SAMPLES = ["DS","933-7","p.lacustris"] #,"JBC07", "JBM10", "JBNZ41"

SAMPLES3 = set(axenic_SAMPLES+non_axenic_SAMPLES)

rule all:
    input:
        #expand("../results/analysis/{sample}_merged.csv.trim", sample=config["samples"]),
        #expand("../results/analysis/{sample}_complete.csv.trim", sample=config["samples"]),
        #"../results/analysis/kegg/KO_id_2_hierarchy.txt",
        #"../results/analysis/kegg/KO_id_2_hierarchy.trim.tsv",
        expand("../results/analysis/{sample}_ko_functional_groups.list", sample=SAMPLES3),
        #expand("../results/analysis/{sample}_count_B.tsv", sample=SAMPLES3),
        #"/home/stephan/test/results/analysis//sec_papercompare_meta_3.csv"

rule pre_assign_gene2funtion:
    input: "../results/analysis/tmp.file"
    output:
        expand("../results/new_extracted_genes/{sample}_alignment_gene_analysis_functions.txt", sample=SAMPLES3),
        "../results/analysis/p.lacustris_ko.list"
    script:"scripts/extract_gene_functions_assign2alignment_paf_analysis2.py"

rule extract_KO_id:
    input:
        "../results/new_extracted_genes/{sample}_alignment_gene_analysis_functions.txt"
    output:
        "../results/analysis/{sample}_ko.list"
    run:
        with open(input[0],"r") as infile, open(output[0],"w") as outfile:
            for line in infile:
                ko_id = line.split("\t")[2]
                outfile.write(ko_id+"\t"+ko_id+"\n")

rule get_functional_group:
    input:
        strain_ko_ids = "../results/analysis/{sample}_ko.list",
        KO_assignment = "../results/analysis/kegg/KO_id_2_hierarchy.trim.tsv"
    output:
        "../results/analysis/{sample}_ko_functional_groups.list"
    run:
        ko_metadata = {}
        ko_functional_group_count = defaultdict(int)
        ko_functional_group = {}
        with open(str(input.KO_assignment),"r") as infile:
            for line in infile:
                ko_id = line.split()[0]
                ko_metadata[ko_id] = line.strip()
        with open(str(input.strain_ko_ids), "r") as infile:
            for line in infile:
                ko_id = line.split()[0]
                ko_functional_group_count[ko_id] += 1
                try:
                    ko_functional_group[ko_id] = ko_functional_group_count[ko_id], ko_metadata[ko_id]
                except KeyError:
                    print("KeyError: "+ko_id)
                    # KO_id has no assignment to a hierarchy
        with open(output[0],"w") as outfile:
            for entry in ko_functional_group.keys():
                i, string = ko_functional_group[entry]
                outfile.write("\t".join((str(i),string,"\n")))

rule get_KEGG_module_names:
    input:
        "../raw/modules_kegg.txt"
    output:
        "../results/analysis/kegg/KO_id_2_hierarchy.txt" #temp()
    run:
        # each KEGG hierachry gets an assingment
        with open(input[0],"r") as infile, open(output[0],"w") as outfile:
            for line in infile:
                line = line.strip()
                if line.startswith("B"):
                    if len(line) > 1:
                        name = " ".join(line.split()[1:])
                        level_B = re.sub("</?b>","",name)
                elif line.startswith("C"):
                    level_C = " ".join(line.split()[1:])
                elif line.startswith("D"):
                    level_D = " ".join(line.split()[2:])
                    level_D = re.sub("\[.*?\]", "", level_D)
                elif line.startswith("E"):
                    ko_id = line.split()[1]
                    name = " ".join(line.split()[2:])
                    outfile.write("\t".join((ko_id,level_B,level_C,level_D,name))+"\n")

rule remove_pointless_KEGG_pathways:
    input: "../results/analysis/kegg/KO_id_2_hierarchy.txt"
    output: "../results/analysis/kegg/KO_id_2_hierarchy.trim.tsv"
    run:
        level_B_pointless = ["Organismal Systems","Human Diseases","Drug Development","Gene set"]
        level_C_pointless = ["Bacterial secretion system"]
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            for line in infile:
                data = line.split("\t")
                if data[1] not in level_B_pointless:
                    if data[2] not in level_C_pointless:
                        outfile.write(line)

rule count_hierarchy:
    input: "../results/analysis/{sample}_ko_functional_groups.list"
    output:
        a = "../results/analysis/{sample}_count_B.tsv",
    run:
        B_dict = defaultdict(int)
        C_dict = defaultdict(int)
        D_dict = defaultdict(int)
        with open(input[0], "r") as infile:
            for line in infile:
                data = line.strip().split("\t")
                count = int(data[0])
                B_hierarchy = data[2]
                C_hierarchy = data[3]
                D_hierarchy = data[4]
                B_dict[B_hierarchy] += count
                C_dict[C_hierarchy] += count
                D_dict[D_hierarchy] += count
        with open("../results/analysis/"+str(wildcards.sample)+"_count_B.tsv", "w") as outfile:
            for entry in B_dict.keys():
                outfile.write(entry+"\t"+str(B_dict[entry])+"\n")
        with open("../results/analysis/"+str(wildcards.sample)+"_count_C.tsv", "w") as outfile:
            for entry in C_dict.keys():
                outfile.write(entry + "\t" + str(C_dict[entry]) + "\n")
        with open("../results/analysis/"+str(wildcards.sample)+"_count_D.tsv", "w") as outfile:
            for entry in D_dict.keys():
                outfile.write(entry + "\t" + str(D_dict[entry]) + "\n")
        #shell("touch ../results/analysis/tmp.file")

rule combine_rules:
    input: expand("../results/analysis/{sample}_count_B.tsv", sample=SAMPLES3)
    output: temp("../results/analysis/tmp.file")
    shell: "touch ../results/analysis/tmp.file"

rule merge_hierarchy_counts:
    input: "../results/analysis/tmp.file"
    output: "/home/stephan/test/results/analysis/sec_paper/compare_meta_2.csv" #temp()
    script:"../prep_metabolsim_R_sec_paper.py"

rule order_hierarchy:
    input:
        file = "/home/stephan/test/results/analysis/sec_paper/compare_meta_2.csv",
        template = "../raw/modules_kegg.txt"
    output: "/home/stephan/test/results/analysis/sec_paper/compare_meta_3.csv"
    run:
        dict_set = {}
        with open(str(input.file),"r") as infile:
            for line in infile:
                name = line.split("\t")[0]
                dict_set[name] = line
        with open(str(input.template),"r") as template, open(output[0],"w") as outfile:
            outfile.write("0\t"+dict_set["pathway"])
            i = 1
            for line in template:
                if line.startswith("C"):
                    level_C = " ".join(line.split()[1:])
                    try:
                        outfile.write(str(i)+"\t"+dict_set[level_C])
                        del dict_set[level_C]
                        i += 1
                    except KeyError:
                        pass




rule kegg_modules:
    input:
        module= expand("../results/analysis/test/KEGG_web/{sample}_kegg_module.csv", sample=SAMPLES3)
    output:
        "../results/analysis/test.txt"
    run:
        sample_dict = {}
        #func_dict = {}
        dict_1 = {}
        for sample in input.module:
            #dict_1 = {}
            #dict_2 = {}
            #dict_3 = {}
            name = sample.split("/")[-1].split("_kegg_module.csv")[0]
            print(name)
            #create dict of modules per organism
            '''
            with open(sample,"r") as infile:
                for line in infile:
                    data = line.split(",")
                    if not data[0] =="":
                        now_1 = " ".join(data).strip()
                        #dict_1[now_1]=None
                    elif data[2] !="":
                        now_2 =" ".join(data).strip()
                    elif data[3] =="" and data[4] !="":
                        now_3 = " ".join(data[5:]).strip()
                        dict_1[data[4]] =(now_1,now_2,now_3)
            sample_dict[name]=dict_1
            '''
            # compare dict:
            with open(sample, "r") as infile:
                for line in infile:
                    data = line.split(",")
                    if not data[0] == "":
                        now_1 = " ".join(data).strip()
                        # dict_1[now_1]=None
                    elif data[2] != "":
                        now_2 = " ".join(data).strip()
                    elif data[3] == "" and data[4] != "":
                        now_3 = " ".join(data[5:]).strip()
                        dict_1[data[4]] = (now_1, now_2, now_3)
                        if not data[4] in sample_dict.keys():
                            sample_dict[data[4]]=[name]
                        else:
                            sample_dict[data[4]].append(name)

            #print(sample_dict)
        with open(output[0],"w") as outfile_1:
            for entry in sample_dict.keys():
                try:
                    outfile_1.write(entry+"\t"+"_".join(sample_dict[entry])+"\t"+"\t".join(dict_1[entry])+"\n")
                except KeyError as error:
                    print(error)

rule metabolism_compare:
    input:
        expand("../results/analysis/{sample}_metabolism_count_hierarchy_2.csv", sample=SAMPLES3)
    output:
        "test.txt"
    script:
        "R/compare_metabolism.R"



'''

'''