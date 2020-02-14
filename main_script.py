'''
this script identify and removes prokaryotic contigs out of metagenomic assemblies

input:
	fastq files with illumina & pacbio reads
	assembled scaffolds or contigs
functions:
	Assembly
	gene prediction
	KEGG annotation
'''
from Bio import SeqIO
import glob
import shutil
import os, sys
#import subprocess
#from multiprocessing import Process, Manager
from _collections import defaultdict
#from itertools import zip_longest

# set your samples
SAMPLES = ["1006","199hm", "AR4D6",  "Chromulina",  "Epipyxis",  "FU18K-A",  "JBC27",  "JBMS11",  "LO226KS",  "LO234KE",  "WA18K-M","DS","933-7"]  #list_of_spezies_names # "1006","199hm",  "AR4D6",  "Chromulina",  "Epipyxis",  "FU18K-A",  "JBC27",  "JBMS11",  "LO226KS",  "LO234KE",  "WA18K-M"
non_axenic_SAMPLES = ["1006","199hm", "AR4D6",  "Chromulina",  "Epipyxis",  "FU18K-A",  "JBC27",  "JBMS11",  "LO226KS",  "LO234KE",  "WA18K-M"] #list_of_spezies_names contaminated with bacteria
axenic_SAMPLES = ["DS","933-7","JBC07", "JBM10", "JBNZ41"]
SAMPLES3 = set(axenic_SAMPLES+SAMPLES)

Illumina_mRNA= "../raw/transcriptome/" #folder with transcriptome data
Illumina = "../raw/illumina/"

# get all info from KEGG files
#include: "scripts/get_kegg_info_2.py"
#include: "scripts/SnM_trim_uniprot.py"

rule all:
    input:
        #expand("{sample}_complete.csv", sample=SAMPLES),
        #expand("../results/kraken/{sample}_bin_statistic.tsv", sample=SAMPLES),
        #expand("../results/bowtie/{sample}/merge.2.fq", sample=SAMPLES),
        #expand("../results/bowtie/{sample}/euk.backmapped.1.fq", sample=SAMPLES),
        #expand("../results/MetaProb/{sample}/composition.txt", sample=SAMPLES),
        expand("../results/{sample}/_length_statistic_percentile.txt", sample=SAMPLES), #todo
        ##expand("../results/maxbin/{sample}/bin2.list", sample=SAMPLES),
        expand("../results/kraken/{sample}/pacbio_kraken.log", sample=SAMPLES),
        expand("../results/new_spades/{sample}/scaffolds.fasta.filter",  sample=SAMPLES),
        expand("../results/augustus/{sample}_augustus_protein_seq.fa",  sample=SAMPLES),
        expand("../results/smudge/{sample}_touch",  sample=SAMPLES),
        "../results/assembly_coverage.tsv",
        "../results/cd_hit/cluster_consensus.fa"


rule canu_assembly:
    input:
        size = "../raw/genome_size.csv"
    output:
        out = "../results/pacbio_sequence.txt"
    run:
        fastq_list = subprocess.getoutput("find /media/tintri/raw/pacbio -name *merge.fastq")
        with open(output.out,"w") as outfile:
            outfile.write(fastq_list)
            outfile.write("\n")

        genome_size = {}
        with open(input.size,"rb") as infile:
            for line in infile:
                try:
                    data = line.split("\t")
                    try:
                        size = int(data[4].split(".")[0])
                    except:
                        size = 100
                        species_name = data[0].replace(" ", "")
                        species_name = species_name.replace("|","")
                        #species_name = species_name.replace("-", "")
                    genome_size[species_name] = size
                except:
                    pass

        for tmp_sample in SAMPLES:
            for element in fastq_list.split("\n"):
                if tmp_sample in element:
                    try:
                        size = genome_size[tmp_sample]
                    except:
                        size = 100
                    if not os.path.exists("/media/tintri/results/canu/{}/sm_log.txt".format(tmp_sample)):
                        shell('echo "canu genomeSize={size}m correctedErrorRate=0.105 -p {tmp_sample} -d /media/tintri/results/canu/{tmp_sample} -pacbio-raw {element}" > /media/tintri/results/canu/{tmp_sample}_sm_log.txt')
                        shell("/media/tintri/software/canu/Linux-amd64/bin/canu genomeSize={size}m correctedErrorRate=0.105 -p {tmp_sample} -d /media/tintri/results/canu/{tmp_sample} -pacbio-raw {element}")

rule find_illumina_data:
    input:
        size = "../raw/genome_size.csv"
    output:
        temp("../results/illumina_files.txt")
    run:
        for tmp_sample in SAMPLES:
            shell("find /media/tintri/raw/illumina/{tmp_sample}*/*1.fq* >> {output[0]}")

rule SPAdes_assembly:
    input:
        #ill = expand("{illumina}{sample}.{replicate}.fq.gz", sample=SAMPLES, illumina = Illumina, replicate = [1, 2])
        #ill = "../raw/illumina/{dataset}.1.fq.gz",
        ill = "../results/illumina_files.txt",
        pacbio = "../results/pacbio_sequence.txt"
    output:
        "../results/spades/{sample}/scaffolds.fasta"
    #wildcard_constraints:
    #    dataset="{sample}-\w+/\w+"
    threads: 20
    run:
        tmp_sample = wildcards.sample
        if tmp_sample in non_axenic_SAMPLES:
            meta = "--meta "
        else:
            meta = ""
        pacbio_file_parameter =""
        with open(input.pacbio,"r") as infile:
            for line in infile:
                if tmp_sample in line:
                    pacbio_file_parameter =" --pacbio " + line.strip()
                    print(pacbio_file_parameter)
        if os.path.exists("/media/tintri/results/canu/{}/sm_log.txt".format(tmp_sample)):
            pacbio_file_parameter = " --untrusted-contigs " + "/media/tintri/results/canu/{}/{}.contigs.fasta".format(tmp_sample,tmp_sample)
            if meta == "--meta ":
                with open("../results/pacbio_sequence.txt","r") as pacbio_fastq:
                    for line in pacbio_fastq:
                        if tmp_sample in line:
                            pacbio_file_parameter = " --pacbio {}".format(line.strip())

        with open(input.ill, "r") as infile:
            for line in infile:
                if tmp_sample in line:
                    ill1 = line.strip()
                    ill2 = ill1.replace("1.fq","2.fq")
                    try:
                        shell("mkdir ../results/spades/tmp_{tmp_sample}")
                    except:
                        pass
                    shell("/media/tintri/software/spades/SPAdes-3.13.0-Linux/bin/spades.py -m 370 -1 {ill1} -2 {ill2} -t {threads}{pacbio_file_parameter} {meta}-o ../results/spades/tmp_{tmp_sample}")
                    shell("mv ../results/spades/tmp_{tmp_sample}/{{scaffolds.fasta,params.txt}} ../results/spades/{tmp_sample}/.")
                    try:
                        shell("mv ../results/spades/tmp_{tmp_sample}/warnings.log ../results/new_spades/{tmp_sample}/.")
                    except:
                        pass
                    shell("rm -r ../results/spades/tmp_{tmp_sample}/")

# filter contigs < 500 bps
rule filter_contigs:
    input:
        "../results/spades/{sample}/scaffolds.fasta"
    output:
        "../results/spades/{sample}/scaffolds.fasta.filter"
    script:
        "scripts/extract_high.py"

#toDo: implement tool
rule get_assembly_statistic:
    input:
        "../results/spades/{sample}/scaffolds.fasta.filter"
    output:
        "../results/{sample}/_length_statistic_percentile.txt"
    script:
        "scripts/biopython_len_fasta_modul.py" # find tool


rule MaxBin2_metagenomic_binning:
    input:
        ill = "../results/illumina_files_modify.txt",
        contigs = "../results/spades/{sample}/scaffolds.fasta.filter",
    output:
        "../results/maxbin/{sample}/bin2.list",
        #touch("{sample}maxbin_check.txt")
    threads: 18
    run:
        tmp_sample = wildcards.sample
        with open(input.ill, "r") as infile:
            for line in infile:
                if tmp_sample in line:
                    ill1 = line.strip()
                    ill2 = ill1.replace("1.fq", "2.fq")
        if tmp_sample in non_axenic_SAMPLES:
            # try:
            #     shell("gzip -d {ill1} {ill2}") # only needed for maxbin Version < 2.4
            #     ill1 = ill1.replace(".fq.gz", ".fq")
            #     ill2 = ill2.replace(".fq.gz", ".fq")
            # except:
            #     pass
            shell("/media/tintri/software/MaxBin-2.2.5/run_MaxBin.pl -contig {input.contigs} -reads {ill1} -reads2 {ill2} -thread {threads} -reassembly -out ../results/maxbin/{tmp_sample}/{tmp_sample}")
            shell("ls ../results/maxbin/{wildcards.sample}/{tmp_sample}.*.fasta > ../results/maxbin/{tmp_sample}/bin2.list")
        else:
            shell("touch {output[0]")
        # Maxbin creates bin, but fails to reassemble flawless -> pass task
#
rule build_Kraken_DB:
    input:
        "../raw/genome_size.csv"
    output:
        touch("../results/kraken/DB_tmp.txt")
    threads: 16
    run:
        #shell("kraken2-build --standard --db /media/tintri/kraken/2018_11_01_kraken_DB -threads {threads}")
        try:
            shell("mkdir /media/tintri/kraken/add_lib")
        except:
            pass
        for axenic in axenic_SAMPLES:
            try:
                with open("/media/tintri/results/spades/{}/scaffolds.fasta.filter".format(axenic),"r") as infile:
                    with open("/media/tintri/kraken/add_lib/{}.fa".format(axenic),"w") as outfile:
                        for line in infile:
                            if line.startswith(">"):
                                outfile.write(line.strip()+"|kraken:taxid|98651\n") # taxid from chrysophyceae
                            else:
                                outfile.write(line)
                try:
                    shell("kraken2-build -db /media/tintri/kraken/2018_11_01_kraken_DB -threads {threads} --add-to-library /media/tintri/kraken/add_lib/{axenic}.fa")
                except Exception as error:
                    print(error)
            except Exception as error:
                print(error)
        #shell("kraken2-build --build --db /media/tintri/kraken/2018_11_01_kraken_DB -threads {threads}")

rule modify_fastq:
    input:
        "../results/illumina_files.txt",
    output:
        temp("../results/illumina_files_modify.txt")
    threads: 18
    run:
        #todo: only 1.fq in list!!

        fastq =[]
        with open(input[0],"r") as infile:
            for line in infile:
                line = line.strip()
                if "fq.gz" in line:
                    for i in range(1, 3):
                        new_outfile1 = line.replace("1.fq.gz", "{}.fq.gz".format(i))
                        new_outfile2 = new_outfile1.replace("{}.fq.gz".format(i),"{}.fq".format(i))
                        shell("gzip -d -c {new_outfile1} > {new_outfile2}")
                        fastq.append(new_outfile2)
                else:
                    fastq.append(line)
                    fastq.append(fastq[0].replace("1.fq", "2.fq"))

        # for fastq_file in fastq:
        #     with open(fastq_file,"r") as infile, open(str(fastq_file).replace(".fq","")+".fq","w") as outfile:
        #         for line in infile:
        #             if line.startswith("@"):
        #                 new_id = line.split()
        #                 line = new_id[0]+"\n"
        #             outfile.write(line)
        shell("find /media/tintri/raw/illumina/*/*1.fq >> {output[0]}")


rule Kraken_pacbio_bin_classification:
    input:
        "../results/pacbio_sequence.txt",
        "../results/kraken/DB_tmp.txt"
    output:
        "../results/kraken/{sample}/pacbio_kraken.log"
    threads: 18
    run:
        tmp_sample = wildcards.sample
        if tmp_sample in non_axenic_SAMPLES:
            if os.path.exists("/media/tintri/results/canu/{}/sm_log.txt".format(tmp_sample)):
                pacbio= "/media/tintri/results/canu/{}/{}.contigs.fasta".format(tmp_sample, tmp_sample)
                kraken_DB = "/media/tintri/kraken/2018_11_01_kraken_DB/"
                shell("kraken2 --db {kraken_DB} --use-names --threads {threads} --classified-out ../results/kraken/{tmp_sample}/pacbio_{tmp_sample}.fq --unclassified-out ../results/kraken/{tmp_sample}/pacbio_{tmp_sample}_unclassified.fq {pacbio} >> {output[0]}")
                # WARNING: kraken2 fastq outputs for reads 1 & 2 are equal, taxid is wrong

                # repeat for unassembled
                pacbio = "/media/tintri/results/canu/{}/{}.unassembled.fasta".format(tmp_sample, tmp_sample)
                shell("kraken2 --db {kraken_DB} --use-names --threads {threads} --classified-out ../results/kraken/{tmp_sample}/pacbio_{tmp_sample}_unassembled.fa --unclassified-out ../results/kraken/{tmp_sample}/pacbio_{tmp_sample}_unassembled_unclassified.fa {pacbio} >> {output[0]}_unassembled")
            else:
                shell("touch {output[0]}")
        else:
            shell("touch {output[0]}")

rule Kraken_pacbio_bin_processing:
    input:
        log = "../results/kraken/{sample}/pacbio_kraken.log"
    output:
        new_pacbio = "../results/kraken/{sample}/euK_pacbio_{sample}.fa",
    threads: 18
    run:
        def walk_through(taxid, taxonomy, path,line):
            parent = taxonomy[taxid]
            if parent != taxid or parent != "1":
                if taxid in ["2","2759","2157","10239"]:
                    # save the trouble of full path calculation
                    return path, taxid
                else:
                    path = "_".join((parent, path))
                    try:
                        path, taxid = walk_through(parent, taxonomy, path,line)
                    except Exception as error:
                        print(error)
                        print(line)
                    return path, taxid
            return path, taxid

        tmp_sample = wildcards.sample
        if os.stat(input.log).st_size > 0:
            taxonomy = {}
            with open("/media/tintri/kraken/2018_11_01_kraken_DB/taxonomy/nodes.dmp", "r") as taxonomy_graph:
                for line in taxonomy_graph:
                    data = line.split("|")
                    sp1 = data[0].strip()
                    parent = data[1].strip()
                    taxonomy[sp1] = parent
            hits = {}
            hits["U"] = []
            #hit_length = defaultdict(int)
            for log_file in ["../results/kraken/{}/pacbio_kraken.log".format(tmp_sample),"../results/kraken/{}/pacbio_kraken.log_unassembled".format(tmp_sample)]:
                with open(log_file,"r") as kraken_classify:
                    for line in kraken_classify:
                        if line.startswith("C"):
                            data = line.split("\t")
                            species = data[2]
                            taxid = species.split(" ")[-1].strip(")")
                            contig_id = data[1]
                            #seq_length = int(data[3])
                            path = taxid
                            full_path, taxid = walk_through(taxid, taxonomy, path,line)
                            path = full_path.split("_")
                            suffix = "U"
                            for taxid in path:
                                if taxid == "2":
                                    suffix = "BAC"
                                    break
                                elif taxid == "2759":
                                    suffix = "EUK"
                                    break
                                elif taxid == "2157":
                                    suffix = "ARC"  # archea
                                    break
                                elif taxid == "10239":
                                    suffix = "VIR"  # virus
                                    break
                                elif taxid == "9606":
                                    suffix = "human"
                                    break
                            if suffix in hits.keys():
                                hits[suffix].append(contig_id)
                            else:
                                hits[suffix] = [contig_id]
                            #hit_length[suffix] += seq_length
                        elif line.startswith("U"):
                            contig_id = line.strip().split("\t")[1]
                            hits["U"].append(contig_id)
                filtered_file = ["../results/kraken/{}/pacbio_{}.fq".format(tmp_sample,tmp_sample),"../results/kraken/{}/euK_pacbio_{}.fa".format(tmp_sample,tmp_sample)]
                if log_file == "../results/kraken/{}/pacbio_kraken.log_unassembled".format(tmp_sample):
                    filtered_file = ["../results/kraken/{}/pacbio_{}_unassembled.fa".format(tmp_sample, tmp_sample),
                                     "../results/kraken/{}/euK_pacbio_{}_unassembled.fa".format(tmp_sample, tmp_sample)]
                with open(filtered_file[0],"r") as fasta_infile, open(filtered_file[1].format(tmp_sample,tmp_sample),"w") as pac_euk_out:
                    with open("../results/kraken/{}/human_pacbio_{}.fa".format(tmp_sample, tmp_sample), "w") as pac_human_out:
                        write_switch = False
                        hwrite_switch = False
                        for line in fasta_infile:
                            if line.startswith(">"):
                                contig_id = line.split()[0].strip(">")
                                try:
                                    if contig_id in hits["EUK"] or contig_id in hits["U"]:
                                        write_switch = True
                                        pac_euk_out.write(line)
                                except:
                                    if contig_id in hits["U"]:
                                        write_switch = True
                                        pac_euk_out.write(line)
                                    elif contig_id in hits["human"]:
                                        hwrite_switch = True
                                        pac_human_out.write(line)
                            elif write_switch:
                                pac_euk_out.write(line)
                                write_switch = False
                            elif hwrite_switch:
                                pac_human_out.write(line)
                                hwrite_switch = False
        else:
            shell("touch {output.new_pacbio}")



rule Kraken_bin_classification:
    input:
        "../results/maxbin/{sample}/bin2.list",
        "../results/kraken/DB_tmp.txt"
    output:
        #"../results/kraken/{sample}/list.txt",
        "../results/kraken/{sample}_bin_statistic.tsv",
        "../results/maxbin/{sample}/{sample}_bin_analysis.list",
        "../results/kraken/{sample}/{sample}.001.fasta.report"
    threads: 18
    run:
        def walk_through(taxid,taxonomy,path):
            parent = taxonomy[taxid]
            if parent != taxid or parent != "1":
                path = "_".join((parent,path))
                try:
                    path, taxid = walk_through(parent, taxonomy, path)
                except Exception as error:
                    print(error)
                return path, taxid
            return path, taxid

        kraken_DB = "/media/tintri/kraken/2018_11_01_kraken_DB/"
        infile = os.path.abspath(input[0])

        taxonomy = {}
        with open("/media/tintri/kraken/2018_11_01_kraken_DB/taxonomy/nodes.dmp", "r") as taxonomy_graph:
            for line in taxonomy_graph:
                data = line.split("|")
                sp1 = data[0].strip()
                parent = data[1].strip()
                taxonomy[sp1] = parent
        with open(infile,"r") as sample_list:
            for file in sample_list:
                filepath = file.strip()
                name = os.path.basename(file).strip()
                name_prefix = name.split(".")[0]
                shell("kraken2 --threads {threads} --db {kraken_DB} {filepath} --use-names > /tmp/{name_prefix}_result_kraken.tmp")
                #shell("kraken2 --db {kraken_DB} --report ../results/kraken/{name_prefix}/{name}.report /tmp/result_kraken.tmp   && ls ../results/kraken/{name_prefix}/*.report > {output}")
                #need perl > 5.26.2
                # python alternative:
                hits = {}
                hits["U"] =[]
                hit_length = defaultdict(int)
                #U = 0
                with open("/tmp/{}_result_kraken.tmp".format(name_prefix), "r") as kraken_classify:
                    for line in kraken_classify:
                        if line.startswith("C"):
                            data = line.split("\t")
                            species = data[2]
                            taxid = species.split(" ")[-1].strip(")")
                            contig_id = data[1]
                            seq_length = int(data[3])
                            path = taxid
                            full_path, taxid = walk_through(taxid, taxonomy, path)
                            path = full_path.split("_")
                            suffix = "U"
                            for taxid in path:
                                if taxid == "2":
                                    suffix = "BAC"
                                    break
                                elif taxid == "2759":
                                    suffix = "EUK"
                                    break
                                elif taxid == "2157":
                                    suffix = "ARC"  # archea
                                    break
                                elif taxid == "10239":
                                    suffix = "VIR"  # virus
                                    break
                            if suffix in hits.keys():
                                hits[suffix].append(contig_id)
                            else:
                                hits[suffix] = [contig_id]
                            hit_length[suffix] += seq_length
                        elif line.startswith("U"):
                            contig_id = line.strip().split("\t")[1]
                            hits["U"].append(contig_id)

                with open("/media/tintri/results/kraken/{}/{}.report".format(name_prefix, name), "w") as outfile:
                    with open(output[0],"a") as out_statistic:
                        out_statistic.write(name + "\t")
                        total_number = 0
                        total_length = 0
                        for suffix in hits.keys():
                            outfile.write("\t".join((suffix, str(len(hits[suffix])), str(hits[suffix]) + "\n")))
                            out_statistic.write("\t" +suffix + " " + str(len(hits[suffix]))+ " " + str(hit_length[suffix]))
                            total_number += len(hits[suffix])
                            total_length += hit_length[suffix]
                        total_dict ={}
                        for supergroup in ["EUK","BAC","U"]:
                            try:
                                percentage = float(len(hits[supergroup]))/total_number*100
                                percentage_len = float(hit_length[supergroup])/total_length*100
                                out_statistic.write("\n" +"percentage_{}".format(supergroup) + "\t{:10.4f}".format(percentage))
                                out_statistic.write("\t" + "percentage_len_{}".format(supergroup) + "\t{:10.4f}".format(percentage_len))
                                total_dict[supergroup] = [percentage,percentage_len]
                            except Exception as error:
                                print(error)
                        out_statistic.write("\n")

                # classify if bin is BAC or EUK
                with open(output[1], "a") as outfile2:
                    is_eukaryote = False
                    try:
                        if total_dict["EUK"][1]*2 > total_dict["BAC"][1]:
                            is_eukaryote = True
                        elif total_dict["U"][0] > total_dict["BAC"][0]*2:
                            is_eukaryote = True
                        # possibly total_dict["EUK"][1] < 1%, exclude possibly wrong assigned, but number is irrelevant
                    except Exception as error:
                        print(error)
                    outfile2.write("\t".join((name,str(is_eukaryote)+"\n")))

rule process_classified_bins:
    input:
        boolean = "../results/maxbin/{sample}/{sample}_bin_analysis.list",
        report = "../results/kraken/{sample}/{sample}.001.fasta.report"
    output:
        euk = "../results/maxbin/{sample}/euk.fa",
        bac = "../results/maxbin/{sample}/bac.fa"
    threads: 1
    run:
        import ast
        tmp_sample = wildcards.sample
        report_dict = {}
        with open(input.boolean,"r") as infile:
            for line in infile:
                data = line.strip().split("\t")
                fasta_file = data[0]
                with open("../results/kraken/{}/{}.report".format(tmp_sample,fasta_file),"r") as report_in:
                    report_number = fasta_file.split(".")[1]
                    report_dict[report_number] = defaultdict(dict)
                    for nline in report_in:
                        ndata = nline.strip().split("\t")
                        report_dict[report_number][ndata[0]] = ast.literal_eval(ndata[2])
                    report_dict[report_number]["is_euk"] = data[1]

        fasta_dict = {}
        with open(output.euk, "w") as euk_out, open(output.bac, "w") as bac_out:
            for fasta_file_number in report_dict.keys():
                fasta_dict[fasta_file_number] = defaultdict(str)
                with open("../results/maxbin/{}/{}.{}.fasta".format(tmp_sample, tmp_sample,fasta_file_number), "r") as fasta_infile:
                    for line in fasta_infile:
                        #file_switch = bac_out
                        if line.startswith(">"):
                            contig_id = line.strip().strip(">")
                        else:
                            fasta_dict[fasta_file_number][contig_id] += line
                #print(fasta_file_number)
                #print(len(fasta_dict[fasta_file_number]))

                outfile_dict = {"BAC":bac_out,"EUK":euk_out}
                if report_dict[fasta_file_number]["is_euk"] == "True":
                    outfile_dict["U"] = euk_out
                elif report_dict[fasta_file_number]["is_euk"] == "False":
                    outfile_dict["U"] = bac_out
                else:
                    print("something wrong:")
                    print(report_dict[fasta_file_number]["is_euk"])

                for myclass in report_dict[fasta_file_number].keys():
                    if myclass == "is_euk":
                        continue
                    if myclass not in outfile_dict.keys():
                        outfile_dict[myclass] = bac_out
                    for contig in report_dict[fasta_file_number][myclass]:
                        outfile_dict[myclass].write(">{}\n".format(contig))
                        outfile_dict[myclass].write(fasta_dict[fasta_file_number][contig])

rule build_bowtie_classified_bins:
    input:
        euk = "../results/maxbin/{sample}/euk.fa",
        bac = "../results/maxbin/{sample}/bac.fa"
    output:
        euk = "../results/maxbin/{sample}/euk.1.bt2",
        bac = "../results/maxbin/{sample}/bac.1.bt2"
    threads: 18
    run:
        suffix = ["euk","bac"]
        for i in range(0,2):
            new_input = input[i]
            current_suffix = suffix[i]
            shell('echo "bowtie2-build --threads {threads} -f {new_input} ../results/maxbin/{wildcards.sample}/{current_suffix}"')
            shell("bowtie2-build --threads {threads} -f {new_input} ../results/maxbin/{wildcards.sample}/{current_suffix}")

rule bowtie_classified_bins1:
    input:
        #euk  = "../results/maxbin/{sample}/euk.1.bt2",
        bac = "../results/maxbin/{sample}/bac.1.bt2",
        fastq ="../results/illumina_files_modify.txt"
    output:
        bac = "../results/bowtie/{sample}/bak.backmapped.1.fq",
        bac2 = "../results/bowtie/{sample}/bak.backmapped.2.fq",
        include = "../results/bowtie/{sample}/include.backmapped.1.fq",
    threads: 18
    run:
        tmp_sample = wildcards.sample
        bac_out = str(output.bac).replace(".1.fq",".fq")
        include_out = bac_out.replace("bak.backmapped","include.backmapped")
        file_list = []
        with open(input.fastq, "r") as infile:
            for line in infile:
                if tmp_sample in line:
                    file_list.append(line.strip())
                    file_list.append(file_list[0].replace("1.fq", "2.fq"))
                    break
        index = str(input.bac).split(".1.bt2")[0]
        shell("bowtie2 -N 1 --no-unal -x {index} -1 {file_list[0]} -2 {file_list[1]} -p {threads} -S ../results/maxbin/{tmp_sample}/bak.sam --al-conc {bac_out} --un-conc {include_out}")

rule bowtie_classified_bins2:
    input:
        euk  = "../results/maxbin/{sample}/euk.1.bt2",
        bac1 = "../results/bowtie/{sample}/bak.backmapped.1.fq",
        bac2 = "../results/bowtie/{sample}/bak.backmapped.2.fq",
    output:
        euk = "../results/bowtie/{sample}/euk.backmapped.1.fq",
        exclude = "../results/bowtie/{sample}/exclude.1.fq",
    threads: 18
    run:
        tmp_sample = wildcards.sample
        euk_out = str(output.euk).replace(".1.fq",".fq")
        exclude_out = str(output.exclude).replace(".1.fq",".fq")
        index = str(input.euk).split(".1.bt2")[0]
        shell("bowtie2 -N 1 --no-unal -x {index} -1 {input.bac1} -2 {input.bac2} -S ../results/maxbin/{tmp_sample}/euk.sam -p {threads} --al-conc {euk_out} --un-conc {exclude_out}")

rule merge_new_fastq:
    input:
        euk = "../results/bowtie/{sample}/euk.backmapped.1.fq",
        include = "../results/bowtie/{sample}/include.backmapped.1.fq",
    output:
        new1 = "../results/bowtie/{sample}/merge.1.fq",
        new2 = "../results/bowtie/{sample}/merge.2.fq",
    threads: 18
    run:
        # todo: parallel
        shell("cat {input.euk} {input.include} >> {output.new1}")
        #shell('echo "cat {input.euk} {input.include} >> {output.new1}"')
        euk2 = str(input.euk).replace("1.fq","2.fq")
        inc2 = str(input.include).replace("1.fq","2.fq")
        shell("cat {euk2} {inc2} >> {output.new2}")
        #shell('echo "cat {euk2} {inc2} >> {output.new2}"')

rule new_SPAdes_assembly:
    input:
        new1 = "../results/bowtie/{sample}/merge.1.fq",
        new2 = "../results/bowtie/{sample}/merge.2.fq",
        #pacbio = "../results/pacbio_sequence.txt"
        new_pacbio = "../results/kraken/{sample}/euK_pacbio_{sample}.fa",
    output:
        out = "../results/new_spades/{sample}/scaffolds.fasta"
    threads: 20
    run:
        tmp_sample = wildcards.sample
        #new_tmp_sample = "new_" + str(wildcards.sample)

        pacbio_file_parameter =""
        # replace if with check input.new_pacbio size:
        if os.stat(input.new_pacbio).st_size > 0:
            pacbio_file_parameter = " --untrusted-contigs ../results/kraken/{}/euK_pacbio_{}.fa --pacbio ../results/kraken/{}/euK_pacbio_{}_unassembled.fa".format(tmp_sample,tmp_sample,tmp_sample,tmp_sample)

        try:
            shell("mkdir ../results/spades/tmp_{tmp_sample}")
        except:
            pass
        shell("/media/tintri/software/spades/SPAdes-3.13.0-Linux/bin/spades.py -m 370 -1 {input.new1} -2 {input.new2} -t {threads}{pacbio_file_parameter} -o ../results/spades/tmp_{tmp_sample}")
        shell("mv ../results/spades/tmp_{tmp_sample}/{{scaffolds.fasta,params.txt}} ../results/new_spades/{tmp_sample}/.")
        try:
            shell("mv ../results/spades/tmp_{tmp_sample}/warnings.log ../results/new_spades/{tmp_sample}/.")
        except:
            pass
        shell("rm -r ../results/spades/tmp_{tmp_sample}/")

rule filter_contigs_again:
    input:
        "../results/new_spades/{sample}/scaffolds.fasta",
    output:
        "../results/new_spades/{sample}/scaffolds.fasta.filter",
    script:
        "scripts/extract_high.py"

rule calculate_coverage:
    input:
        expand("../results/new_spades/{sample}/scaffolds.fasta.filter",  sample=SAMPLES),
    output:
        "../results/assembly_coverage.tsv",
    script:
        "scripts/cov_calculate.py"

rule pre:
    input:
        "../results/illumina_files.txt",
        #todo: new trigger
    output:
        expand("../results/smudge/{sample}_k21.hist", sample=SAMPLES),
    threads: 16
    run:
        for tmp_sample in axenic_SAMPLES:
            if not os.path.exists("/media/tintri/results/smudge/{}_k21.hist".format(tmp_sample)):
                try:
                    #shell("ls /media/tintri/raw/illumina/{tmp_sample}*.fq.gz > FILES_{tmp_sample}")
                    #shell("/media/tintri/software/KMC/bin/kmc -k21 -m300 -ci1 -cs10000 -t14 @FILES_{tmp_sample} ../results/smudge/{tmp_sample}_kmer_counts /media/tintri/results/ploidy/tmp/")
                    shell("/media/tintri/software/KMC/bin/kmc_tools transform ../results/smudge/{tmp_sample}_kmer_counts histogram ../results/smudge/{tmp_sample}_k21.hist -cx10000")
                    #shell("gzip {tmp_sample}_kmer_counts.kmc_suf")
                except:
                    print("Did not find raw data:")
                    print(tmp_sample)
        for tmp_sample in non_axenic_SAMPLES:
            if not os.path.exists("/media/tintri/results/smudge/{}_k21.hist".format(tmp_sample)):
                shell("ls /media/tintri/results/bowtie/{tmp_sample}/merge.*.fq > FILES_{tmp_sample}")
                shell("/media/tintri/software/KMC/bin/kmc -k21 -m300 -ci1 -cs10000 -t14 @FILES_{tmp_sample} ../results/smudge/{tmp_sample}_kmer_counts /media/tintri/results/ploidy/tmp/")
                shell("/media/tintri/software/KMC/bin/kmc_tools transform ../results/smudge/{tmp_sample}_kmer_counts histogram ../results/smudge/{tmp_sample}_k21.hist -cx10000")
                #shell("gzip {tmp_sample}_kmer_counts.kmc_suf")

# rule smu:
#     input:
#         hist = "../results/smudge/{sample}_k21.hist",
#     output:
#         "../results/smudge/{sample}_kmer_pairs_coverages_2.tsv"
#     threads: 1
#     run:
#         tmp_sample = wildcards.sample
#         #shell("L=$(kmer_cov_cutoff.R ../results/smudge/{tmp_sample}_k21.hist L) && U=$(kmer_cov_cutoff.R ../results/smudge/{tmp_sample}_k21.hist U) && /media/tintri/software/KMC/bin/kmc_dump -ci$L -cx$U ../results/smudge/{tmp_sample}_kmer_counts ../results/smudge/{tmp_sample}_kmer_k21.dump")
#         shell("LA=5 && L=$(kmer_cov_cutoff.R ../results/smudge/{tmp_sample}_k21.hist L) && LCOV=$(($L-$LA)) && U=$(kmer_cov_cutoff.R ../results/smudge/{tmp_sample}_k21.hist U) && /media/tintri/software/KMC/bin/kmc_dump -ci$LCOV -cx$U ../results/smudge/{tmp_sample}_kmer_counts ../results/smudge/{tmp_sample}_kmer_k21.dump")
#         shell("hetkmers.py -k 21 -t {threads} -o ../results/smudge/{tmp_sample}_kmer_pairs < ../results/smudge/{tmp_sample}_kmer_k21.dump")

rule visualize:
    input:
        "../results/smudge/{sample}_kmer_pairs_coverages_2.tsv",
    output:
        touch("../results/smudge/{sample}_touch")
    threads: 1
    run:
        tmp_sample = wildcards.sample
        shell("smudgeplot.R -i ../results/smudge/{tmp_sample}_kmer_pairs_coverages_2.tsv -o ../results/smudge/{tmp_sample}")

rule gene_prediction:
    input:
        "../results/new_spades/{sample}/scaffolds.fasta.filter",
    output:
        "../results/augustus/{sample}_augustus.out",
    run:
        shell("augustus --species=arabidopsis {input[0]} --gff3=on --progress=true --singlestrand=true --AUGUSTUS_CONFIG_PATH=/media/tintri/software/Augustus/config/ --UTR=off > {output[0]}")
        #todo: remove comments from gff

rule gene_extraction:
    input:
        "../results/augustus/{sample}_augustus.out",
    output:
        "../results/augustus/{sample}_augustus_protein_seq.fa",
    threads: 1
    run:
        #protein_dict = defaultdict(str)
        with open(input[0],"r",encoding='utf-8') as infile,open(output[0],"w") as outfile:
            switch = False
            id = 0
            for line in infile:
                if line.startswith("# protein sequence = ["):
                    switch=True
                    id +=1
                    outfile.write(">g{}\n".format(id))
                elif line.startswith("# end gene"):
                    switch=False
                    outfile.write("\n")
                if switch:
                    seq = line.lstrip("# protein sequence = [").strip().rstrip("]")
                    outfile.write(seq)

# rule gene_annotation_interproscan:
#     input:
#         "../results/augustus/{sample}_augustus_protein_seq.fa",
#     output:
#         touch("../results/augustus/{sample}_interprotouch.txt")
#     threads: 5
#     run:
#         tmp_sample = wildcards.sample
#         folder = os.path.abspath("../results/augustus/{}_augustus_protein_seq.fa".format(tmp_sample))
#         try:
#             shell("mkdir ../results/interproscan")
#         except:
#             pass
#         shell("/media/tintri/software/interproscan-5.32-71.0/interproscan.sh -d /media/tintri/results/interproscan -dp --formats TSV -goterms -i {folder} -ms 25 -pa")

# todo: rule Interproscan_classify.py

rule augustus2dna:
    input:
        genome = "../results/new_spades/{sample}/scaffolds.fasta.filter",
        augustus_file = "../results/augustus/{sample}_augustus.out"
    output:
        augustus_dna = "../results/augustus/{sample}_augustus_CDS.fa",
    threads: 1
    run:
        contig_dict={}
        gene_dict=defaultdict(str)
        for seq_record in SeqIO.parse(input.genome, "fasta"):
            contig_dict[seq_record.id]=str(seq_record.seq)

        with open(input.augustus_file, "r", encoding="utf-8") as infile: #,
            for line in infile:
                if not line.startswith("#"):
                    data = line.split("\t")
                    if data[2] == "transcript":
                        gene = data[-1].strip().split("Parent=")[1]

                    elif data[2] == "CDS":
                        start = int(data[3]) - 1
                        end = int(data[4]) -1
                        gene_dict[gene] += contig_dict[data[0]][start:end]

        with open(output.augustus_dna, "w") as outfile:
            for gene in gene_dict.keys():
                outfile.write(">{}_{}\n".format(wildcards.sample,gene))
                outfile.write(gene_dict[gene] + "\n")

rule pre_cd_hit:
    input:
        augustus_file = expand("../results/augustus/{sample}_augustus_CDS.fa", sample=SAMPLES3),
    output:
        pool = "../results/cd_hit/pool.fa",
    threads: 1
    run:
        shell("echo {input.augustus_file}")
        shell("cat {input.augustus_file} > {output.pool}")
        # manually added JBM10,JBC07 and JBNZ41


rule cd_hit:
    input:
        pool = "../results/cd_hit/pool.fa",
    output:
        cluster = "../results/cd_hit/cd_hit-cluster3.clstr",
    threads: 18
    run:
        shell("cd-hit -i {input.pool} -o {output.cluster} -M 11000 -T {threads}")

rule cd_hit_consensus:
    input:
        cluster = "../results/cd_hit/cd_hit-cluster3.clstr",
        ref = "../results/cd_hit/pool.fa",
    output:
        "../results/cd_hit/cluster_consensus.fa",
    threads: 5
    run:
        shell("cdhit-cluster-consensus {input.cluster} {input.ref} {output[0]}")
# todo: minimap
# KEGG analysis
include: "scripts/SnM_analysis_5.py"