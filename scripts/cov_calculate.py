# calculate: coverage of assembly
# -coverage of assembly
# -percentage euk filter
# -N50

import numpy as np
import glob
assembly_fliles = "../results/new_spades/*/scaffolds.fasta.filter"

files = glob.glob(assembly_fliles)

def sum_up_fastq(fastq_file,my_sum2):
    with open(fastq_file, "r") as infile:
        switch = False
        for line in infile:
            if line.startswith("@"):
                switch = True
            elif switch:
                my_sum2 += len(line.strip())
                switch = False
    return my_sum2

#genome_dict = {}
with open("../results/assembly_coverage.tsv","w") as outfile:
    outfile.write("\t".join(("strain","coverage","%used_read_bases","N50","sum_contigs","number of contigs","\n")))
    for fasta_file in files:
        genome_name = fasta_file.split("/")[-2]
        my_sum = 0
        number_contigs = 0
        length_list = []
        with open(fasta_file,"r") as infile:
            for line in infile:
                line = line.strip()
                if not line.startswith(">"):
                    my_sum += len(line)
                    length_list.append(len(line))
                else:
                    number_contigs += 1

        a = np.array(length_list)
        unique_elements, counts_elements = np.unique(a, return_counts=True)
        # print("50 percentile:")
        # print(np.percentile(a,50))
        half_of_all_bases_length = sum(length_list) / 2
        length_array = len(unique_elements)
        current_length = 0
        for i in range(0, length_array - 1):
            current_length_tmp = current_length + unique_elements[i] * counts_elements[i]
            if current_length_tmp < half_of_all_bases_length:
                current_length = current_length_tmp
                N50 = unique_elements[i]
            elif current_length + unique_elements[i] < half_of_all_bases_length:
                current_length = current_length_tmp
                N50 = unique_elements[i]
            else:
                break

        # calculate coverage and percentage used
        my_sum2 = 0
        try:
            for i in range(1,3):
                fastq_file = "../results/bowtie/{}/merge.{}.fq".format(genome_name,i)
                my_sum2 = sum_up_fastq(fastq_file, my_sum2)
            coverage = float(my_sum2) / my_sum

            my_sum3 = 0
            original_fastq = glob.glob("../raw/illumina/{}*/*.fq".format(genome_name))
            for o_fastq in original_fastq:
                my_sum3 = sum_up_fastq(o_fastq, my_sum3)
            percentage_used_reads = float(my_sum2) / my_sum3
        except FileNotFoundError:
            coverage = "n.A."
            percentage_used_reads = 100

        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(genome_name,coverage,percentage_used_reads,N50,my_sum,number_contigs))
        print("{}\t{}\t{}\t{}\t{}\t{}\n".format(genome_name,coverage,percentage_used_reads,N50,my_sum,number_contigs))