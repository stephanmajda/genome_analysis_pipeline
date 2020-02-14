from Bio import SeqIO
#import collections
import numpy as np

sample =snakemake.sample
print(sample)
length_list = []
for seq_record in SeqIO.parse(snakemake.input[0], "fasta"):
#for seq_record in SeqIO.parse("/home/stephan/test/results/spades/"+sample+"/scaffolds.fasta", "fasta"):
	length_list.append(len(seq_record))

a = np.array(length_list)
unique_elements, counts_elements = np.unique(a, return_counts=True)
with open("../results/"+sample+"/_length_statistic_percentile.txt","w") as outfile_1:
	outfile_1.write("50,75,95 percentile:")
	outfile_1.write(np.percentile(a,50))
	outfile_1.write(np.percentile(a,75))
	outfile_1.write(np.percentile(a,95))

length_array = len(unique_elements)
with open("../results/"+sample+"/contig_length_distribution_statistic.csv","w") as outfile:
	for i in range(0,length_array -1):
		outfile.write("\t".join([str(unique_elements[i]),str(counts_elements[i])+"\n"]))

