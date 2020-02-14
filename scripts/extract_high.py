from Bio import SeqIO
#import argparse
import sys
import numpy as np

#print(sys.argv[1])
#parser = argparse.ArgumentParser()
#parser.add_argument("sample", help="fasta file to extract high quality contigs",type=str)
sample =  str(snakemake.input[0])
#sample = sys.argv[1] #vars(parser.parse_args())
print(sample)

length_list = []
with open(sample+".filter","w") as outfile:
    for seq_record in SeqIO.parse(sample, "fasta"):
        _split = str(seq_record.id).split("_")
        length = int(_split[3])
        cov = float(_split[5])
        if length > 500 and cov > 1:
            length_list.append(len(seq_record))
            outfile.write(">"+seq_record.id+"\n")
            outfile.write(str(seq_record.seq)+"\n")

a = np.array(length_list)
unique_elements, counts_elements = np.unique(a, return_counts=True)
print("50,75,95 percentile:")
print(np.percentile(a,50))
print(np.percentile(a,75))
print(np.percentile(a,95))

length_array = len(unique_elements)
with open(sample+".filter.count","w") as outfile:
	for i in range(0,length_array -1):
		outfile.write("\t".join([str(unique_elements[i]),str(counts_elements[i])+"\n"]))