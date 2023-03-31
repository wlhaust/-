# -*- coding: utf-8 -*-
import os          
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
import argparse
import sys

'''
python get_seq_from_proteinortho.py -i myproject.proteinortho.tsv -o output2 -p 0.3 -c 1
'''

#解析参数
#参数分为必须参数(required)和可选参数(additional)
parser = argparse.ArgumentParser(description="Options for get_seq_from_proteinortho.py", 
                                 add_help=True)
required = parser.add_argument_group("Required arguments")
required.add_argument('-i', '--infile', action="store", metavar='\b', type=str, 
                      required=True, default="myproject.proteinortho.tsv", 
                      help="Name of the hyde proteinortho file")  
required.add_argument('-o', '--output', action="store", metavar='\b', type=str, 
                      required=True, default="output2", help="Name of the output folder")
required.add_argument('-p', '--missing', action="store", metavar='\b', 
                      type=float, required=True, default=0.3, 
                      help="max allow missing species proportion")  
required.add_argument('-c', '--copy', action="store", metavar='\b', type=int, 
                      required=True, default=1, help="max allow low copy num")

args                                        = parser.parse_args()
proteinortho_output_file                    = args.infile
output_seq_dir                              = args.output
max_allow_missing_species_proportion        = args.missing
max_allow_low_copy_num                      = args.copy



# make output dir
try:
    os.makedirs(os.getcwd() + os.sep + output_seq_dir)
except:
    print("The output directory \"" + output_seq_dir + "\" already exists, script stop.")
    sys.exit()

# Import the FASTA sequence into memory
file_dict = {}  
for each in os.listdir(os.getcwd()):        
    if ".pep" in each:
        file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
    elif ".cds" in each:
        file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
    else:
        pass

# get sequences from ".cds" or ".pep" files.            
def get_nucleotide(IDs): 
    return str(file_dict[IDs.split("++")[0][:-3] + "cds"][IDs.split("++")[1]].seq)     
def get_protein(IDs): 
    return str(file_dict[IDs.split("++")[0][:-3] + "pep"][IDs.split("++")[1]].seq)    

# main function
def get_seq():
    protein_table = pd.read_csv(proteinortho_output_file, sep="\t")
    species_num = protein_table.shape[1] - 3
    min_species_num = species_num * (1 - max_allow_missing_species_proportion)
    print(min_species_num)
    columns_values = protein_table.columns.values
    i = 0
    for index, row in protein_table.iterrows():
        i = i + 1
        if row.iloc[0] < min_species_num:
            pass
        else:
            nucleotide_fasta = open(os.getcwd() + os.sep + output_seq_dir + os.sep + "ortho" + str(i) + "_cds.fasta", "a")
            protein_fasta = open(os.getcwd() + os.sep + output_seq_dir + os.sep + "ortho" + str(i) + "_pep.fasta", "a")
            for n in range(3, protein_table.shape[1]):
                if row.iloc[n] != "*":
                    if "," not in row.iloc[n]:
                        nucleotide_fasta.write(">" + columns_values[n] + "\n")
                        protein_fasta.write(">" + columns_values[n] + "\n")
                        nucleotide_fasta.write(get_nucleotide(columns_values[n] + "++" + row.iloc[n]) + "\n")
                        protein_fasta.write(get_protein(columns_values[n] + "++" + row.iloc[n]) + "\n")
                    elif len(row.iloc[n].split(",")) > max_allow_low_copy_num:
                        nucleotide_fasta.close()
                        protein_fasta.close()
                        os.remove(os.getcwd() + os.sep + output_seq_dir + os.sep + "ortho" + str(i) + "_cds.fasta")
                        os.remove(os.getcwd() + os.sep + output_seq_dir + os.sep + "ortho" + str(i) + "_pep.fasta")
                        break
                    else:
                        a = "a"
                        for each_seq in row.iloc[n].split(","):
                            seq = get_nucleotide(columns_values[n] + "++" + each_seq)
                            if len(seq) > len(a):
                                a = seq
                                b = each_seq
                        nucleotide_fasta.write(">" + columns_values[n] + "\n")
                        protein_fasta.write(">" + columns_values[n] + "\n")
                        nucleotide_fasta.write(get_nucleotide(columns_values[n] + "++" + b) + "\n")
                        protein_fasta.write(get_protein(columns_values[n] + "++" + b) + "\n")
                else:
                    pass

get_seq()