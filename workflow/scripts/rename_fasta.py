import os
import pandas as pd
from Bio import SeqIO

#!/usr/bin/env python
# Author: Xavi Guitart
#   change names in fasta. output new named fasta and a table describing the contig/read name changes.
name_change_df = pd.DataFrame(columns = ['original_fasta', 'original_name', 'new_fasta' , 'new_name'])
with open(snakemake.input.samp_fasta) as original_file:
    records = [ i for i in SeqIO.parse(original_file, 'fasta') ] #generate list of fasta sequences
    for i,seq in enumerate(records):
        name_change_df.loc[i] = [ snakemake.input.samp_fasta , seq.id, snakemake.output.new_fa ,f"{snakemake.wildcards.sm}__{snakemake.wildcards.h}__{i}"]
        seq.id = name_change_df.loc[i, "new_name"]
        seq.description = "" #description added to end of header. was old name previously.

    if(i > 0): # if more than one record then append to existing fasta. Otherwise if just one record write new file.
        with open(snakemake.output.new_fa, "a") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
    else:
        with open(snakemake.output.new_fa, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

os.system(f"samtools faidx {snakemake.output.new_fa}")
name_change_df.to_csv( snakemake.output.old_new_name_map , sep = "\t", header = True, index = False)
