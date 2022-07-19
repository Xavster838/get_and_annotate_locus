import pandas as pd

manifest = pd.read_csv(config.manifest, sep = "\t")
manifest = manifest.set_index(["sample" , "hap"] , drop = False)

def get_fasta(wc):
    '''given wc region, sample, hap: return fasta of query from manifest.'''
    return manifest.loc[[wc.sm , wc.h]['fasta']

####################                                Rhodonite: RepeatMasker and DupMasker
module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/master/workflow/Snakefile"
    config:
        config = {
            "samples" : dict( zip( list(manifest["sample"] + "__" + manifest["hap"] ) ,  list(manifest["fasta"]) ) )
        }


# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*

# use RepeatMasker rules from Rhodonite
use rule RepeatMasker from Rhodonite as Rhodonite_RepeatMasker with:
    output:
        # You rename the output to anything you want
        # but maintain the order and keep "{sample}" in the name and the (.gz).
        # Every rule will have the standard output of the program
        # and a bed output that is sorted in the order of the ref.
        out="Mask/{sample}.rm.out",
        bed="Mask/{sample}.rm.bed.gz",


# use DupMasker rules from Rhodonite
use rule DupMasker from Rhodonite as Rhodonite_DupMasker with:
    output:
        extra="Mask/{sample}.duplicons.extra",
        bed="Mask/{sample}.duplicons.bed.gz",
################################################################################                                       

rule get_locus_annotation:
    '''given a sequence fasta, align and find places that locus maps to into he query haplotypes (for genes, duplicons, etc.).'''
    input:
        locus_fa = ,
        hap_fa = get_fasta,
    output:
        bed = "{sm}_{h}_mappings.bed"
    shell:"""
minimap2 -ax asm20 --secondary=yes -p 0.3 -N 10000 --eqx -r 500 -K 500M {input.hap_fa} {input.locus_fa} | \
    samtools view -b - | samtools sort | \
    bedtools bamtobed -i - | awk 'BEGIN{{OFS="\t"}}{{print $0, "{wildcards.r}__{wildcards.sm}_{wildcards.h}__" NR }}' > {output.bed} || touch {output.bed}
"""