from snakemake.utils import min_version
import pandas as pd

min_version("6.0")

module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/v0.7-alpha/workflow/Snakefile"
    config:
        config

manifest = pd.read_csv(config["manifest"], sep = "\t")
manifest = manifest.set_index(["sample" , "hap"] , drop = False)

config["samples"] = dict( zip( list(manifest["sample"] + "__" + manifest["hap"] ) ,  list( "rename_fastas/fastas/" + manifest["sample"] + "_" + manifest["hap"] + ".fasta" ) ) ) #list(manifest["fasta"]) ) )

def get_fasta(wc):
    '''return fasta from manifest'''
    return manifest.loc[wc.sm, wc.h]['fasta']

rule _rename_fasta:
    '''rename sample contigs to only contain: sample_hap_tig1,2,... and merge to single combined fasta'''
    input:
      samp_fasta = get_fasta
    output:
      old_new_name_map = "rename_fastas/tigName_maps/{sm}_{h}_tig_name_changes.tbl" ,
      new_fa = temp("rename_fastas/fastas/{sm}_{h}.fasta"),
    conda: 
      "envs/env.yml"
    script:
      "scripts/rename_fasta.py"

# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*

# use RepeatMasker rules from Rhodonite
rule only_dup:
    input:
        expand(rules.Rhodonite_DupMasker.output, sample = config["samples"].keys() )

# rule get_locus_annotation:
#     '''given a sequence fasta, align and find places that locus maps to into he query haplotypes (for genes, duplicons, etc.).'''
#     input:
#         locus_fa = config["locus_fa"] ,
#         hap_fa = get_fasta,
#     output:
#         bed = "{sm}_{h}_mappings.bed"
#     shell:"""
# minimap2 -ax asm20 --secondary=yes -p 0.3 -N 10000 --eqx -r 500 -K 500M {input.hap_fa} {input.locus_fa} | \
#     samtools view -b - | samtools sort | \
#     bedtools bamtobed -i - | awk 'BEGIN{{OFS="\t"}}{{print $0, "{wildcards.r}__{wildcards.sm}_{wildcards.h}__" NR }}' > {output.bed} || touch {output.bed}
# """

rule all:
    input:
        expand(rules.Rhodonite_DupMasker.output, sample=config["samples"].keys())

