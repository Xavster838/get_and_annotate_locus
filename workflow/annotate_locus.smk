from snakemake.utils import min_version
import pandas as pd

min_version("6.0")

manifest = pd.read_csv(config["manifest"], sep = "\t")
manifest = manifest.set_index(["sample" , "hap"] , drop = False)

def get_fasta(wc):
    '''return fasta from manifest'''
    return manifest.loc[wc.sm, wc.h]['fasta']

rule _rename_fasta:
    '''rename sample contigs to only contain: sample_hap_tig1,2,... and merge to single combined fasta'''
    input:
      samp_fasta = get_fasta
    output:
      old_new_name_map = "rename_fastas/tigName_maps/{sm}_{h}_tig_name_changes.tbl" ,
      new_fa = "rename_fastas/fastas/{sm}_{h}.fasta",
      new_fai = "rename_fastas/fastas/{sm}_{h}.fasta.fai"
    conda: 
      "envs/env.yml"
    script:
      "scripts/rename_fasta.py"

rule _get_rhodonite_manifest:
    '''generate manifest needed to structure config for running dupmasker and repeat masker in mitchell's rhodonite snakemake'''
    input:
      old_new_name_maps = expand("rename_fastas/tigName_maps/{samp}_tig_name_changes.tbl" , samp = ["_".join(x) for x in list( zip(manifest["sample"], manifest["hap"]) )] ) #rules._rename_fasta.output.old_new_name_map
    output:
      dup_man = "Dupmasker/dup_manifest.tbl"
    params:
      sms = manifest["sample"] + "_" + manifest["hap"]
    run:
      new_name_map_df = pd.concat( [pd.read_csv(x, sep = "\t") for x in input.old_new_name_maps ] , axis = 0)
      out_df = pd.DataFrame(columns = ["sample", "fasta"])
      for i,sm in enumerate(params.sms):
          fasta = [f for f in new_name_map_df['new_fasta'] if sm in f]
          fasta = list(set(fasta))
          assert len(fasta) == 1 , f"rule _get_rhodonite_manifest failed for {sm} sample.  got ether zero or more than one fasta with {sm} in name: \n{row}"
          out_df.loc[i] = [sm ,fasta[0] ]
      out_df.to_csv(output.dup_man , sep = "\t", index = False, header = True)
        
      
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
# rule all:
#     input:
#         expand(rules.Rhodonite_DupMasker.output, sample=config["samples"].keys())

