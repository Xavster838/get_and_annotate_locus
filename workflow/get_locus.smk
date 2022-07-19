rule index:
    input:
        fasta=get_fasta,
    output:
        mmi=temp("temp/index/{sm}.{h}.mmi"),
    conda:
        env.yaml
    threads: 8
    shell:"minimap2 {config.MM_OPTS} -t {threads} -d {output.mmi} {input.fasta}"

rule ref_rgn_query:
    input:
        ref = config.ref,
        fai = f"{config.ref}.fai",
        bed = config.bed
    output:
        ref_locus_fa = temp("temp/reference_locus.fasta")
        ref_locus_fai = temp("temp/reference_locus.fasta.fai")
    conda:
        env.yaml
    shell:"""
bedtools getfasta -name+ -bed <(bedtools flank -b {config.flank} -g {input.fai} -i {input.bed}) -fi {input.ref} > {output.ref_locus_fa}
samtools faidx {output.ref_locus_fa}
"""

# rule get_rgn_paf:
#     input:
#         q_ref = rules.ref_rgn_query.output.ref_locus_fa,
#         hap_fa = get_fasta 
#     output:
#         paf = temp("temp/{sm}.{h}.paf")
#     conda:
#         env.yaml
#     threads: 16

# rule get_rgn:
#     input:
#         paf = rules.get_rgn_paf.output.paf,
#     output:
#         fa = "{r}/{r}__{s}_{h}.fasta"
#         fai = "{r}/{r}__{s}_{h}.fasta.fai"
#     conda:
#         env.yaml
#     threads:
#         4
