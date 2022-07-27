from snakemake.utils import min_version
import pandas as pd

min_version("6.0")

module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/v0.7-alpha/workflow/Snakefile"
    config:
        config

if( "samples" in config ):
    continue
elif("dup_manifest" in config):
    dup_manifest = pd.read_csv(config["dup_manifest"], sep = "\t")
    config["samples"] = dict( zip( list(dup_manifest["sample"]) ,  list(dup_manifest["fasta"]) ) )

# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*

# use DupMasker rules from Rhodonite
use rule DupMasker from Rhodonite as Rhodonite_DupMasker with:
    output:
        extra="Dupmasker/{sample}.duplicons.extra",
        bed="Dupmasker/{sample}.duplicons.bed.gz",

# use RepeatMasker rules from Rhodonite
rule only_dup:
    input:
        expand(rules.Rhodonite_DupMasker.output, sample = config["samples"].keys() )


