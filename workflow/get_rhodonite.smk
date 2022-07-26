from snakemake.utils import min_version
import pandas as pd

min_version("6.0")

module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/v0.7-alpha/workflow/Snakefile"
    config:
        config

manifest = pd.read_csv(config["manifest"], sep = "\t")
config["samples"] = dict( zip( list(manifest["sample"]) ,  list(manifest["fasta"]) ) ) #list(manifest["fasta"]) ) )

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
