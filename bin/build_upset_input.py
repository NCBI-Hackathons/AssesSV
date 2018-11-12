from glob import glob
import pandas as pd
import os
import gzip
import sys

path_to_vcf_files = sys.argv[1]
true_variants = sys.argv[2]
vcf_files = glob(path_to_vcf_files + "/giab*.txt")

all_variants = []
summary = {}

## Build master list of True Variants
with gzip.open(true_variants, 'rb') as f:
    for line in f:
        if line.startswith(b'#'):
        	continue
        all_variants.append(line.split()[2].decode("utf-8"))

## Iterate through each of the condition VCFs
for vcf in vcf_files:
    sample_variants = []
    with open(vcf, 'r') as v:
        sample_name = os.path.splitext(os.path.basename(vcf))[0]
        ## Build list of varints which exist in this condition
        for line in v:
            if line.startswith('#'):
                continue
            sample_variants.append(line.split()[2])
    ## Check and score variants that exist in this condition
    for variant in all_variants:
        if variant not in summary: summary[variant] = {}
        if variant in sample_variants:
            summary[variant][sample_name] = 1
        else:
            summary[variant][sample_name] = 0

## Convert summary dict to data frame and write to file
pd.DataFrame(summary).transpose().to_csv('upset_input.txt', sep=';', index_label='SV')
