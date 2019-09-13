#!/usr/bin/env python

import os
import scipy.stats as stats
import numpy as np
import pandas as pd
import json

# The schema of the cpg genotype file is the following and the header is included
# snp_id
# ref_reads: number of reads on the REF
# alt_reads: number of reads on the ALT
# ref: array of 'methyl', which is the methyl % of each REF reaf
# alt: array of 'methyl', which is the methyl % of each ALT reaf


CPG_GENOTYPE = os.environ['CPG_GENOTYPE']
OUTPUT_FILE = os.environ['CPG_ASM']

import newlinejson as nlj

data =  nlj.open('gm12878_snp_for_dmr.json')

with open('gm12878_snp_for_dmr.json') as f:
    print(f.read()))

with nlj.open('sample-data/dictionaries.json') as src, \
    for line in 5:
        stc.write(line)

 as smp, \
        with nlj.open('out.json', 'w') as dst:
    for line in src:
        dst.write(line)

with open('out.json') as f:
    print(f.read()))

# Import CpG file into Python
df = pd.read_csv(CPG_GENOTYPE) 

# Function to extract Fisher p-value (5-digit rounding)
def fisher_pvalue(row):
  _, pvalue =  stats.fisher_exact([[row['ref_meth'], \
                                row['alt_meth']], \
                              [row['ref_cov']-row['ref_meth'], \
                                row['alt_cov']-row['alt_meth'] \
                            ]])
  return round(pvalue,5)

# Create a column with the Fisher p-value
df['fisher_pvalue'] = df.apply(fisher_pvalue, axis = 1)

# Save to CSV
df.to_csv (OUTPUT_FILE, index = None, header = False)




#######

stats.mannwhitneyu(data1, data2


bashCommand = "cwm --rdf test.rdf --ntriples > test.nt"
import subprocess
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
