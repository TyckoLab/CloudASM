#!/usr/bin/env python

import os
import scipy.stats as stats
import numpy as np
import pandas as pd
import ndjson
from pandas.io.json import json_normalize



# The schema of the cpg genotype file is the following and the header is included
# snp_id
# chr
# dmr_inf
# dmr_sup 
# ref_reads
# alt_reads
# effect
# ref: array of 'methyl', which is the methyl % of each REF reaf
# alt: array of 'methyl', which is the methyl % of each ALT reaf
# nb_cpg
# nb_sig_cpg
# cpg: array of pos, effect, fisher_pvalue, ref_cov, and alt_cov

INPUT_FILE = os.environ['DMR']
OUTPUT_FILE = os.environ['DMR_PVALUE']


# load from file-like objects
with open(INPUT_FILE) as f:
    data = ndjson.load(f)

# Converte the JSON file in dataframe
df = json_normalize(data)


################################## Calculate p-value of DMR

# Function to extract Wilcoxon p-value (5-digit rounding)
def wilcoxon_pvalue(row):
    try:
        _, pvalue = stats.mannwhitneyu(json_normalize(row['ref']), json_normalize(row['alt']))
        return round(pvalue,5)
    # If the ref and alt datasets are equal or one is included in the other one:
    except ValueError:
        return 1


# Create a column with the p-value
df['wilcoxon_pvalue'] = df.apply(wilcoxon_pvalue, axis = 1)

################################## Calculate number of significant consecutive CpGs in the same direction.

def consecutive_cpg(row):
  if int(row['nb_sig_cpg']) > 1 :
      flat_cpg = json_normalize(row['cpg'])
      found = 0
      for index, row in flat_cpg.iterrows():
          if index > 0:
              if flat_cpg.iloc[index-1].fisher_pvalue < 0.05 and row.fisher_pvalue < 0.05 and np.sign(flat_cpg.iloc[index-1].effect) == np.sign(row.effect):
                  found = found + 1
      return found
  else:
      return 0

# Create a column with the number of consecutive CpGs that have significant ASM in the same direction
df['nb_consecutive_asm'] = df.apply(consecutive_cpg, axis = 1)

################################## Save file in JSON format

# Save to JSON
df.to_json(OUTPUT_FILE, orient = "records", lines = True)