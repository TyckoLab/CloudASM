#!/usr/bin/env python

import os
import scipy.stats as stats
import numpy as np
import pandas as pd
import ndjson
from pandas.io.json import json_normalize



# The schema of the cpg genotype file is the following and the header is included
# snp_id
# ref_reads: number of reads on the REF
# alt_reads: number of reads on the ALT
# ref: array of 'methyl', which is the methyl % of each REF reaf
# alt: array of 'methyl', which is the methyl % of each ALT reaf

INPUT_FILE = os.environ['DMR']
OUTPUT_FILE = os.environ['DMR_PVALUE']


# load from file-like objects
with open(INPUT_FILE) as f:
    data = ndjson.load(f)

# Converte the JSON file in dataframe
df = json_normalize(data)

# Function to extract Wilcoxon p-value (5-digit rounding)
def wilcoxon_pvalue(row):
  _, pvalue =  stats.mannwhitneyu(json_normalize(row['ref']), \
                                json_normalize(row['alt']), \
                                )
  return round(pvalue,5)

# Create a column with the p-value
df['wilcoxon_pvalue'] = df.apply(wilcoxon_pvalue, axis = 1)

# Save to JSON
df.to_json(OUTPUT_FILE, orient = "records", lines = True)
