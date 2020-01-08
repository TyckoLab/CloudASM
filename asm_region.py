#!/usr/bin/env python

import os
import scipy.stats as stats
import numpy as np
import pandas as pd
import ndjson
from pandas.io.json import json_normalize
import statsmodels.stats.multitest as mt


# The schema of the cpg genotype file is the following and the header is included
# snp_id
# chr
# asm_region_inf
# asm_region_sup 
# ref_reads
# alt_reads
# effect
# ref: array of 'methyl', which is the methyl % of each REF reaf
# alt: array of 'methyl', which is the methyl % of each ALT reaf
# nb_cpg
# nb_sig_cpg
# cpg: array of pos, effect, fisher_pvalue, ref_cov, and alt_cov

INPUT_FILE = os.environ['ASM_REGION']
OUTPUT_FILE = os.environ['ASM_REGION_PVALUE']
P_VALUE = float(os.environ['P_VALUE'])
BH_THRESHOLD = float(os.environ['BH_THRESHOLD'])

# load from file-like objects
with open(INPUT_FILE) as f:
    data = ndjson.load(f)

# Converte the JSON file in dataframe
df = json_normalize(data)


################################## Calculate p-value of asm_region

# Function to extract Wilcoxon p-value (5-digit rounding)
def wilcoxon_pvalue(row):
    try:
        _, pvalue = stats.mannwhitneyu(json_normalize(row['ref']), json_normalize(row['alt']), alternative = "two-sided")
        return round(pvalue,5)
    # If the ref and alt datasets are equal or one is included in the other one:
    except ValueError:
        return 1

# Create a column with the p-value
df['wilcoxon_pvalue'] = df.apply(wilcoxon_pvalue, axis = 1)


################################## Calculate p-value corrected for multiple testing using Benjamini–Hochberg

df['wilcoxon_corr_pvalue'] = mt.multipletests(df['wilcoxon_pvalue'], alpha = BH_THRESHOLD, method = 'fdr_bh')[1]
df['wilcoxon_corr_pvalue'] = df['wilcoxon_corr_pvalue'].round(5)

################################## Calculate number of significant consecutive CpGs in the same direction.

# Find consecutive significant ASM CpGs that are negative
def consecutive_neg_cpg(row):
  if int(row['nb_sig_cpg']) > 1 :
      flat_cpg = json_normalize(row['cpg'])
      max_nb_consec = 0
      current_nb_consec = 0
      for index, row in flat_cpg.iterrows():
          if (index > 0):
              if (flat_cpg.iloc[index-1].fisher_pvalue < P_VALUE and 
                  row.fisher_pvalue < P_VALUE and 
                  np.sign(flat_cpg.iloc[index-1].effect) == -1 and 
                  np.sign(row.effect) == -1):
                  if (current_nb_consec == 0): 
                      current_nb_consec = 2
                  else:
                      current_nb_consec = current_nb_consec + 1
                  max_nb_consec = max(max_nb_consec, current_nb_consec)
              else: 
                  current_nb_consec = 0
      return max_nb_consec
  else:
      return 0

# Find consecutive significant ASM CpGs that are negative
def consecutive_pos_cpg(row):
  if int(row['nb_sig_cpg']) > 1 :
      flat_cpg = json_normalize(row['cpg'])
      max_nb_consec = 0
      current_nb_consec = 0
      for index, row in flat_cpg.iterrows():
          if (index > 0):
              if (flat_cpg.iloc[index-1].fisher_pvalue < P_VALUE and 
                  row.fisher_pvalue < P_VALUE and 
                  np.sign(flat_cpg.iloc[index-1].effect) == 1 and 
                  np.sign(row.effect) == 1):
                  if (current_nb_consec == 0): 
                      current_nb_consec = 2
                  else:
                      current_nb_consec = current_nb_consec + 1
                  max_nb_consec = max(max_nb_consec, current_nb_consec)
              else: 
                  current_nb_consec = 0
      return max_nb_consec
  else:
      return 0

# Create a column with the number of consecutive CpGs that have significant ASM in the same direction
df['nb_consec_pos_sig_asm'] = df.apply(consecutive_pos_cpg, axis = 1)
df['nb_consec_neg_sig_asm'] = df.apply(consecutive_neg_cpg, axis = 1)

# Sort the dataframe per the chromosome column (pushing Y chromosomes first) to avoid BigQuery treat the chr column as integers
df_sorted = df.sort_values(by=['chr'], ascending = False)

################################## Save file in JSON format

# Save to JSON
df.to_json(OUTPUT_FILE, orient = "records", lines = True)