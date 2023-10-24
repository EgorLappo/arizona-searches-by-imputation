import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import os

d = pd.concat([pd.read_csv(f'../computed_matches/{f}')
               for f in os.listdir('../computed_matches')])\
    .sort_values(['rep'])\
    .reset_index(drop=True)\
    .drop(columns=['s1', 's2'])

# medians
print('median for true allele matches:')
print(d['true_allele_matches'].median())
print('median for beagle-called allele matches:')
print(d['called_allele_matches'].median())
print('median for beagle-ap matches:')
print(d['mean_allele_matches'].median())

print('median for true full matches:')
print(d['true_full_matches'].median())
print('median for beagle-called full matches:')
print(d['called_full_matches'].median())
print('median for beagle-ap full matches:')
print(d['mean_full_matches'].median())

print('median for true partial matches:')
print(d['true_partial_matches'].median())
print('median for beagle-called partial matches:')
print(d['called_partial_matches'].median())
print('median for beagle-ap partial matches:')
print(d['mean_partial_matches'].median())

# max numbers
print('max for true allele matches:')
ta_max = d['true_allele_matches'].max()
print(ta_max)
print('rows with max for true allele matches:')
print(d.loc[d['true_allele_matches'] ==
      ta_max, ['rep', 'true_allele_matches']])
print(
    f'total number of these values: {len(d.loc[d["true_allele_matches"] == ta_max])}')

print('max for beagle-called allele matches:')
ca_max = d['called_allele_matches'].max()
print(ca_max)
print('rows with max for beagle-called allele matches:')
print(d.loc[d['called_allele_matches'] ==
      ca_max, ['rep', 'called_allele_matches']])
print(
    f'total number of these values: {len(d.loc[d["called_allele_matches"] == ca_max])}')

print('max for beagle-ap allele matches:')
ma_max = d['mean_allele_matches'].max()
print(ma_max)
print('rows with max for beagle-ap allele matches:')
print(d.loc[d['mean_allele_matches'] ==
            ma_max, ['rep', 'mean_allele_matches']])
print(
    f'total number of these values: {len(d.loc[d["mean_allele_matches"] == ma_max])}')

print('max for true full matches:')
tf_max = d['true_full_matches'].max()
print(tf_max)
print('rows with max for true full matches:')
print(d.loc[d['true_full_matches'] ==
            tf_max, ['rep', 'true_full_matches']])
print(
    f'total number of these values: {len(d.loc[d["true_full_matches"] == tf_max])}')

print('max for beagle-called full matches:')
cf_max = d['called_full_matches'].max()
print(cf_max)
print('rows with max for beagle-called full matches:')
print(d.loc[d['called_full_matches'] ==
            cf_max, ['rep', 'called_full_matches']])
print(
    f'total number of these values: {len(d.loc[d["called_full_matches"] == cf_max])}')

print('max for beagle-ap full matches:')
mf_max = d['mean_full_matches'].max()
print(mf_max)
print('rows with max for beagle-ap full matches:')
print(d.loc[d['mean_full_matches'] ==
            mf_max, ['rep', 'mean_full_matches']])
print(
    f'total number of these values: {len(d.loc[d["mean_full_matches"] == mf_max])}')

print('max for true partial matches:')
tp_max = d['true_partial_matches'].max()
print(tp_max)
print('rows with max for true partial matches:')
print(d.loc[d['true_partial_matches'] ==
            tp_max, ['rep', 'true_partial_matches']])
print(
    f'total number of these values: {len(d.loc[d["true_partial_matches"] == tp_max])}')

print('max for beagle-called partial matches:')
cp_max = d['called_partial_matches'].max()
print(cp_max)
print('rows with max for beagle-called partial matches:')
print(d.loc[d['called_partial_matches'] ==
            cp_max, ['rep', 'called_partial_matches']])
print(
    f'total number of these values: {len(d.loc[d["called_partial_matches"] == cp_max])}')

print('max for beagle-ap partial matches:')
mp_max = d['mean_partial_matches'].max()
print(mp_max)
print('rows with max for beagle-ap partial matches:')
print(d.loc[d['mean_partial_matches'] ==
            mp_max, ['rep', 'mean_partial_matches']])
print(
    f'total number of these values: {len(d.loc[d["mean_partial_matches"] == mp_max])}')
