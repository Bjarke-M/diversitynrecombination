import arviz as az
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm
import bambi as bmb
import xarray as xr
import random

#### Load data
df = pd.read_csv('../results/model/Nested_model/Reduced_Ne_Pi_recomb.csv')
#### Correct for missingness in the data
df = df[df['freq_mean']>0.7]
df = df[df['chr']!='chrX']
df['corrected_pi'] = df['PI']*df['freq_mean']
#####


########################### Basic Model ###########################
########################### Basic Model ###########################
########################### Basic Model ###########################
if sys.argv[1] == 'basic':
    ### Linear model of pi ~Recombination rate, with no effect of Ne
    ### Normalize the data
    ### Step 1: calculate means and standard diviations the data
    grouped_df = df.groupby('full_species')
    means = pd.DataFrame()
    means[['mean_pi_species','mean_cm_per_mb_species']] = grouped_df[['corrected_pi', 'cm_per_mb']].mean()
    means[['sd_pi_species','sd_cm_per_mb_species']] = grouped_df[['corrected_pi', 'cm_per_mb']].std()
    means['primate_mean_ne'] = df['NE_MEDIAN'].mean()
    means['primate_sd_ne'] = df['NE_MEDIAN'].std()
    means = means[['mean_pi_species','mean_cm_per_mb_species', 'sd_pi_species','sd_cm_per_mb_species', 'primate_mean_ne', 'primate_sd_ne']]
    ### Step 2: Normalize the data
