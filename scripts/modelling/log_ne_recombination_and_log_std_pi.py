import arviz as az
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm
import bambi as bmb
import xarray as xr
import random
import os
df = pd.read_csv('../../results/model/Nested_model/Reduced_Ne_Pi_recomb.csv')
df = df[df['freq_mean']>0.5]
df = df[df['chr']!='chrX']
grouped_df = df.groupby('full_species')
means = pd.DataFrame()
means['mean_pi_species'] = grouped_df['PI'].mean()
means['sd_pi_species'] = grouped_df['PI'].std()
means = means[['mean_pi_species', 'sd_pi_species']]
# Merge dataframes on the 'full_species' column
merged_df = pd.merge(df, means, on='full_species')
# Reset the index to remove the default index column
merged_df = merged_df.reset_index(drop=True)
# remove real zeros
merged_df = merged_df[merged_df['cm_per_mb']!=0]
merged_df['z_pi'] = (merged_df['PI']-merged_df['mean_pi_species'])/merged_df['sd_pi_species']
merged_df['log10_cm_per_mb'] = np.log10(merged_df['cm_per_mb'])
merged_df['log10_Ne'] = np.log10(merged_df['NE_MEDIAN'])
z_df = merged_df[['full_species','log10_Ne','z_pi', 'log10_cm_per_mb']]
unique_Species = z_df['full_species'].unique()
species_lookup = dict(zip(unique_Species, range(len(unique_Species))))
Ne =  (pd.DataFrame([z_df['full_species'], z_df['log10_Ne']]).transpose()).drop_duplicates()
pi = z_df['z_pi'].values
recombinationrate = z_df['log10_cm_per_mb'].values
species = z_df['full_species'].replace(species_lookup).values

sd_pi_logne_log_cm_model = pm.Model(coords = {"Species": unique_Species, 
                                        "obs_id": np.arange(len(recombinationrate))})
non_centered = True

with sd_pi_logne_log_cm_model:
# Data
    recomb = pm.ConstantData('recomb', recombinationrate, dims = 'obs_id')
    pi = pm.ConstantData('pi', pi, dims = 'obs_id')
    sp = pm.ConstantData('sp',species, dims = 'obs_id')
    Ne =  pm.ConstantData('Ne', Ne['log10_Ne'], dims = 'Species')

# Hyperpriors:
    g0 = pm.Normal("g0", mu=0, sigma=1)
    g1 = pm.Normal("g1", mu=0, sigma=1)
    h0 = pm.Normal("h0", mu=0, sigma=1)
    h1 = pm.Normal("h1", mu=0, sigma=1)

    mu_a = g0+g1*Ne
    mu_b = h0+h1*Ne
    sigma_a = pm.Exponential("sigma_a", 1)
    sigma_b = pm.Exponential("sigma_b", 1)
    if non_centered == True:
    # Varying intercepts:
        a_offset = pm.Normal("a_offset", 0, sigma=1, dims="Species")
        a = pm.Deterministic("a", mu_a + a_offset * sigma_a, dims="Species")
        # Varying slopes:
        b_offset = pm.Normal("b_offset", 0, sigma=1, dims="Species")
        b = pm.Deterministic("b", mu_b + b_offset * sigma_b, dims="Species")
# Expected value per species:
    y = a[sp] + b[sp] * recomb
# Model error
    sigma = pm.Exponential("sigma", 0.01)
    Pi = pm.Normal("Pi", mu=y, sigma=sigma, observed=pi, dims="obs_id")


trace_file = '../../results/model/Nested_model/log_ne_log_cm_sd_pi/log_ne_log_cm_sd_pi_model_24_04_2024.nc'
if os.path.exists(trace_file):
    with sd_pi_logne_log_cm_model:
         sd_pi_logne_log_cm_model_idata = az.from_netcdf(trace_file)
else:
    with sd_pi_logne_log_cm_model:
        sd_pi_logne_log_cm_model_idata = pm.sample(2000, target_accept=0.95, return_inferencedata=True,
                                 progressbar=True, cores=4, chains=4)
        sd_pi_logne_log_cm_model_idata.to_netcdf(trace_file)