
import arviz as az
import pandas as pd
import numpy as np
import pymc as pm
import sys
import os

#### Load data
df = pd.read_csv('../../results/model/Nested_model/Reduced_Ne_Pi_recomb.csv')
pg_names = pd.read_csv('../../data/genus_phylogenetic_group_metadata.txt', sep="\t")
#### Correct for missingness in the data
df = df[df['freq_mean']>0.7]
df = df[df['chr']!='chrX']
df['corrected_pi'] = df['PI']*df['freq_mean']
df = pd.merge(df, pg_names, on='genus')
df = df.reset_index()

# add the phylogenetic group names to the dataframe
grouped_df = pd.merge(df, pg_names, on='genus')
grouped_df = grouped_df.reset_index()

# group the dataframe by species
grouped_df = df.groupby('full_species')



####
trace_file = sys.argv[2]
date = sys.argv[3]
####
non_centered = True

########################################################################################################################
############################################# ---- HIERARCHICAL MODEL ---- ###############################################
########################################################################################################################
# if sys.argv[1] == 'hierarchical_model':
#     if os.path.exists(trace_file):
#         print("Trace exists, here {trace_file}")
#     else:
#         ### Linear model of pi ~Recombination rate, with effect of Ne
#         ### Normalize the data
#         ### Step 1: calculate means and standard diviations the data
#         means = pd.DataFrame()
#         means[['mean_pi_species','mean_cm_per_mb_species']] = grouped_df[['corrected_pi', 'cm_per_mb']].mean()
#         means[['sd_pi_species','sd_cm_per_mb_species']] = grouped_df[['corrected_pi', 'cm_per_mb']].std()
#         means['primate_mean_ne'] = df['NE_MEDIAN'].mean()
#         means['primate_sd_ne'] = df['NE_MEDIAN'].std()
#         means = means[['mean_pi_species','mean_cm_per_mb_species', 'sd_pi_species','sd_cm_per_mb_species', 'primate_mean_ne', 'primate_sd_ne']]
#         ### Step 2: Normalize the data
#         # Merge dataframes on the 'full_species' column
#         merged_df = pd.merge(df, means, on='full_species')
#         #Reset the index to remove the default index column
#         merged_df = merged_df.reset_index(drop=True)
#         merged_df['z_pi'] = (merged_df['corrected_pi']-merged_df['mean_pi_species'])/merged_df['sd_pi_species']
#         merged_df['z_cm_per_mb'] = (merged_df['cm_per_mb']-merged_df['mean_cm_per_mb_species'])/merged_df['sd_cm_per_mb_species']
#         merged_df['z_ne'] = (merged_df['NE_MEDIAN']-merged_df['primate_mean_ne'])/merged_df['primate_sd_ne']
#         z_df = merged_df[['full_species','NE_MEDIAN','z_pi', 'z_cm_per_mb', 'z_ne']]
#         ### Step 3: make data in a pymc3 format
#         unique_Species = z_df['full_species'].unique()
#         species_lookup = dict(zip(unique_Species, range(len(unique_Species))))
#         Ne =  (pd.DataFrame([z_df['full_species'], z_df['z_ne']]).transpose()).drop_duplicates()
#         pi = z_df['z_pi'].values
#         recombinationrate = z_df['z_cm_per_mb'].values
#         species = z_df['full_species'].replace(species_lookup).values
#         ### step 4: formulate the model
#         hierarchical_model = pm.Model(coords = {"Species": unique_Species, 
        #                                 "obs_id": np.arange(len(recombinationrate))})
        # with hierarchical_model:
        # # Data
        #     recomb = pm.ConstantData('recomb', recombinationrate, dims = 'obs_id')
        #     pi = pm.ConstantData('pi', pi, dims = 'obs_id')
        #     sp = pm.ConstantData('sp',species, dims = 'obs_id')
        #     Ne =  pm.ConstantData('Ne', Ne['z_ne'], dims = 'Species')

        # # Hyperpriors:
        #     g0 = pm.Normal("g0", mu=0, sigma=1)
        #     g1 = pm.Normal("g1", mu=0, sigma=1)
        #     h0 = pm.Normal("h0", mu=0, sigma=1)
        #     h1 = pm.Normal("h1", mu=0, sigma=1)

        #     mu_a = g0+g1*Ne
        #     mu_b = h0+h1*Ne
        #     sigma_a = pm.Exponential("sigma_a", 1)
        #     sigma_b = pm.Exponential("sigma_b", 1)
        #     if non_centered == True:
        #     # Varying intercepts:
        #         a_offset = pm.Normal("a_offset", 0, sigma=1, dims="Species")
        #         a = pm.Deterministic("a", mu_a + a_offset * sigma_a, dims="Species")
        #         # Varying slopes:
        #         b_offset = pm.Normal("b_offset", 0, sigma=1, dims="Species")
        #         b = pm.Deterministic("b", mu_b + b_offset * sigma_b, dims="Species")
        #     else:
            # Varying intercepts:
        #         a = pm.Normal("a", mu=mu_a, sigma=sigma_a, dims="Species")
        #         # Varying slopes:
        #         b = pm.Normal("b", mu=mu_b, sigma=sigma_b, dims="Species")
        #      # Expected value per species:
        #     y = a[sp] + b[sp] * recomb
        # # Model error
        #     sigma = pm.Exponential("sigma", 0.01)
        #     Pi = pm.Normal("Pi", mu=y, sigma=sigma, observed=pi, dims="obs_id")
        # # Sample posterior
        # with hierarchical_model:
        #     hier_reg_idata = pm.sample(2000, tune=2000, target_accept=0.99, return_inferencedata=True,
        #                                progressbar=True, cores=8, chains=4)
        # ### Save the trace and posterior
        # az.to_netcdf(hier_reg_idata, trace_file)
        # g0 = hier_reg_idata.posterior.g0.to_dataframe()
        # g0.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/g0_{date}", sep="\t")
        # g1 = hier_reg_idata.posterior.g1.to_dataframe()
        # g1.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/g1_{date}", sep="\t")
        # h0 = hier_reg_idata.posterior.h0.to_dataframe()
        # h0.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/h0_{date}", sep="\t")
        # h1 = hier_reg_idata.posterior.h1.to_dataframe()
        # h1.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/h1_{date}", sep="\t")
        # sigma_a = hier_reg_idata.posterior.sigma_a.to_dataframe()
        # sigma_a.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_a_{date}", sep="\t")
        # sigma_b = hier_reg_idata.posterior.sigma_b.to_dataframe()
        # sigma_b.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_b_{date}", sep="\t")
        # a = hier_reg_idata.posterior.a.to_dataframe()
        # a.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/a_{date}", sep="\t")
        # b = hier_reg_idata.posterior.b.to_dataframe()
        # b.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/b_{date}", sep="\t")
        # sgm = hier_reg_idata.posterior.sigma.to_dataframe()
        # sgm.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_model_{date}", sep="\t")
        # if non_centered is True:
        #     b_offset = hier_reg_idata.posterior.b_offset.to_dataframe()
        #     b_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/b_offset_{date}", sep="\t")
        #     a_offset = hier_reg_idata.posterior.a_offset.to_dataframe()
        #     a_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/a_offset_{date}", sep="\t")


########################################################################################################################
##################### ---- LOG TRANSFORMED NE AND RECOMBINATION RATE HIERARCHICAL MODEL ---- ###########################
########################################################################################################################
# if sys.argv[1] == 'hierarchical_log_scaled_ne_cmpermb_model':
#     if os.path.exists(trace_file):
#         print("Trace exists, here {trace_file}")
#     else:
#     #TAKING THE LOG TO BOTH NE AND CM_PER_MB
#         means_log = pd.DataFrame()
#         # Step 2: Calculate Mean and Standard Deviation
#         means_log['mean_pi_species'] = grouped_df['corrected_pi'].mean()
#         means_log['mean_cm_per_mb_species'] = np.log10(grouped_df['cm_per_mb'].mean())
#         means_log['sd_pi_species'] = grouped_df['corrected_pi'].std()
#         means_log['sd_cm_per_mb_species'] = np.log10(grouped_df['cm_per_mb'].std())
#         means_log['primate_mean_ne'] = np.log10(df['NE_MEDIAN'].mean())
#         means_log['primate_sd_ne'] = np.log10(df['NE_MEDIAN'].std())
#         means_log = means_log[['mean_pi_species','mean_cm_per_mb_species', 'sd_pi_species','sd_cm_per_mb_species', 'primate_mean_ne', 'primate_sd_ne']]
#         # Merge dataframes on the 'full_species' column
#         merged_log_df = pd.merge(df, means_log, on='full_species')
#         # Reset the index to remove the default index column
#         merged_log_df = merged_log_df.reset_index(drop=True).dropna()
#         #remove real zeros
#         merged_log_df = merged_log_df[merged_log_df['cm_per_mb']!=0]
#         # Step 2: Normalize the data
#         merged_log_df['z_pi'] = (merged_log_df['corrected_pi']-merged_log_df['mean_pi_species'])/merged_log_df['sd_pi_species']
#         merged_log_df['z_cm_per_mb'] = (np.log10(merged_log_df['cm_per_mb'])-merged_log_df['mean_cm_per_mb_species'])/merged_log_df['sd_cm_per_mb_species']
#         merged_log_df['z_ne'] = (np.log10(merged_log_df['NE_MEDIAN'])-merged_log_df['primate_mean_ne'])/merged_log_df['primate_sd_ne']
#         z_log_df = merged_log_df[['full_species','NE_MEDIAN','z_pi', 'z_cm_per_mb', 'z_ne']]
#         ### Step 3: make data in a pymc3 format 
#         unique_Species_log = z_log_df['full_species'].unique()
#         species_lookup_log = dict(zip(unique_Species_log, range(len(unique_Species_log))))
#         Ne_log =  (pd.DataFrame([z_log_df['full_species'], z_log_df['z_ne']]).transpose()).drop_duplicates()
#         pi_log = z_log_df['z_pi'].values
#         recombinationrate_log = z_log_df['z_cm_per_mb'].values
#         species_log = z_log_df['full_species'].replace(species_lookup_log).values
#         ### step 4: formulate the model
#         hierarchical_log_scaled_ne_cmpermb_model = pm.Model(coords = {"Species": unique_Species_log, 
        #                                 "obs_id": np.arange(len(recombinationrate_log))})
        # with hierarchical_log_scaled_ne_cmpermb_model:
        #     # Data
        #     recomb = pm.ConstantData('recomb', recombinationrate_log, dims = 'obs_id')
        #     pi = pm.ConstantData('pi', pi_log, dims = 'obs_id')
        #     sp = pm.ConstantData('sp',species_log, dims = 'obs_id')
        #     Ne =  pm.ConstantData('Ne', Ne_log['z_ne'], dims = 'Species')

        #     #Hyperpriors:
        #     g0 = pm.Normal("g0", mu=0, sigma=1)
        #     g1 = pm.Normal("g1", mu=0, sigma=1)
        #     h0 = pm.Normal("h0", mu=0, sigma=1)
        #     h1 = pm.Normal("h1", mu=0, sigma=1)

        #     mu_a = g0+g1*Ne
        #     mu_b = h0+h1*Ne
        #     sigma_a = pm.Exponential("sigma_a", 1)
        #     sigma_b = pm.Exponential("sigma_b", 1)
        #     if non_centered == True:
        #     # Varying intercepts:
        #         a_offset = pm.Normal("a_offset", 0, sigma=1, dims="Species")
        #         a = pm.Deterministic("a", mu_a + a_offset * sigma_a, dims="Species")
        #         # Varying slopes:
        #         b_offset = pm.Normal("b_offset", 0, sigma=1, dims="Species")
        #         b = pm.Deterministic("b", mu_b + b_offset * sigma_b, dims="Species")
        #     # Expected value per species:
        #     y = a[sp] + b[sp] * recomb
            # Model error
        #     sigma = pm.Exponential("sigma", 0.01)
        #     Pi = pm.Normal("Pi", mu=y, sigma=sigma, observed=pi, dims="obs_id")
        # # Sample posterior
        # with hierarchical_log_scaled_ne_cmpermb_model:
        #     hierarchical_log_scaled_ne_cmpermb_model_idata = pm.sample(2000, tune=2000, target_accept=0.99, return_inferencedata=True,
        #                                                                progressbar=True, cores=25, chains=4)
        # ### Save the trace and posterior
        # az.to_netcdf(hierarchical_log_scaled_ne_cmpermb_model_idata, '../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/test_hierarchical_log_scaled_ne_cmpermb_model_18_12_2023.nc')
        # az.to_netcdf(hierarchical_log_scaled_ne_cmpermb_model_idata, trace_file)
        # g0 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.g0.to_dataframe()
        # g0.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/g0_{date}", sep="\t")
        # g1 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.g1.to_dataframe()
        # g1.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/g1_{date}", sep="\t")
        # h0 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.h0.to_dataframe()
        # h0.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/h0_{date}", sep="\t")
        # h1 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.h1.to_dataframe()
        # h1.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/h1_{date}", sep="\t")
        # sigma_a = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.sigma_a.to_dataframe()
        # sigma_a.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_a_{date}", sep="\t")
        # sigma_b = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.sigma_b.to_dataframe()
        # sigma_b.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_b_{date}", sep="\t")
        # a = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.a.to_dataframe()
        # a.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/a_{date}", sep="\t")
        # b = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.b.to_dataframe()
        # b.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/b_{date}", sep="\t")
        # sgm = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.sigma.to_dataframe()
        # sgm.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_model_{date}", sep="\t")
        # if non_centered is True:
            # b_offset = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.b_offset.to_dataframe()
            # b_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/b_offset_{date}", sep="\t")
            # a_offset = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.a_offset.to_dataframe()
            # a_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/a_offset_{date}", sep="\t")
        


########################################################################################################################
############################################# ---- BASIC MODEL ---- ####################################################
########################################################################################################################
# if sys.argv[1] == 'basic_model':
#     if os.path.exists(trace_file):
#         print("Trace exists, here {trace_file}")
#     # else:
    #     ### Linear model of pi ~Recombination rate, with no effect of Ne
    #     ### Normalize the data
    #     ### Step 1: calculate means and standard diviations the data
    #     means = pd.DataFrame()
    #     means[['mean_pi_species','mean_cm_per_mb_species']] = grouped_df[['corrected_pi', 'cm_per_mb']].mean()
    #     means[['sd_pi_species','sd_cm_per_mb_species']] = grouped_df[['corrected_pi', 'cm_per_mb']].std()
    #     means['primate_mean_ne'] = df['NE_MEDIAN'].mean()
    #     means['primate_sd_ne'] = df['NE_MEDIAN'].std()
    #     means = means[['mean_pi_species','mean_cm_per_mb_species', 'sd_pi_species','sd_cm_per_mb_species', 'primate_mean_ne', 'primate_sd_ne']]
    #     ### Step 2: Normalize the data
    #     # Merge dataframes on the 'full_species' column
    #     merged_df = pd.merge(df, means, on='full_species')
    #     #Reset the index to remove the default index column
    #     merged_df = merged_df.reset_index(drop=True)
    #     merged_df['z_pi'] = (merged_df['corrected_pi']-merged_df['mean_pi_species'])/merged_df['sd_pi_species']
    #     merged_df['z_cm_per_mb'] = (merged_df['cm_per_mb']-merged_df['mean_cm_per_mb_species'])/merged_df['sd_cm_per_mb_species']
    #     merged_df['z_ne'] = (merged_df['NE_MEDIAN']-merged_df['primate_mean_ne'])/merged_df['primate_sd_ne']
    #     z_df = merged_df[['full_species','NE_MEDIAN','z_pi', 'z_cm_per_mb', 'z_ne']]
    #     ### Step 3: make data in a pymc3 format
    #     unique_Species = z_df['full_species'].unique()
    #     species_lookup = dict(zip(unique_Species, range(len(unique_Species))))
    #     Ne =  (pd.DataFrame([z_df['full_species'], z_df['z_ne']]).transpose()).drop_duplicates()
    #     pi = z_df['z_pi'].values
    #     recombinationrate = z_df['z_cm_per_mb'].values
    #     species = z_df['full_species'].replace(species_lookup).values
    #     ### step 4: formulate the model
    #     basic_model = pm.Model(coords = {"Species": unique_Species, 
    #                                     "obs_id": np.arange(len(recombinationrate))})
    #     with basic_model:
    #     # Data
    #         recomb = pm.ConstantData('recomb', recombinationrate, dims = 'obs_id')
    #         pi = pm.ConstantData('pi', pi, dims = 'obs_id')
    #         sp = pm.ConstantData('sp',species, dims = 'obs_id')

        #     sigma_a = pm.Exponential("sigma_a", 1)
        #     sigma_b = pm.Exponential("sigma_b", 1)
        #     if non_centered == True:
        #     # Varying intercepts:
        #         a_offset = pm.Normal("a_offset", 0, sigma=1, dims="Species")
        #         a = pm.Deterministic("a", 0 + a_offset * sigma_a, dims="Species")
        #         # Varying slopes:
        #         b_offset = pm.Normal("b_offset", 0, sigma=1, dims="Species")
        #         b = pm.Deterministic("b", 0 + b_offset * sigma_b, dims="Species")
        #     else:
        #     # Varying intercepts:
        #         a = pm.Normal("a", mu=0, sigma=sigma_a, dims="Species")
        #         # Varying slopes:
        #         b = pm.Normal("b", mu=0, sigma=sigma_b, dims="Species")
        #      # Expected value per species:
        #     y = a[sp] + b[sp] * recomb
        # # Model error
        #     sigma = pm.Exponential("sigma", 0.01)
        #     Pi = pm.Normal("Pi", mu=y, sigma=sigma, observed=pi, dims="obs_id")
        # # Sample posterior
        # with basic_model:
        #     basic_model = pm.sample(2000, tune=2000, target_accept=0.99, return_inferencedata=True,
        #                                progressbar=True, cores=10, chains=4)
        # ### Save the trace and posterior
        # az.to_netcdf(basic_model, trace_file)
        # g0 = basic_model.posterior.g0.to_dataframe()
        # g0.to_csv("../../results/model/Nested_model/basic_model/g0_{date}", sep="\t")
        # g1 = basic_model.posterior.g1.to_dataframe()
        # g1.to_csv("../../results/model/Nested_model/basic_model/g1_{date}", sep="\t")
        # h0 = basic_model.posterior.h0.to_dataframe()
        # h0.to_csv("../../results/model/Nested_model/basic_model/h0_{date}", sep="\t")
        # h1 = basic_model.posterior.h1.to_dataframe()
        # h1.to_csv("../../results/model/Nested_model/basic_model/h1_{date}", sep="\t")
        # sigma_a = basic_model.posterior.sigma_a.to_dataframe()
        # sigma_a.to_csv("../../results/model/Nested_model/basic_model/sigma_a_{date}", sep="\t")
        # sigma_b = basic_model.posterior.sigma_b.to_dataframe()
        # sigma_b.to_csv("../../results/model/Nested_model/basic_model/sigma_b_{date}", sep="\t")
        # a = basic_model.posterior.a.to_dataframe()
        # a.to_csv("../../results/model/Nested_model/basic_model/a_{date}", sep="\t")
        # b = basic_model.posterior.b.to_dataframe()
        # b.to_csv("../../results/model/Nested_model/basic_model/b_{date}", sep="\t")
        # sgm = basic_model.posterior.sigma.to_dataframe()
        # sgm.to_csv("../../results/model/Nested_model/basic_model/sigma_model_{date}", sep="\t")
        # if non_centered is True:
        #     b_offset = basic_model.posterior.b_offset.to_dataframe()
        #     b_offset.to_csv("../../results/model/Nested_model/basic_model/b_offset_{date}", sep="\t")
        #     a_offset = basic_model.posterior.a_offset.to_dataframe()
        #     a_offset.to_csv("../../results/model/Nested_model/basic_model/a_offset_{date}", sep="\t")

########################################################################################################################
##################### ---- LOG TRANSFORMED RECOMBINATION RATE AND nonlog pi/Ne HIERARCHICAL MODEL ---- ###########################
# ########################################################################################################################
# if sys.argv[1] == 'hierarchical_log_scaled_cmpermb_model':
#     if os.path.exists(trace_file):
#         print("Trace exists, here {trace_file}")
#     else:
#     #TAKING THE LOG TO BOTH NE AND CM_PER_MB
#         means_log = pd.DataFrame()
#         # Step 2: Calculate Mean and Standard Deviation
#         means_log['mean_pi_species'] = grouped_df['corrected_pi'].mean()
#         means_log['mean_cm_per_mb_species'] = np.log10(grouped_df['cm_per_mb'].mean())
#         means_log['sd_pi_species'] = grouped_df['corrected_pi'].std()
#         means_log['sd_cm_per_mb_species'] = np.log10(grouped_df['cm_per_mb'].std())
#         means_log['primate_mean_ne'] = np.log10(df['NE_MEDIAN'].mean())
#         means_log['primate_sd_ne'] = np.log10(df['NE_MEDIAN'].std())
#         means_log = means_log[['mean_pi_species','mean_cm_per_mb_species', 'sd_pi_species','sd_cm_per_mb_species', 'primate_mean_ne', 'primate_sd_ne']]
#         # Merge dataframes on the 'full_species' column
#         merged_log_df = pd.merge(df, means_log, on='full_species')
#         # Reset the index to remove the default index column
#         merged_log_df = merged_log_df.reset_index(drop=True).dropna()
#         #remove real zeros
#         merged_log_df = merged_log_df[merged_log_df['cm_per_mb']!=0]
#         # Step 2: Normalize the data
#         merged_log_df['z_pi'] = (merged_log_df['corrected_pi']-merged_log_df['mean_pi_species'])/merged_log_df['sd_pi_species']
#         merged_log_df['z_cm_per_mb'] = (np.log10(merged_log_df['cm_per_mb'])-merged_log_df['mean_cm_per_mb_species'])/merged_log_df['sd_cm_per_mb_species']
#         merged_log_df['z_ne'] = (np.log10(merged_log_df['NE_MEDIAN'])-merged_log_df['primate_mean_ne'])/merged_log_df['primate_sd_ne']
#         z_log_df = merged_log_df[['full_species','NE_MEDIAN','z_pi', 'z_cm_per_mb', 'z_ne']]
#         ### Step 3: make data in a pymc3 format 
#         unique_Species_log = z_log_df['full_species'].unique()
#         species_lookup_log = dict(zip(unique_Species_log, range(len(unique_Species_log))))
#         Ne_log =  (pd.DataFrame([z_log_df['full_species'], z_log_df['z_ne']]).transpose()).drop_duplicates()
#         pi_log = z_log_df['z_pi'].values
#         recombinationrate_log = z_log_df['z_cm_per_mb'].values
#         species_log = z_log_df['full_species'].replace(species_lookup_log).values
#         ### step 4: formulate the model
#         hierarchical_log_scaled_cmpermb_model = pm.Model(coords = {"Species": unique_Species_log, 
#                                         "obs_id": np.arange(len(recombinationrate_log))})
#         with hierarchical_log_scaled_cmpermb_model:
#             # Data
#             recomb = pm.ConstantData('recomb', recombinationrate_log, dims = 'obs_id')
#             pi = pm.ConstantData('pi', pi_log, dims = 'obs_id')
#             sp = pm.ConstantData('sp',species_log, dims = 'obs_id')
#             Ne =  pm.ConstantData('Ne', Ne_log['z_ne'], dims = 'Species')

#             #Hyperpriors:
#             g0 = pm.Normal("g0", mu=0, sigma=1)
#             g1 = pm.Normal("g1", mu=0, sigma=1)
#             h0 = pm.Normal("h0", mu=0, sigma=1)
#             h1 = pm.Normal("h1", mu=0, sigma=1)

#             mu_a = g0+g1*Ne
#             mu_b = h0+h1*Ne
#             sigma_a = pm.Exponential("sigma_a", 1)
#             sigma_b = pm.Exponential("sigma_b", 1)
#             if non_centered == True:
#             # Varying intercepts:
#                 a_offset = pm.Normal("a_offset", 0, sigma=1, dims="Species")
#                 a = pm.Deterministic("a", mu_a + a_offset * sigma_a, dims="Species")
#                 # Varying slopes:
#                 b_offset = pm.Normal("b_offset", 0, sigma=1, dims="Species")
#                 b = pm.Deterministic("b", mu_b + b_offset * sigma_b, dims="Species")
#             # Expected value per species:
#             y = a[sp] + b[sp] * recomb
#             # Model error
#             sigma = pm.Exponential("sigma", 0.01)
#             Pi = pm.Normal("Pi", mu=y, sigma=sigma, observed=pi, dims="obs_id")
#         # Sample posterior
#         with hierarchical_log_scaled_cmpermb_model:
#             hierarchical_log_scaled_cmpermb_model_idata = pm.sample(2000, tune=2000, target_accept=0.99, return_inferencedata=True,
#                                                                        progressbar=True, cores=8, chains=4)
#         ### Save the trace and posterior
#         az.to_netcdf(hierarchical_log_scaled_cmpermb_model_idata, '../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/test_hierarchical_log_scaled_cmpermb_model_18_12_2023.nc')
#         az.to_netcdf(hierarchical_log_scaled_cmpermb_model_idata, trace_file)
#         g0 = hierarchical_log_scaled_cmpermb_model_idata.posterior.g0.to_dataframe()
#         g0.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/g0_{date}", sep="\t")
#         g1 = hierarchical_log_scaled_cmpermb_model_idata.posterior.g1.to_dataframe()
#         g1.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/g1_{date}", sep="\t")
#         h0 = hierarchical_log_scaled_cmpermb_model_idata.posterior.h0.to_dataframe()
#         h0.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/h0_{date}", sep="\t")
#         h1 = hierarchical_log_scaled_cmpermb_model_idata.posterior.h1.to_dataframe()
#         h1.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/h1_{date}", sep="\t")
#         sigma_a = hierarchical_log_scaled_cmpermb_model_idata.posterior.sigma_a.to_dataframe()
#         sigma_a.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/sigma_a_{date}", sep="\t")
#         sigma_b = hierarchical_log_scaled_cmpermb_model_idata.posterior.sigma_b.to_dataframe()
#         sigma_b.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/sigma_b_{date}", sep="\t")
#         a = hierarchical_log_scaled_cmpermb_model_idata.posterior.a.to_dataframe()
#         a.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/a_{date}", sep="\t")
#         b = hierarchical_log_scaled_cmpermb_model_idata.posterior.b.to_dataframe()
#         b.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/b_{date}", sep="\t")
#         sgm = hierarchical_log_scaled_cmpermb_model_idata.posterior.sigma.to_dataframe()
#         sgm.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/sigma_model_{date}", sep="\t")
#         if non_centered is True:
#             b_offset = hierarchical_log_scaled_cmpermb_model_idata.posterior.b_offset.to_dataframe()
#             b_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/b_offset_{date}", sep="\t")
#             a_offset = hierarchical_log_scaled_cmpermb_model_idata.posterior.a_offset.to_dataframe()
#             a_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_cmpermb_model/a_offset_{date}", sep="\t")
        

########################################################################################################################
##################### ---- PG_Group specific hyperpriors LOG TRANSFORMED NE AND RECOMBINATION RATE HIERARCHICAL MODEL ---- ###########################
########################################################################################################################
if sys.argv[1] == 'group_specific_priors_hierarchical_log_scaled_ne_cmpermb_model':
    if os.path.exists(trace_file):
        print("Trace exists, here {trace_file}")
    else:
    #TAKING THE LOG TO BOTH NE AND CM_PER_MB
        means_log = pd.DataFrame()
        # Step 2: Calculate Mean and Standard Deviation
        means_log['mean_pi_species'] = grouped_df['corrected_pi'].mean()
        #print(means_log)
        means_log['mean_cm_per_mb_species'] = np.log10(grouped_df['cm_per_mb'].mean())
        means_log['sd_pi_species'] = grouped_df['corrected_pi'].std()
        means_log['sd_cm_per_mb_species'] = np.log10(grouped_df['cm_per_mb'].std())
        means_log['primate_mean_ne'] = np.log10(df['NE_MEDIAN'].mean())
        means_log['primate_sd_ne'] = np.log10(df['NE_MEDIAN'].std())

        
        means_log = means_log[['mean_pi_species','mean_cm_per_mb_species', 'sd_pi_species','sd_cm_per_mb_species', 'primate_mean_ne', 'primate_sd_ne']]


        # Merge dataframes on the 'full_species' column
        merged_log_df = pd.merge(df, means_log, on='full_species')

        # Reset the index to remove the default index column
        merged_log_df = merged_log_df.reset_index(drop=True).dropna()
        #remove real zeros
        merged_log_df = merged_log_df[merged_log_df['cm_per_mb']!=0]
        # # Step 2: Normalize the data
        merged_log_df['z_pi'] = (merged_log_df['corrected_pi']-merged_log_df['mean_pi_species'])/merged_log_df['sd_pi_species']
        merged_log_df['z_cm_per_mb'] = (np.log10(merged_log_df['cm_per_mb'])-merged_log_df['mean_cm_per_mb_species'])/merged_log_df['sd_cm_per_mb_species']
        merged_log_df['z_ne'] = (np.log10(merged_log_df['NE_MEDIAN'])-merged_log_df['primate_mean_ne'])/merged_log_df['primate_sd_ne']

       
        z_log_df = merged_log_df[['full_species','NE_MEDIAN','z_pi', 'z_cm_per_mb', 'z_ne','pg_name']]

        ### Step 3: make data in a pymc3 format 
        unique_Species_log = z_log_df['full_species'].unique()
        species_lookup_log = dict(zip(unique_Species_log, range(len(unique_Species_log))))
        species_log = z_log_df['full_species'].replace(species_lookup_log).values

        Ne_log =  (pd.DataFrame([z_log_df['full_species'], z_log_df['z_ne']]).transpose()).drop_duplicates()

        pi_log = z_log_df['z_pi'].values

        recombinationrate_log = z_log_df['z_cm_per_mb'].values


        unique_pg_log = z_log_df['pg_name'].unique()
        pg_lookup_log = dict(zip(unique_pg_log, range(len(unique_pg_log))))
        pg_log = z_log_df['pg_name'].replace(pg_lookup_log).values


        print('STARTING THE MODEL')
  
        ### step 4: formulate the model
        hierarchical_log_scaled_ne_cmpermb_model = pm.Model(coords = 
                                                            {"Species": np.unique(np.array(species_log)), 
                                                             "obs_id": np.arange(len(recombinationrate_log)),
                                                             "PhyloGroup": np.unique(np.array(pg_log))
                                                            })
        
        with hierarchical_log_scaled_ne_cmpermb_model:
            # Data
            recomb = pm.ConstantData('recomb', recombinationrate_log, dims = 'obs_id')
            pi = pm.ConstantData('pi', pi_log, dims = 'obs_id')
            sp = pm.ConstantData('sp',species_log, dims = 'obs_id')
            Ne =  pm.ConstantData('Ne', z_log_df['z_ne'], dims = 'obs_id') #can be changed back to Ne_log['z_ne'] dims = Species
            pg = pm.ConstantData('pg', pg_log, dims = 'obs_id')

            #Hyperpriors:
            g0 = pm.Normal("g0", mu=0, sigma=1, dims="PhyloGroup")
            g1 = pm.Normal("g1", mu=0, sigma=1, dims = "PhyloGroup")
            h0 = pm.Normal("h0", mu=0, sigma=1, dims="PhyloGroup")  # Add dimension
            h1 = pm.Normal("h1", mu=0, sigma=1, dims="PhyloGroup")

            mu_a = g0[pg]+g1[pg]*Ne
            mu_b = h0[pg]+h1[pg]*Ne

            sigma_a = pm.Exponential("sigma_a", 1)
            sigma_b = pm.Exponential("sigma_b", 1)


            # Varying intercepts:
            a_offset = pm.Normal("a_offset", 0, sigma=1, dims="Species")
            a = pm.Deterministic("a",mu_a + a_offset[sp] * sigma_a, dims="obs_id")
                # Varying slopes:
            b_offset = pm.Normal("b_offset", 0, sigma=1, dims="Species")
            b = pm.Deterministic("b", mu_b + b_offset[sp] * sigma_b, dims="obs_id")


            # Expected value per species:
            y = a[sp] + b[sp] * recomb

            # Model error
            sigma = pm.Exponential("sigma", 0.01)
            # Likehood
            Pi = pm.Normal("Pi", mu=y, sigma=sigma, observed=pi, dims="obs_id")

        # Sample posterior
        # with hierarchical_log_scaled_ne_cmpermb_model:
            hierarchical_log_scaled_ne_cmpermb_model_idata = pm.sample(500, tune=500, target_accept=0.99, return_inferencedata=True,
                                                                       progressbar=True, cores=20, chains=4)



        # # ### Save the trace and posterior
        # # # az.to_netcdf(hierarchical_log_scaled_ne_cmpermb_model_idata, '../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/test_hierarchical_log_scaled_ne_cmpermb_model_18_12_2023.nc')
        az.to_netcdf(hierarchical_log_scaled_ne_cmpermb_model_idata, "../../results/model/Nested_model/group_specific_priors_hierarchical_log_scaled_ne_cmpermb_model/group_specific_priors_hierarchical_log_scaled_ne_cmpermb_model_190324.nc")


        # # g0 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.g0.to_dataframe()
        # # g0.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/g0_{date}", sep="\t")
        # # g1 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.g1.to_dataframe()
        # # g1.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/g1_{date}", sep="\t")
        # h0 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.h0.to_dataframe()
        # h0.to_csv("../../results/model/Nested_model/test/h0_{date}", sep="\t")
        # h1 = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.h1.to_dataframe()
        # h1.to_csv("../../results/model/Nested_model/test/h1_{date}", sep="\t")
        # # sigma_a = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.sigma_a.to_dataframe()
        # # sigma_a.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_a_{date}", sep="\t")
        # # sigma_b = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.sigma_b.to_dataframe()
        # # sigma_b.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_b_{date}", sep="\t")
        # a = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.a.to_dataframe()
        # a.to_csv("../../results/model/Nested_model/test/a_{date}", sep="\t")
        # b = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.b.to_dataframe()
        # b.to_csv("../../results/model/Nested_model/test/b_{date}", sep="\t")
        # # sgm = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.sigma.to_dataframe()
        # sgm.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/sigma_model_{date}", sep="\t")
        # if non_centered is True:
            # b_offset = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.b_offset.to_dataframe()
            # b_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/b_offset_{date}", sep="\t")
            # a_offset = hierarchical_log_scaled_ne_cmpermb_model_idata.posterior.a_offset.to_dataframe()
            # a_offset.to_csv("../../results/model/Nested_model/hierarchical_log_scaled_ne_cmpermb_model/a_offset_{date}", sep="\t")