import pandas as pd
import pymc as pm
import numpy as np

# Step 1: Prepare the data
df = pd.DataFrame({
    'species': ['A', 'A', 'A', 'A', 'B', 'B','B', 'B', 'C', 'C', 'C', 'C'],     
    'recomb': [1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3], 
    'pi': [0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.7, 0.9, 0.9, 0.3, 0.4, 0.6], 
    'ne': [100,100,100,100,200,200,200,200,300,300,300,300], 
    'pg': ['pg1', 'pg1', 'pg1', 'pg1', 'pg2', 'pg2', 'pg2', 'pg2', 'pg1', 'pg1', 'pg1', 'pg1']
})
print(df[df['species_idx']==0]['recomb'])

# Step 2: Convert categorical variables to numerical indices
group_idx, groups = pd.factorize(df["pg"], sort=True)
species_idx, species = pd.factorize(df["species"], sort=True)



print(species_idx)
# Step 3: Define the model
coords = {
    "group": groups,
    "obs_id": np.arange(df.shape[0]),
    "replicate": species,
}

model = pm.Model(coords=coords)

with model:
    # Hyperpriors
    g0 = pm.Normal("g0", mu=0, sigma=1, dims="group")
    g1 = pm.Normal("g1", mu=0, sigma=1, dims="group")
    h0 = pm.Normal("h0", mu=0, sigma=1, dims="group")
    h1 = pm.Normal("h1", mu=0, sigma=1, dims="group")

    # Linear models
    mu_a = g0[group_idx] + g1[group_idx] * df['ne']
    mu_b = h0[group_idx] + h1[group_idx] * df['ne']

    a = pm.Deterministic("a", mu_a, dims="replicate")
    b = pm.Deterministic("b", mu_b, dims="replicate")

    # Expected value
    y = a[species_idx] + b[species_idx] * df['recomb']

    # Model error
    sigma = pm.Exponential("sigma", 0.01)

    # Likelihood
    Pi = pm.Normal("Pi", mu=y, sigma=sigma, observed=df['pi'], dims="obs_id")

    idata = pm.sample(draws=2000, tune=1000, chains=1, return_inferencedata=True) 






#     import pandas as pd
# import numpy as np

# # True values
# true_g0 = np.array([0.1, 0.2])  # for 'pg1' and 'pg2'
# true_g1 = np.array([0.3, 0.4])  # for 'pg1' and 'pg2'
# true_h0 = np.array([0.5, 0.6])  # for 'pg1' and 'pg2'
# true_h1 = np.array([0.7, 0.8])  # for 'pg1' and 'pg2'

# # Create a dataframe
# df = pd.DataFrame({
#     'species': ['A', 'A', 'A', 'A', 'B', 'B','B', 'B', 'C', 'C', 'C', 'C'],     
#     'recomb': [1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3], 
#     'ne': [100,100,100,100,200,200,200,200,300,300,300,300], 
#     'pg': ['pg1', 'pg1', 'pg1', 'pg1', 'pg2', 'pg2', 'pg2', 'pg2', 'pg1', 'pg1', 'pg1', 'pg1']
# })

# # Convert categorical variables to numerical indices
# group_idx, groups = pd.factorize(df["pg"], sort=True)
# species_idx, species = pd.factorize(df["species"], sort=True)

# # Calculate 'a' and 'b' using the true values and 'ne'
# a = true_g0[group_idx] + true_g1[group_idx] * df['ne']
# b = true_h0[group_idx] + true_h1[group_idx] * df['ne']

# # Calculate 'pi' using 'a', 'b' and 'recomb'
# df['pi'] = a[species_idx] + b[species_idx] * df['recomb'] + np.random.normal(0, 0.01, df.shape[0])  # add some noise


