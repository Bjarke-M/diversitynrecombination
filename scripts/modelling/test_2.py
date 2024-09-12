import pymc as pm
import pandas as pd
import numpy as np



df = pd.DataFrame({
    'obs_id': [1, 2, 3, 4, 5, 6, 7, 8],
    'group': [1, 1, 1, 2, 2, 3, 3, 3],
    'subject': [1, 2, 2, 2, 2, 1, 1, 2],
    'replicate': [1, 2, 1, 2, 1, 2, 1, 2]
})


group_idx, groups = pd.factorize(df["group"], sort=True)
replicate_idx, replicates = pd.factorize(df["replicate"], sort=True)

coords = {
    "group": groups,
    "obs_id": np.arange(df.shape[0]),
    "replicate":replicates,
}

#print(coords['group'])

model = pm.Model(coords=coords)

with model:
    
    mu_a = pm.Normal("μ_a", mu=0, sigma=10, dims="group")
    sigma_a = pm.Normal("σ_a", mu=1, sigma=10, dims="group")
    
    sigma_b = pm.Normal("σ_b",mu=1 ,sigma=10, dims="group")

    mu = pm.Normal("μ", mu=mu_a, sigma=sigma_a, dims=("replicate","group"))
    sigma = pm.Normal("σ",mu=1, sigma=sigma_b, dims=("replicate", "group"))

    obs = pm.Normal(
        "obs",
        mu=mu[replicate_idx, group_idx],
        sigma=sigma[replicate_idx, group_idx],
        observed=df["obs_id"],
        dims="obs_id",
    )
    idata = pm.sample(draws=2_000, tune=1_000, chains=4, return_inferencedata=True)
