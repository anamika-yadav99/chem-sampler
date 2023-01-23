from chemsampler import ChemSampler
from chemsampler import example

smiles_list = example()
sampler = ChemSampler()
sampled_smiles = sampler.sample(smiles_list, num_samples=100, sim_ub=0.95, sim_lb=0.6, distribution="ramp")
