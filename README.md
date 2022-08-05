# ChemSampler
Sample small molecules, both from large chemical libraries as well as generative models

## Usage

```python
from chemsampler import ChemSampler
from chemsampler import example

smiles_list = example()
sampler = ChemSampler()
sampled_smiles = sampler.sample(smiles_list, num_samples=1000, sim_ub=0.95, sim_lb=0.6, distribution="ramp")
```