import random
from timeit import default_timer as timer

from .samplers.chembl.sampler import ChemblSampler
from .samplers.pubchem.sampler import PubChemSampler
from .samplers.smallworld.sampler import SmallWorldSampler
from .samplers.stoned.sampler import StonedSampler
from .samplers.mollib.sampler import MollibSampler


SAMPLERS_LIST = [ChemblSampler, PubChemSampler, SmallWorldSampler, StonedSampler, MollibSampler]


class ChemSampler(object):
    def __init__(self, inflation=10, max_greedy_iterations=10, samplers_list=None):
        if samplers_list is None:
            self.samplers_list = SAMPLERS_LIST
        else:
            self.samplers_list = samplers_list
        self.inflation = inflation
        self.max_greedy_iterations = max_greedy_iterations

    def _one_sampler(self, smiles_list):
        random.shuffle(smiles_list)
        Sampler = random.sample(self.samplers_list, 1)[0]
        print(Sampler)
        if Sampler == ChemblSampler:
            print("ChemblSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=smiles_list,
                n=self.num_samples,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == PubChemSampler:
            print("PubChemSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=smiles_list,
                n=self.num_samples,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == SmallWorldSampler:
            print("SmallWorldSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=smiles_list,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == StonedSampler:
            print("StonedSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=smiles_list,
                n=self.num_samples,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == MollibSampler:
            print("MollibSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=smiles_list,
                n=max(self.num_samples,100) # TODO check
            )

    def _greedy_sample(self, smiles_list, num_samples, time_budget_sec):
        t0 = timer()
        sampled_smiles = set()
        for _ in range(self.max_greedy_iterations):
            t1 = timer()
            if len(sampled_smiles) > (num_samples*self.inflation):
                return sampled_smiles
            if (t1-t0) > time_budget_sec:
                return sampled_smiles
            sampled_smiles.update(self._one_sampler(smiles_list))
            print(len(sampled_smiles))
        return sampled_smiles

    def sample(self, smiles_list, num_samples, time_budget_sec):
        self.time_budget_sec = int(time_budget_sec) + 1
        self.one_sampler_time_budget_sec = (
            int(self.time_budget_sec / len(self.samplers_list)) + 1
        )
        self.num_samples = num_samples
        sampled_smiles = self._greedy_sample(smiles_list, num_samples, time_budget_sec)
        sampled_smiles = list(sampled_smiles)
        random.shuffle(sampled_smiles)
        return sampled_smiles