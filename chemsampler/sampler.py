import random
from timeit import default_timer as timer
import numpy as np
from tqdm import tqdm

from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import DataStructs

from .samplers.chembl.sampler import ChemblSampler
from .samplers.pubchem.sampler import PubChemSampler
from .samplers.smallworld.sampler import SmallWorldSampler
from .samplers.stoned.sampler import StonedSampler
from .samplers.mollib.sampler import MollibSampler
from .samplers.fasmifra.sampler import FasmifraSampler


SAMPLERS_LIST = [
    ChemblSampler,
    PubChemSampler,
    SmallWorldSampler,
    StonedSampler,
    MollibSampler,
]


class ChemSampler(object):
    def __init__(self, inflation=10, max_greedy_iterations=10, samplers_list=None):
        if samplers_list is None:
            self.samplers_list = SAMPLERS_LIST
        else:
            self.samplers_list = samplers_list
        self.inflation = inflation
        self.max_greedy_iterations = max_greedy_iterations
        self.small_list_size = 10
        self.sampled_smiles = None

    def _one_sampler(self, smiles_list):
        random.shuffle(smiles_list)
        small_smiles_list = smiles_list[: self.small_list_size]
        Sampler = random.sample(self.samplers_list, 1)[0]
        print(Sampler)
        if Sampler == ChemblSampler:
            print("ChemblSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=small_smiles_list,
                n=self.num_samples,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == PubChemSampler:
            print("PubChemSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=small_smiles_list,
                n=self.num_samples,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == SmallWorldSampler:
            print("SmallWorldSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=small_smiles_list,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == StonedSampler:
            print("StonedSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=small_smiles_list,
                n=self.num_samples,
                time_budget_sec=self.one_sampler_time_budget_sec,
            )
        if Sampler == MollibSampler:
            print("MollibSampler")
            sampler = Sampler()
            return sampler.sample(
                smiles_list=small_smiles_list,
                n=max(self.num_samples, 100),  # TODO check
            )
        if Sampler == FasmifraSampler:
            print("FasmifraSampler")
            sampler = Sampler()
            return sampler.sample(smiles_list=small_smiles_list, n=self.num_samples)

    def _greedy_sample(self, smiles_list, num_samples, time_budget_sec):
        t0 = timer()
        sampled_smiles = set()
        for _ in range(self.max_greedy_iterations):
            t1 = timer()
            if len(sampled_smiles) > (num_samples * self.inflation):
                return sampled_smiles
            if (t1 - t0) > time_budget_sec:
                return sampled_smiles
            sampled_smiles.update(self._one_sampler(smiles_list))
            print(len(sampled_smiles))
        return sampled_smiles

    def more(self, smiles_list, num_samples=100, time_budget_sec=60):
        self.time_budget_sec = int(time_budget_sec) + 1
        self.one_sampler_time_budget_sec = (
            int(self.time_budget_sec / len(self.samplers_list)) + 1
        )
        self.num_samples = num_samples
        sampled_smiles = self._greedy_sample(smiles_list, num_samples, time_budget_sec)
        sampled_smiles = list(sampled_smiles)
        if self.sampled_smiles is None:
            self.sampled_smiles = set(sampled_smiles)
        else:
            self.sampled_smiles.update(sampled_smiles)

    def subsample(
        self,
        smiles_list,
        sampled_smiles,
        sim_ub,
        sim_lb,
        distribution="normal",
        num_per_sample=10,
    ):
        steps = 100
        offset = 0.1
        num_atoms_proportion_difference = 0.5
        assert sim_ub > sim_lb
        if distribution == "normal":
            mean = (sim_ub + sim_lb) / 2
            std = (sim_ub - sim_lb) / 5
            values = np.random.normal(loc=mean, scale=std, size=steps)
            idxs = np.argsort(values)
            values = values[idxs]
            weights = None
        elif distribution == "ramp":
            values = np.linspace(sim_lb, sim_ub, steps)
            values = values
            weights = (values - np.min(values)) / (np.max(values) - np.min(values))
            weights = weights / np.sum(weights)
        elif distribution == "uniform":
            values = np.linspace(sim_lb, sim_ub, steps)
            values = values
            weights = None
        else:
            raise Exception("distribution must be 'normal', 'ramp' 'uniform'")
        input_fps = [
            AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 3)
            for smi in smiles_list
        ]
        sampled_fps = [
            AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 3)
            for smi in sampled_smiles
        ]
        R = []
        for input_idx, input_fp in tqdm(enumerate(input_fps)):
            sims = DataStructs.BulkTanimotoSimilarity(input_fp, sampled_fps)
            input_num_atoms = Chem.MolFromSmiles(smiles_list[input_idx]).GetNumAtoms()
            selected_idxs = set()
            for _ in range(num_per_sample):
                ref_sim = np.random.choice(values, size=None, p=weights)
                distance = np.abs(sims - ref_sim)
                idx = np.argmin(distance)
                if distance[idx] < offset:
                    sampled_num_atoms = Chem.MolFromSmiles(sampled_smiles[idx]).GetNumAtoms()
                    if (
                        (abs(sampled_num_atoms - input_num_atoms) / input_num_atoms)
                        < num_atoms_proportion_difference
                    ):
                        selected_idxs.update([idx])
            if len(selected_idxs) != 0:
                selected_idxs = sorted(selected_idxs)
                selected_sims = DataStructs.BulkTanimotoSimilarity(
                    input_fp, [sampled_fps[i] for i in selected_idxs]
                )
                selected_smiles = [sampled_smiles[i] for i in selected_idxs]
                ordered_idxs = np.argsort(selected_sims)[::-1]
                selected_smiles = [selected_smiles[i] for i in ordered_idxs]
            else:
                selected_smiles = []
            R += [selected_smiles]
        return R

    def sample(
        self,
        smiles_list,
        num_samples=100,
        sim_ub=0.95,
        sim_lb=0.4,
        distribution="ramp",
        time_budget_sec=60,
        flatten=False,
    ):
        num_per_sample = max(3, int(num_samples / len(smiles_list)))
        if self.sampled_smiles is None:
            self.more(
                smiles_list=smiles_list,
                num_samples=num_samples,
                time_budget_sec=time_budget_sec,
            )
        self.time_budget_sec = int(time_budget_sec) + 1
        self.one_sampler_time_budget_sec = (
            int(self.time_budget_sec / len(self.samplers_list)) + 1
        )
        R = self.subsample(
            smiles_list=smiles_list,
            sampled_smiles=list(self.sampled_smiles),
            sim_ub=sim_ub,
            sim_lb=sim_lb,
            distribution=distribution,
            num_per_sample=num_per_sample,
        )
        data = R
        if flatten:
            flat_data = list(set([x for d in data for x in d]))
            random.shuffle(flat_data)
            return flat_data
        else:
            return data
