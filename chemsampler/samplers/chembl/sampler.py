import requests
import random
from tqdm import tqdm
from timeit import default_timer as timer
from ratelimit import sleep_and_retry, limits
from standardiser import standardise
from rdkit import Chem


@sleep_and_retry
@limits(calls=2, period=30)
def run_chembl_sampler(origin_smiles: str, similarity: float = 0.4):
    """Function adapted from Andrew White's Exmol"""
    similarity = int(similarity * 100)
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{origin_smiles}/{similarity}?format=json"
    try:
        reply = requests.get(url, timeout=10)
    except requests.exceptions.Timeout:
        print("ChEMBL seems to be down right now")
        return []
    try:
        data = reply.json()
    except:
        return []
    smiles = [d["molecule_structures"]["canonical_smiles"] for d in data["molecules"]]
    smiles = list(set(smiles))

    return smiles


class ChemblSampler(object):
    def __init__(self):
        self.inflation = 2

    def _sample(self, smiles):
        smiles = run_chembl_sampler(origin_smiles=smiles)
        std_smiles = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                mol = standardise.run(mol)
            except:
                continue
            std_smiles += [Chem.MolToSmiles(mol)]
        return std_smiles

    def sample(self, smiles_list, n, time_budget_sec=600):
        sampled_smiles = []
        t0 = timer()
        for smi in tqdm(smiles_list):
            sampled_smiles += self._sample(smi)
            t1 = timer()
            if (t1 - t0) > time_budget_sec:
                break
            if len(set(sampled_smiles)) > n:
                break
        sampled_smiles = list(set(sampled_smiles))
        random.shuffle(sampled_smiles)
        return sampled_smiles
