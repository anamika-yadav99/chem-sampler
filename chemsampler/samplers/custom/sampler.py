import tempfile
import os, sys
import csv
import shutil
import random
import time
import numpy as np
from tqdm import tqdm

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "tools", "fpsim2"))

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem

from FPSim2.io import create_db_file
from FPSim2 import FPSim2Engine


class CustomLibrarySampler(object):
    def __init__(
        self, min_similarity=0.6, max_similarity=0.9, time_budget_sec=60, n_workers=4
    ):
        self.min_similarity = min_similarity
        self.max_similarity = max_similarity
        self.n_workers = n_workers
        self.library_dir = tempfile.mkdtemp()
        self.fp_filename = os.path.join(self.library_dir, "fp.h5")
        self.db_smiles_filename = self.get_db_filename()
        self.engine = None
        self.db_smiles = None
        self.time_budget_sec = time_budget_sec
        self.elapsed_time = 0
        self.finished = False

    def get_db_filename(self):
        return self.fp_filename[-3:] + ".csv"

    def read_db_smiles(self):
        smiles_list = []
        with open(self.db_smiles_filename, "r") as f:
            reader = csv.reader(f)
            for r in reader:
                smiles_list += [r[0]]
        return smiles_list

    def write_db_smiles(self, smiles_list):
        with open(self.db_smiles_filename, "w") as f:
            writer = csv.writer(f)
            for smi in smiles_list:
                writer.writerow([smi])

    def _load_engine(self):
        self.engine = FPSim2Engine(self.fp_filename)
        self.db_smiles = self.read_db_smiles()

    def fit(self, smiles_list):
        self.write_db_smiles(smiles_list)
        smiles_list = [[smi, i] for i, smi in enumerate(smiles_list)]
        create_db_file(
            smiles_list, self.fp_filename, "Morgan", {"radius": 2, "nBits": 2048}
        )
        self._load_engine()

    def search_by_similarity_cutoff(self, smiles, cutoff):
        results = self.engine.similarity(smiles, cutoff, n_workers=self.n_workers)
        results = [(r[0], self.db_smiles[r[0]], r[1]) for r in results]
        return results

    def _select_by_similarity(self, smiles):
        sel_smiles = []
        for smi in tqdm(smiles):
            mol = Chem.MolFromSmiles(smi)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            sims = DataStructs.BulkTanimotoSimilarity(fp, self.seed_fps)
            sim = np.max(sims)
            if sim < self.min_similarity or sim > self.max_similarity:
                continue
            sel_smiles += [smi]
        return sel_smiles

    def _sample(self, smiles_list, cutoff):
        random.shuffle(smiles_list)
        sampled_smiles = []
        sampled_sim = []
        for smi in tqdm(smiles_list):
            t0 = time.time()
            results = self.search_by_similarity_cutoff(self, smi, cutoff)
            sampled = [r[1] for r in results]
            sampled_smiles += sampled[0]
            sampled_sim += sampled[1]
            t1 = time.time()
            dt = t1 - t0
            self.elapsed_time += dt
            if self.elapsed_time > self.time_budget_sec:
                self.finished = True
                break
        smiles = []
        for smi, sim in zip(sampled_smiles, sampled_sim):
            if sim < self.min_similarity or sim > self.max_similarity:
                continue
            smiles += [smi]
        smiles = self._select_by_similarity(smiles)
        return set(smiles)

    def sample(self, smiles_list, n):
        self.seed_smiles = list(smiles_list)
        self.seed_fps = [
            AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 2)
            for smi in self.seed_smiles
        ]
        smiles = set(smiles_list)
        sampled_smiles = set()
        for i in range(n):
            new_smiles = self._sample(list(smiles))
            sampled_smiles.update(new_smiles)
            smiles.update(new_smiles)
            if self.finished:
                break
        smiles = list(sampled_smiles)
        smiles = self._select_by_similarity(smiles)
        self.elapsed_time = 0
        self.finished = False
        return smiles

    def save(self, save_dir):
        if os.path.exists(save_dir):
            shutil.rmtree(save_dir)
        new_fp_filename = os.path.join(save_dir, "fp.h5")
        shutil.move(self.fp_filename)
        self.fp_filename = new_fp_filename
        new_db_smiles_filename = self.get_db_filename()
        shutil.move(self.db_smiles_filename, new_db_smiles_filename)
        self.db_smiles_filename = new_db_smiles_filename
        self.library_dir = os.path.abspath(save_dir)

    def load(self, library_dir):
        self.library_dir = os.path.abspath(library_dir)
        self.fp_filename = os.path.join(self.library_dir, "fp.h5")
        self.db_smiles_filename = self.get_db_filename()
        self._load_engine()
