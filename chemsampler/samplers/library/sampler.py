import tempfile
import os, sys
import csv
import shutil
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "tools", "fpsim2"))

from FPSim2.io import create_db_file
from FPSim2 import FPSim2Engine


class LibrarySampler(object):

    def __init__(self, n_workers=4):
        self.n_workers = n_workers
        self.library_dir = tempfile.mkdtemp()
        self.fp_filename = os.path.join(self.library_dir, "fp.h5")
        self.db_smiles_filename = self.get_db_filename()
        self.engine = None
        self.db_smiles = None

    def get_db_filename(self):
        return self.fp_filename[-3:]+".csv"

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

    def search_by_similarity_cutoff(self, smiles, cutoff=0.7):
        results = self.engine.similarity(smiles, cutoff, n_workers=self.n_workers)
        results = [(r[0], self.db_smiles[r[0]], r[1]) for r in results]
        return results

    def sample(self, smiles_list, n):
        smiles = set(smiles_list)
        sampled_smiles = set()
        for i in range(n):
            new_smiles = self._sample(list(smiles), n)
            sampled_smiles.update(new_smiles)
            smiles.update(new_smiles)
            if self.finished:
                break
        smiles = list(sampled_smiles)
        self.elapsed_time = 0
        self.finished = False
        return smiles

    def save(self, save_dir):
        if os.path.exists(save_dir):
            shutil.rmtree(save_dir)
        new_fp_filename = os.path.join(save_dir, "fp.h5")
        shutil.move(self.fp_filename, )
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