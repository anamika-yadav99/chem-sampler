import os
import shutil
import configparser
import tempfile
import pandas as pd
import numpy as np

from ...utils.conda import SimpleConda

root = os.path.dirname(os.path.abspath(__file__))

MOLER_CONDA_ENVIRONMENT = "moler"

class moler_sampler:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.CHECKPOINT_DIR = os.path.join(root, "molecule_generation", "MODEL_DIR" )

    def _sample(self, n):
        cmd = "cd {0}; bash run_sample.sh {1} {2}; cd {3}".format(
            self.exec_folder, self.CHECKPOINT_DIR, n, self.cwd)
        SimpleConda().run_commandlines(MOLER_CONDA_ENVIRONMENT, cmd)
        
    def _read_molecules(self):
        output_folder = os.path.join(
            self.exec_folder, "samples_generated.csv"
        )
        smiles_df = pd.read_csv(output_folder)
        smiles_list = smiles_df['smiles']
        return smiles_list

    def _clean(self):
        shutil.rmtree(os.path.join(self.exec_folder, "results", self.data_tag))
        shutil.rmtree(os.path.join(self.exec_folder, "results", "data", self.data_tag))
        os.remove(self.data_file)

    def sample(self, smiles_list, n):
        print(self.cwd)
        print(self.exec_folder)
        self._sample(n)
        molecules = self._read_molecules()
        # self._clean()
        return molecules
