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
        self.CHECKPOINT_DIR = os.path.join(root, "MODEL_DIR" )

    def sample(self, n_molecules):
        cmd = "cd {0}; bash run_sample.sh {1} {2}; cd {3}".format(
            self.exec_folder, self.CHECKPOINT_DIR, n_molecules self.cwd
            )
        SimpleConda().run_commandlines(MOLER_CONDA_ENVIRONMENT, cmd)
