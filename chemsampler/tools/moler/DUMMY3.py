#!/usr/bin/env python3
from .molecule_generation.wrapper import load_model_from_directory
from .molecule_generation.utils.cli_utils import (
    setup_logging,
    supress_tensorflow_warnings,
)
import os
import shutil
import configparser

from ...utils.conda import SimpleConda

root = os.path.dirname(os.path.abspath(__file__))

MOLLIB_CONDA_ENVIRONMENT = "moler"
DATA_TAG = "my_molecules"


class MoLerSampler:
    def __init__(self):
        self.cwd = os.getcwd()
        #path here
        self.model_dir = './molecule_generation.MODEL_DIR'
        self.exec_folder = os.path.join(root, "", "experiments")
        self.data_tag = DATA_TAG
        self.data_file_relative = os.path.join("..", "data", self.data_tag + ".txt")
        self.data_file = os.path.join(self.exec_folder, self.data_file_relative)

    def sample(self, num_samples: int):
        supress_tensorflow_warnings()
        setup_logging()
        print('task started')
        with load_model_from_directory(self.model_dir) as model:
            samples = model.sample(num_samples)
        print("\n".join(samples))

        
