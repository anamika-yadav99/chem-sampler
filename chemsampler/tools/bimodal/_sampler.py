import os
import shutil
import configparser
from model.sample import Sampler
import sys

root = os.path.dirname(os.path.abspath(__file__))



class _sampler:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)
        self.temperature = 0.7
        self.fold = [1]
        self.epoch = [9]
        self.store_dir = os.path.join(root, 'evaluation')

    def sampler_ForwardRNN(self, n):
        experiment_name = 'ForwardRNN_512_fixed'
        s = Sampler(experiment_name)
        s.sample(N=n, stor_dir= self.store_dir, T=self.temperature, fold=self.fold, epoch= self.epoch,
                 valid=True, novel=True, unique=True, write_csv=False)

    # def sampler_BIMODAL(self, n):
    #     experiment_name = 'BIMODAL_random_1024'
    #     s = Sampler(experiment_name)
    #     s.sample(N=n, stor_dir= self.store_dir, T=self.temperature, fold=self.fold, epoch= self.epoch,
    #              valid=True, novel=True, unique=True, write_csv=False)

if __name__ == "__main__":
    n = sys.arg[1]
    _sampler.sampler_ForwardRNN(n)

