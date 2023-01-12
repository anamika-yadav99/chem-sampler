from sample import Sampler
import torch
import csv
import pandas as pd
experiment_name = 'ForwardRNN_512'
s = Sampler(experiment_name)
samples = s.sample(N=10, stor_dir='../evaluation', T=0.7, fold=[1], epoch=[9], valid=True, novel=True, unique=True, write_csv=True)
sample_li = []
sample_li.append(samples)
dict_1 = {"smiles": sample_li}
af = pd.DataFrame(dict_1)
af.to_csv('smiles.csv')
