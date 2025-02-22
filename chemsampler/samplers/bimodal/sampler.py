from ...tools.bimodal.sampler import BimodalSampler
from ...tools.fpsim2.searcher import SimilaritySearcher, RandomSearcher
import os

class MolerSampler:
    def __init__(self):
        #TODO: Come up with a better way to store precalculated files
        self.fp_filename = '../chemsampler/data/pre_calculated/bimodal/bimodal_10k.h5'
        self.db_smiles_filename = '../chemsampler/data/pre_calculated/bimodal/bimodal_10k.csv'
        
        if os.path.exists(self.fp_filename) is False:
            smiles_list = SimilaritySearcher.read_db_smiles()
            SimilaritySearcher(self.fp_filename).fit(smiles_list)

    def Sampler(self, smiles_list, n,  search_pre_calculated=True, cutoff = 0.7):
        if search_pre_calculated==True:
            samples = []
            for smile in smiles_list:
                _samples = SimilaritySearcher(self.fp_filename).search(smile, cutoff=0.7)
                if len(_samples) == 0:
                    print('Coudn\'t find samples in the cutoff range : {}. Generating random samples:'.format(cutoff))
                    _samples = RandomSearcher(self.db_smiles_filename).search(n)
                samples += _samples 
        
        else:
            sampler = BimodalSampler()
            samples = sampler.sample(n)

        return samples

