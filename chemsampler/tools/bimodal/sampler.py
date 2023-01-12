from ...utils.conda import SimpleConda
import os
root = os.path.dirname(os.path.abspath(__file__))

BIMODAL_CONDA_ENVIRONMENT = "bimodal"

class _fit:
    def __init__(self):
        self.cwd = os.getcwd()
        self.exec_folder = os.path.join(root)


        
    def sample(self, n):
        cmd = "cd {0} ; python _sampler.py {1}; cd {3}".format(
            self.exec_folder, n ,self.cwd
        )
        SimpleConda().run_commandlines(BIMODAL_CONDA_ENVIRONMENT, cmd)
    
        
