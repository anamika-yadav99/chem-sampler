import csv
import os
import tempfile
import subprocess
import shutil

root = os.path.abspath(os.path.dirname(__file__))


class Sascorer(object):
    def __init__(self):
        self.framework_dir = os.path.join(root, "framework")

    def score(self, smiles_list):
        tmp_dir = tempfile.mkdtemp(prefix="chemsampler")
        input_file = os.path.join(tmp_dir, "input.csv")
        output_file = os.path.join(tmp_dir, "output.csv")
        log_file = os.path.join(tmp_dir, "log.log")
        with open(input_file, "w") as f:
            writer = csv.writer(f)
            writer.writerow(["smiles"])
            for smi in smiles_list:
                writer.writerow([smi])
        run_file = os.path.join(tmp_dir, "run.sh")
        with open(run_file, "w") as f:
            lines = [
                "bash {0}/run_predict.sh {0} {1} {2}".format(
                    self.framework_dir, input_file, output_file
                )
            ]
            f.write(os.linesep.join(lines))
        cmd = "bash {0}".format(run_file)
        with open(log_file, "w") as fp:
            subprocess.Popen(
                cmd, stdout=fp, stderr=fp, shell=True, env=os.environ
            ).wait()
        with open(output_file, "r") as f:
            reader = csv.reader(f)
            next(reader)
            scores = []
            for r in reader:
                scores += [float(r[0])]
        assert len(scores) == len(smiles_list)
        shutil.rmtree(tmp_dir)
        return scores
