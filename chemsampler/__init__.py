import os
import csv

root = os.path.abspath(os.path.dirname(__file__))

def example():
    smiles_list = []
    with open(os.path.join(root, "data", "known_hits.csv"), "r") as f:
        reader = csv.reader(f)
        for r in reader:
            smiles_list += [r[0]]
    return smiles_list