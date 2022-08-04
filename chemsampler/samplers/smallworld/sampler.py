import collections
import urllib
import json
from tqdm import tqdm
import time
import random
import pandas as pd
from smallworld_api import SmallWorld
import warnings

warnings.filterwarnings("ignore")


def get_available_maps():
    url = "https://sw.docking.org/search/maps"
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    return data


def get_maps():
    data = get_available_maps()
    labels = ["REAL", "WuXi", "MCule", "Zinc"]
    found_maps = collections.defaultdict(list)
    for k, v in data.items():
        for l in labels:
            l_ = l.lower()
            k_ = k.lower()
            if l_ in k_:
                found_maps[l] += [k]
    found_maps_ = {}
    for k, v in found_maps.items():
        if len(v) == 1:
            v_ = v[0]
            if data[v_]["enabled"] and data[v_]["status"] == "Available":
                found_maps_[k] = v_
        else:
            v_sel = None
            w_sel = None
            for v_ in v:
                if not data[v_]["enabled"] or data[v_]["status"] != "Available":
                    continue
                w_ = data[v_]["numEntries"]
                if v_sel is None:
                    v_sel = v_
                    w_sel = w_
                else:
                    if w_ > w_sel:
                        v_sel = v_
                        w_sel = w_
            if v_sel is not None:
                found_maps_[k] = v_sel
    result = []
    for l in labels:
        if l in found_maps_:
            result += [(l, found_maps_[l])]
    return result


class SmallWorldSampler(object):
    def __init__(self, dist=10, length=100):
        self.maps = get_maps()
        self.sw = SmallWorld()
        self.dist = dist
        self.length = length
        self.seconds_per_query = 3

    def _sample(self, smiles, time_budget_sec):
        t0 = time.time()
        sampled_smiles = []
        for m in self.maps:
            try:
                db_name = m[1]
                results: pd.DataFrame = self.sw.search(
                    smiles, dist=self.dist, db=db_name, length=self.length
                )
            except:
                print(smiles, m, "did not work...")
                results = None
            if results is not None:
                sampled_smiles += list(results["smiles"])
            t1 = time.time()
            if (t1 - t0) > time_budget_sec:
                break
            t0 = time.time()
        return sampled_smiles

    def sample(self, smiles_list, time_budget_sec=600):
        time_budget_sec_per_query = (
            int(time_budget_sec / (self.seconds_per_query * len(smiles_list))) + 1
        )
        sampled_smiles = []
        for smi in tqdm(smiles_list):
            sampled_smiles += self._sample(smi, time_budget_sec_per_query)
        sampled_smiles = list(set(sampled_smiles))
        random.shuffle(sampled_smiles)
        return sampled_smiles
