import os
import pandas as pd
from configparser import ConfigParser
from tqdm import tqdm
from ast import literal_eval
from glob import glob
import multiprocessing as mp

def load_summary_samples(result_paths, all_results=None, force=False,
                         data_key="samples"):
    all_results = all_results or {}
    for run, rundir in result_paths.items():
        if run not in all_results or force:
            path = os.path.join(rundir, "summary_samples.hdf5")
            if os.path.exists(path):
                all_results[run] = pd.read_hdf(path, data_key)
                all_results[run]["run"] = run
            else:
                print(f"WARNING: did not find {path}")
    return all_results

def load_summary_stats(result_paths, amp_stats=None, force=False, ci=0.68,
                       n=1):
    amp_stats = amp_stats or {}
    for run, rundir in result_paths.items():
        if run not in amp_stats or force:
            path = os.path.join(rundir, f"summary_stats_A{n}_ci{ci}.dat")
            if os.path.exists(path):
                amp_stats[run] = pd.read_csv(path, sep=' ')
                amp_stats[run]["run"] = run
            else:
                print(f"WARNING: did not find {path}")
    return amp_stats

def try_parse(x): 
    try: 
        return float(x) 
    except (TypeError,ValueError):
        try:
            return literal_eval(x)
        except (TypeError,ValueError,SyntaxError):
            return x
        
def _get_injection(p):
    config = ConfigParser()
    config.read(p)
    return {k: try_parse(v) for k,v in config['injection'].items()}

def get_injections(result_paths, nproc=8):
    injkws = None
    for run, rundir in tqdm(result_paths.items()):
        paths = glob(os.path.join(rundir, '*/engine/*/config.ini'))
        with mp.Pool(nproc) as pool:
            injs = pool.map(_get_injection, paths)
        for p, ikws in zip(paths, injs):
            injkws = injkws or ikws
            if injkws != ikws and ikws is not None:
                print(f"WARNING: different injections! {p}")
                print(ikws)
                break
    return injkws
