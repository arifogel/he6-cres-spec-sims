import numpy as np
import os
import pandas as pd
import sys

# Additional settings.
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

#path to local imports
#automatically get directory of current file (in config_files)
config_dir = os.path.dirname(os.path.abspath(__file__))
print("spec-sims config directory: "  +str(config_dir))
parent_dir = os.path.dirname(config_dir)
src_dir = parent_dir + "/src/"
print("spec-sims src directory: "  +str(src_dir))
#add to paths
sys.path.append(src_dir)


# Local imports.
import he6_cres_spec_sims.experiment as exp
fields = np.array([1.0])
traps = np.around(fields*1.8/3.25,6)
rand_seeds = np.array(fields*1213, dtype = int)
base_config_path = config_dir + "/example.yaml"

experiment_params = {
    "experiment_name": "ne_051424",
    "base_config_path": base_config_path,
    "events_to_simulate": 10,
    "betas_to_simulate": 10,
    "isotope": "Ne19",
    "rand_seeds": rand_seeds,
    "fields_T" : fields.tolist(),
    "traps_A": traps.tolist()
}

for key, val in experiment_params.items():
    print("{}: {}".format(key, val))

exp.Experiment(experiment_params)
