from dataclasses import dataclass
import json
from natsort import natsorted
import numpy as np
import pandas as pd
from pathlib import Path
from shutil import copyfile, rmtree
from typing import List
import yaml

import he6_cres_spec_sims.simulation as sim

import he6_cres_spec_sims.simulation_blocks as sim_blocks
import he6_cres_spec_sims.simulation_blocks.config

#this function runs everything, previously in run_local_experiment.py script,
# putting it here allows you to more easily run experiments interactively
def run_local_experiment(dict_path):

    print(f"\n\n\n Beginning local simulation.\n\n\n")

    experiment_name = Path(dict_path).stem

    sim_experiment_params = json.load(open(dict_path))
    sim_experiment_params["experiment_name"] = experiment_name

    for key, val in sim_experiment_params.items():
        print("{}: {}".format(key, val))

    sim_experiment = Experiment(sim_experiment_params)

    print(f"\n\n\n Done running simulation. {sim_experiment_params}")

    return None

# Utility function:
def get_experiment_dir(experiment_params: dict) -> Path:
    base_config_path = Path(experiment_params["base_config_path"])
    parent_dir = base_config_path.parents[0]
    if "output_base_path" in experiment_params:
        parent_dir = Path(experiment_params["output_base_path"])
    experiment_name = experiment_params["experiment_name"]
    experiment_dir = parent_dir / experiment_name
    return experiment_dir

# Utility function:
def get_config_paths(experiment_params: dict) -> List[Path]:

    experiment_dir = get_experiment_dir(experiment_params)

    suffix = "T.yaml"
    config_paths = [
        x
        for x in experiment_dir.glob("**/*{}".format(suffix))
        if (x.is_file() and ("ipynb" not in str(x)))
    ]

    if len(config_paths) == 0:
        raise ValueError(
            "No config files found in experiment_dir: {}. \n\
            Have you run the experiment?\
            (run: exp = exp.Experiment(experiment_params))".format(
                experiment_dir
            )
        )

    # The natsort ensures that you get a good ordering in terms of consecutive integers
    # within the config file name.
    config_paths = natsorted(config_paths, key=str)

    return config_paths


class Experiment:
    def __init__(self, experiment_params: dict) -> None:

        self.experiment_params = experiment_params

        self.run_experiment(self.experiment_params)

    def run_experiment(self, experiment_params: dict) -> None:

        # Make all the config files specified by the experiment dictionary.
        self.create_configs_for_experiment(experiment_params)

        # Make the experiment config file.
        self.create_experiment_config_file(experiment_params)

        # Collect all the config path names.
        self.config_paths = get_config_paths(experiment_params)

        # Run all of those simulations.
        self.run_sims(self.config_paths)

        return None

    def create_configs_for_experiment(self, experiment_params: dict) -> None:

        base_config_path = Path(experiment_params["base_config_path"])
        experiment_dir = get_experiment_dir(experiment_params)

        # Make the experiments_dir if it doesn't exist. May want to delete contents here?
        if not experiment_dir.is_dir():
            experiment_dir.mkdir()
            print("Created directory: {} ".format(experiment_dir))

        else:
            raise ValueError("Directory already exists: {} ".format(experiment_dir))

        # Grab data from experiment_params dict.
        isotope = experiment_params["isotope"]
        events_to_simulate = experiment_params["events_to_simulate"]
        betas_to_simulate = experiment_params["betas_to_simulate"]
        seeds = experiment_params["rand_seeds"]
        fields = experiment_params["fields_T"]
        traps = experiment_params["traps_A"]

        for i, (seed, field, trap) in enumerate(zip(seeds, fields, traps)):

            # Round the field because there are often small rounding errors.
            field = np.around(field, 6)
            config_path = experiment_dir / "{}_field_{}T{}".format(
                i, field, base_config_path.suffix
            )

            copyfile(base_config_path, config_path)

            # Open the config file and grab the contents.
            with open(config_path, "r") as f:
                config_dict = yaml.load(f, Loader=yaml.FullLoader)

            # Make the appropriate alterations to the config_dict
            # For seed = None, rng pulls from hardware entropy
            if seed is not None:
                config_dict["Settings"]["rand_seed"] = int(seed)
            else:
                config_dict["Settings"]["rand_seed"] = None
            config_dict["Physics"]["events_to_simulate"] = int(events_to_simulate)
            config_dict["Physics"]["betas_to_simulate"] = int(betas_to_simulate)
            #config_dict["Physics"]["energy"]["isotope"] = str(isotope)
            config_dict["EventBuilder"]["main_field"] = float(field)
            config_dict["EventBuilder"]["trap_current"] = float(trap)

            with open(config_path, "w") as f:
                yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

        return None

    def create_experiment_config_file(self, experiment_params: dict) -> None:

        experiment_dir = get_experiment_dir(experiment_params)

        # Define path for experiment config.
        experiment_config = experiment_dir / ( experiment_params["experiment_name"] + "_exp.yaml")

        with open(experiment_config, "w") as f:
            yaml.dump(experiment_params, f, default_flow_style=False, sort_keys=False)

        return None

    def run_sims(self, config_paths: List[Path]) -> None:

        for i, config_path in enumerate(config_paths):
            print("+++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
            print("Running simulation {} / {}\n\n".format(i, len(config_paths)))
            print("+++++++++++++++++++++++++++++++++++++++++++++++++")
            simulation = sim.Simulation(config_path)
            simulation.run_full()

        return None


@dataclass
class ExpResults:

    experiment_params: dict
    # base_config: object
    config_paths: List[Path]
    tracks: pd.DataFrame

    @classmethod
    def load(cls, experiment_config_path: str = None):

        # Path to the exp_config.yaml file.
        experiment_config_path = Path(experiment_config_path)

        # Open the config file and grab the contents.
        with open(experiment_config_path, "r") as f:
            experiment_params = yaml.load(f, Loader=yaml.FullLoader)

        if experiment_config_path is None:

            raise ValueError("Need to provide experiment_config_path.")

        exp_results_dict = {
            "experiment_params": experiment_params,
            # "base_config": Config(experiment_params["base_config_path"]),
            "config_paths": None,
            "tracks": None,
        }

        # Then collect all the config path names.
        # Note that the path in the experiment_params may be misleading if the
        # experiment was run on rocks.
        config_paths = get_config_paths_results(experiment_config_path)
        exp_results_dict["config_paths"] = config_paths

        tracks_list = []
        fields = []

        # Figure out how many betas were sampled; depends on the mode (beta_num or event_num).
        # Note that if this is -1 then you get the entire array (event_mode).
        beta_num = experiment_params["betas_to_simulate"]

        for i, config_path in enumerate(config_paths):

            print("+++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
            print("Loading simulation {} / {}\n\n".format(i, len(config_paths)))
            print("+++++++++++++++++++++++++++++++++++++++++++++++++")

            # Get the simulation parameters from the config.
            config = sim_blocks.config.Config(config_path)
            field = config.eventbuilder.main_field
            trap_current = config.eventbuilder.trap_current
            print("\nSet field: {}, Trap current: {}\n".format(field, trap_current))
            results = sim.Results.load(config_path)
            tracks = results.dmtracks
            tracks["simulation_num"] = i
            tracks["field"] = field
            tracks["trap_current"] = trap_current
            tracks_list.append(tracks)

        exp_results_dict["tracks"] = pd.concat(tracks_list)

        exp_results = cls(
            exp_results_dict["experiment_params"],
            # exp_results_dict["base_config"],
            exp_results_dict["config_paths"],
            exp_results_dict["tracks"],
        )

        return exp_results


def get_config_paths_results(experiment_config_path: Path) -> List[Path]:

    experiment_dir = experiment_config_path.parents[0]
    print(experiment_dir)

    suffix = "T.yaml"
    config_paths = [
        x
        for x in experiment_dir.glob("**/*{}".format(suffix))
        if (x.is_file() and ("ipynb" not in str(x)))
    ]

    if len(config_paths) == 0:
        raise ValueError(
            "No config files found in experiment_dir: {}. \n\
            Have you run the experiment?\
            (run: exp = exp.Experiment(experiment_params))".format(
                experiment_dir
            )
        )

    # The natsort ensures that you get a good ordering in terms of consecutive integers
    # within the config file name.
    config_paths = natsorted(config_paths, key=str)

    return config_paths
