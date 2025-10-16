#!/usr/bin/env python3
import json
import sys
import argparse
from pathlib import Path

import he6_cres_spec_sims as spec_sims


def main():
    """
    A script for running an experiment based on a dictionary input (in 
    the form of a json file). See 
    `/he6-cres-spec-sims/config_files/rocks_exp_config_example.json` for 
    an example of what it needs to contain. 

    Args:
        local_dir (str): where to put experiment results. Ideally on a 
            harddrive.  
        sim_exp_name (str): name of the experiment (identical to name of 
            .json that defines experiment).
    """

    # Parse command line arguments.
    par = argparse.ArgumentParser()
    arg = par.add_argument
    arg(
        "-exp",
        "--exp_dict_path",
        type=str,
        help="path (str) to the pickled dict (.txt) that defines the run conditions of the experiment.",
    )

    args = par.parse_args()
    spec_sims.run_local_experiment(args.exp_dict_path)

    return None

if __name__ == "__main__":
    main()
