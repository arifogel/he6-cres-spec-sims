#!/usr/bin/env python3
import sys
import argparse
import logging
from pathlib import Path
import pytz
import datetime

import he6_cres_spec_sims.simulation as sim

logger = logging.getLogger(__name__)

#sys.path.append("/data/eliza4/he6_cres/simulation/he6-cres-spec-sims")

def main():
    """
    DOCUMENT
    """

    # Parse command line arguments.
    par = argparse.ArgumentParser()
    arg = par.add_argument
    arg(
        "-scp",
        "--sim_config_path",
        type=str,
        help="path (str) to the .yaml that defines the simulation parameters.",
    )

    args = par.parse_args()
    run_simulation(args.sim_config_path)

    return None

def run_simulation(sim_config_path):

    logger.info("START %s. Beginning simulation. Path: %s", get_pst_time(), sim_config_path)

    simulation = sim.Simulation(Path(sim_config_path))
    simulation.run_full()

    logger.info("END %s. Done running simulation. Path: %s", get_pst_time(), sim_config_path)

    return None

def get_pst_time():
    tz = pytz.timezone("US/Pacific")
    return datetime.datetime.now(tz).replace(microsecond=0).replace(tzinfo=None)

if __name__ == "__main__":
    main()
