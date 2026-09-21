#!/usr/bin/env python3
import sys
import argparse
import logging
from pathlib import Path
import pytz
import datetime

import he6_cres_spec_sims.simulation as sim
from he6_cres_spec_sims.logging_setup import init_logging

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
    arg(
        "--log-level",
        dest="log_level",
        type=str,
        default="INFO",
        help="root log level (e.g. DEBUG, INFO, WARNING) -- see logging_setup.init_logging",
    )
    arg(
        "--log-override",
        dest="log_override",
        type=str,
        default=None,
        help=(
            "comma-separated logger_name=LEVEL overrides for individual loggers "
            "(e.g. 'he6_cres_spec_sims.simulation_blocks.DAQ=DEBUG') -- see logging_setup.init_logging"
        ),
    )

    args = par.parse_args()
    init_logging(args.log_level, args.log_override)
    run_simulation(args.sim_config_path)

    return None

def run_simulation(sim_config_path):

    logger.info(f"START. Current (PST) time: {get_pst_time()}")
    logger.info(f"\n\n\n Beginning simulation. Path: {sim_config_path} \n\n\n")

    simulation = sim.Simulation(Path(sim_config_path))
    simulation.run_full()

    logger.info(f"\n\n\n Done running simulation. Path: {sim_config_path}\n\n\n")
    logger.info(f"END. Current (PST) time: {get_pst_time()}")
    
    return None

def get_pst_time():
    tz = pytz.timezone("US/Pacific")
    return datetime.datetime.now(tz).replace(microsecond=0).replace(tzinfo=None)

if __name__ == "__main__":
    main()
