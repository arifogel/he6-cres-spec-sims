""" simulation

This module contains a single class (Simulation) that links the simulation blocks together.
One can use the method run_full() to simulate tracks as well as run those tracks through the DAQ,
creating a .spec file. Or one can take a set of downmixed tracks previously created by  run_full()
and saved to a .csv and run them through the DAQ, as it is the calculation of the track properties
(axial_freq, z_max,...) that take the most time.

The general approach is that pandas dataframes, each row describing a single CRES data object (event, segment,
  band, or track), are passed between the blocks, each block adding complexity to the simulation.
 This general structure is broken by the last class (DAQ) which (optionally) creates the binary .spec(k) file
output. This .spec(k) file can then be fed into Katydid just as real data would be.

Classes contained in module:

    * Simulation
    * Results

"""

import time
import logging

import pandas as pd
from numpy import hstack

import he6_cres_spec_sims.simulation_blocks as sim_blocks
import he6_cres_spec_sims.simulation_blocks.config
import he6_cres_spec_sims.simulation_blocks.eventBuilder
import he6_cres_spec_sims.simulation_blocks.trackBuilder
import he6_cres_spec_sims.simulation_blocks.sideBandBuilder
import he6_cres_spec_sims.simulation_blocks.dmTrackBuilder
import he6_cres_spec_sims.simulation_blocks.DAQ

logger = logging.getLogger(__name__)

class Simulation:
    """ Chains together simulation blocks to run full simulation, outputs .csv of Results (defined below)
    """

    def __init__(self, config_path):
        self.config_path = config_path
        self.config = sim_blocks.config.Config(config_path)

    def run_full(self):
        # Initialize all simulation blocks.
        eventbuilder = sim_blocks.eventBuilder.EventBuilder(self.config)
        trackbuilder = sim_blocks.trackBuilder.TrackBuilder(self.config)
        sidebandbuilder = sim_blocks.sideBandBuilder.SideBandBuilder(self.config)
        dmtrackbuilder = sim_blocks.dmTrackBuilder.DMTrackBuilder(self.config)
        if self.config.settings.sim_daq:
            daq = sim_blocks.DAQ.DAQ(self.config)

        # Stage-level timing only -- no computed value or control flow is touched here.
        # Printed as its own clearly-labeled block so it can be grepped straight out of
        # a job's existing log file (see local_spec_sims.py's per-job log_path).
        stage_times = {}

        t0 = time.perf_counter()
        tracks_df = eventbuilder.run()
        stage_times["EventBuilder"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        tracks_df, bands = trackbuilder.run(tracks_df)
        stage_times["TrackBuilder"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        bands = sidebandbuilder.run(tracks_df, bands)
        stage_times["SideBandBuilder"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        downmixed_tracks_df = dmtrackbuilder.run(tracks_df, bands)
        stage_times["DMTrackBuilder"] = time.perf_counter() - t0

        if self.config.settings.sim_daq:
            t0 = time.perf_counter()
            spec_array = daq.run(bands)
            stage_times["DAQ"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        # Save the results of the simulation:
        # For now only write downmixed_tracks to keep things lightweight.
        results = Results(downmixed_tracks_df, bands)
        results.save(self.config_path)
        stage_times["Results.save"] = time.perf_counter() - t0

        total_time = sum(stage_times.values())
        logger.info("\n===== STAGE TIMING BREAKDOWN =====")
        logger.info(f"betas_to_simulate={self.config.physics.betas_to_simulate}, "
              f"events_to_simulate={self.config.physics.events_to_simulate}, "
              f"trapped events (len(tracks_df))={len(tracks_df)}")
        for stage_name, stage_seconds in stage_times.items():
            pct = 100 * stage_seconds / total_time if total_time > 0 else 0
            logger.info(f"  {stage_name:20s} {stage_seconds:9.3f}s  ({pct:5.1f}%)")
        logger.info(f"  {'TOTAL (timed stages)':20s} {total_time:9.3f}s")
        logger.info("===================================\n")

        return None

    def run_daq(self):
        """ Load existing data using Results class (skipping regenerating betas)
        """
        try:
            results = Results.load(self.config_path)
        except Exception as e:
            logger.error("You don't have results to run the daq on.")
            raise e

        # Initialize all necessary simulation blocks.
        daq = sim_blocks.DAQ(self.config)
        specbuilder = sim_blocks.SpecBuilder(self.config, self.config_path)

        # Simulate the action of the DAQ on the loaded dmtracks.
        spec_array = daq.run(results.dmtracks)
        specbuilder.run(spec_array)

        return None

class Results:
    """ Pair of functions (save/ load) that writes the results (currently dmtracks dataFrame)
        to and from a csv with a set name
    """

    def __init__(self, dmtracks, bands=None):
        self.dmtracks = dmtracks
        self.bands = bands

    def get_path_name(self, config_path):
        config_name = config_path.stem
        parent_dir = config_path.parents[0]
        results_dir = parent_dir / "{}".format(config_name)
        return results_dir

    def save(self, config_path):
        results_dict = { "dmtracks": self.dmtracks }

        if self.bands is not None:
            df_bands = pd.DataFrame([band.to_dict() for band in hstack(self.bands) if not band.outside_BW])
            results_dict["bands"] = df_bands

        # First make a results_dir with the same name as the config.
        results_dir = self.get_path_name(config_path)

        # If results_dir doesn't exist, then create it.
        if not results_dir.is_dir():
            results_dir.mkdir()
            logger.info("created directory : %s", results_dir)

        # Now write the results to results_dir:
        for data_name, data in results_dict.items():
            try:
                data.to_csv(results_dir / "{}.csv".format(data_name))
            except Exception as e:
                logger.error("Unable to write {} data.".format(data_name))
                raise e

    def load(self, config_path):
        results_dict = { "dmtracks": None }
        # Load results.
        results_dir = self.get_path_name(config_path)
        for data_name, data in results_dict.items():
            try:
                df = pd.read_csv( results_dir / "{}.csv".format(data_name), index_col=[0])
                results_dict[data_name] = df
            except Exception as e:
                logger.error("Unable to load {} data.".format(data_name))
                raise e

        results = results_dict["dmtracks"]

        return results
