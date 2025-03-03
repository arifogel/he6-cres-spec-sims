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

import pandas as pd

import he6_cres_spec_sims.simulation_blocks as sim_blocks
import he6_cres_spec_sims.simulation_blocks.config
import he6_cres_spec_sims.simulation_blocks.eventBuilder
import he6_cres_spec_sims.simulation_blocks.trackBuilder
import he6_cres_spec_sims.simulation_blocks.sideBandBuilder
import he6_cres_spec_sims.simulation_blocks.dmTrackBuilder
import he6_cres_spec_sims.simulation_blocks.DAQ

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

        tracks_df = eventbuilder.run()
        tracks_df, bands = trackbuilder.run(tracks_df)
        tracks_df, bands = sidebandbuilder.run(tracks_df, bands)
        downmixed_tracks_df = dmtrackbuilder.run(tracks_df, bands)
        print(bands)

        if self.config.settings.sim_daq:
            spec_array = daq.run(downmixed_tracks_df)
        # Save the results of the simulation:
        # For now only write downmixed_tracks to keep things lightweight.
        results = Results(downmixed_tracks_df)
        results.save(self.config_path)

        return None

    def run_daq(self):
        """ Load existing data using Results class (skipping regenerating betas)
        """
        try:
            results = Results.load(self.config_path)
        except Exception as e:
            print("You don't have results to run the daq on.")
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

    def __init__(self, dmtracks):
        self.dmtracks = dmtracks

    def get_path_name(self, config_path):
        config_name = config_path.stem
        parent_dir = config_path.parents[0]
        results_dir = parent_dir / "{}".format(config_name)
        return results_dir

    def save(self, config_path):
        # Only writing these dmtracks to make the simulations more lightweight
        results_dict = { "dmtracks": self.dmtracks }

        # First make a results_dir with the same name as the config.
        results_dir = self.get_path_name(config_path)

        # If results_dir doesn't exist, then create it.
        if not results_dir.is_dir():
            results_dir.mkdir()
            print("created directory : ", results_dir)

        # Now write the results to results_dir:
        for data_name, data in results_dict.items():
            try:
                data.to_csv(results_dir / "{}.csv".format(data_name))
            except Exception as e:
                print("Unable to write {} data.".format(data_name))
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
                print("Unable to load {} data.".format(data_name))
                raise e

        results = results_dict["dmtracks"]

        return results
