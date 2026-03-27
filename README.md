
<p align="center"><img width="30%" src="/demo/readme_imgs/he6-cres_logo.png" /></p>

--------------------------------------------------------------------------------
# he6_cres_spec_sims

A package for simulating cres experiments over a variety of magnetic field values.

--------------------------------------------------------------------------------
### Simulate an experiment then make interactive plots of cres track features!

<p align="center"><img width="40%" src="/demo/readme_imgs/plot_stuff.png" />              <img width="50%" src="/demo/readme_imgs/make_plot.png" /></p>

<p align="center"><img width="40%" src="/demo/readme_imgs/plot_stuff_1.png" />              <img width="50%" src="/demo/readme_imgs/make_plot_1.png" /></p>

--------------------------------------------------------------------------------

## Instructions for running simulations on CENPA cluster (rocks): 

* **Get dependencies**: 
	* *Instructions:* 
		* Log on to rocks. 
		* `cd /data/eliza4/he6_cres/simulation/`
		* `pip3 install -r he6-cres-spec-sims/requirements.txt --user` 
	* *Notes:*
		* May need to upgrade pip for the above to work on rocks. 
			* For Winston and I this worked: `pip3 install --upgrade pip`	
		* The following should contain all necessary python packages but if that isn't the case please let Drew Byron know. 
		* Be sure to add the `module load python-3.7.3` to your enviornment setup file or .bash_profile file so that you have access to python3.
		* The above must be done by each user, as it's the current users python packages that the scripts below will be utilizing. 

* **Simulate an experiment**: 
	* *Instructions:* 
		* Log on to rocks. 
		* `cd /data/eliza4/he6_cres/simulation/he6-cres-spec-sims`
		* Set up: 
			* Before running an experiment one needs a `.json` experiment config and a `.yaml` base config to both be in the following directory on rocks: `/data/eliza4/he6_cres/simulation/sim_results/experiments/`. (TODO: Document the config files and the output fields.)
			* Here is how I copy those over from the examples shown in the repo (`he6-cres-spec-sims/config_files`). You should be able to do the same with minimal adjustment of paths.
				* `!scp /home/drew/He6CRES/he6-cres-spec-sims/config_files/rocks* drewbyron@172.25.100.1:/data/eliza4/he6_cres/simulation/sim_results/experiments`
		* Initial run: 
			* `./he6_cres_spec_sims/run_rocks_experiment.py -exp "/data/eliza4/he6_cres/simulation/sim_results/experiments/rocks_exp_config_example.json"`
		* Clean-up: 
			* `./he6_cres_spec_sims/run_rocks_experiment.py -exp "/data/eliza4/he6_cres/simulation/sim_results/experiments/rocks_exp_config_example.json" -clean True` 
	* *Notes:*
		* Initial run:
			* `-exp` (str): Specify the path to the json file that contains the specific attributes (in the form of a python dictionary) of the simulated experiment. See the docstring for the `run_rocks_experiment.py` module for a complete description of all attributes that the `.json` must contain.
			* Say one made a `.json` experiment config locally based on some queries to the he6 postgreSQL database. You could copy that to the rocks `sim_results/experiments` directory with a command like this: 
				* `!scp /media/drew/T7\ Shield/spec_sims_results/rocks_experiments/exp_demo_nov2022.json drewbyron@172.25.100.1:/data/eliza4/he6_cres/simulation/sim_results/experiments`
			* Here `experiment_copies` number of independent (unique random seeds) but otherwise identical experiments are run over rocks. The experiment attribute `experiment_copies` is specified in the `.json` config file. It is parallelized such that each field specified in each copy is sent to a different node. So for example if the `.json` config had these attributes: `{"experiment_copies": 5, "beta_num": 1000, "fields_T": [1.0, 2.0, 3.0]}`, then 5 copies x 3 fields = 15 nodes would each simulate 1000 betas.   
		* Clean up:
			* In the clean-up phase the different copies of the experiment that are produced by the run are combined into one directory that can then be copied onto a local machine for analysis. 
			* In the example used above where we have 5 copies of an experiment spanning 3 fields each with 1000 betas simulated, all of the resultant `.csvs` containing track info for the 5 copies is combined into one directory. 
		* General: 
			* The logs will get very clogged up with different `.txt` files. Occassionally just delete the contents of `simulation/sim_logs`.
	* **Analyzing simulation results**: 
		* *Intructions:* 
			* Use the class `SimulationResults` from the `results.py` module (need to have a local copy of the repo) to grab the experiment results. 
			* Example code to be run in ipynb or within a script: 
				* `sys.path.append("/home/drew/He6CRES/he6-cres-spec-sims/")` (enables you to run this anywhere, change paths to be your own)
				* `from he6_cres_spec_sims.results import SimulationResults`
				* `local_dir = "/media/drew/T7 Shield/spec_sims_results/rocks_experiments"`
				* `sim_exp_name = "rocks_exp_config_example"`
				* `sim = SimulationResults(local_dir = local_dir, sim_exp_name = sim_exp_name )`
			* Explore all that the `sim` instance contains: `sim.results.__dict__`
			* All simulated tracks (as pd.DataFrame) are here: `sim.results.tracks`
		* *Notes:*
			* See the demo notebook for a full illustration of the above: `he6-cres-spec-sims/demo/rocks_sim_experiment_demo.ipynb`.

--------------------------------------------------------------------------------

## Instructions for running simulations locally: 

* **Install**: 
	* *Instructions:*
 		* *Stable Release*
   			* pip install he6_cres_spec_sims
         * *Development Release*  
			* Navigate into desired parent directory.
			* Clone the repo, or `git pull` if you already have it. A hard reset to the remote may be necessary if you have an old version. You should be on the `develop` branch (default).  
				* `git clone git@github.com:Helium6CRES/he6-cres-spec-sims.git`
			* `pip install .` 
	* *Notes:*
		* You will need pip3 for the above to work. 

* **Simulate an experiment**: 
	* *Instructions:* 
		* Navigate into your `he6-cres-spec-sims` parent directory.
		* Set up: 
			* Before running an experiment one needs a `.json` experiment config and a `.yaml` base config to both be in a directory suited for simulation results. For me this is in an external drive. 
			* Here is how I copy those over from the examples shown in the repo (`he6-cres-spec-sims/config_files`). You should be able to do the same with minimal adjustment of paths.
				* `cp /home/drew/He6CRES/he6-cres-spec-sims/config_files/local* /media/drew/T7\ Shield/spec_sims_results/local_experiments`
			* The `base_config_path` field in the `.json` experiment config needs to be manually changed to point at the `.yaml` base config file. Change this path. 
		* Run experiment:
  			* *Script*
				* `./run_local_experiment.py -exp "/media/drew/T7 Shield/spec_sims_results/local_experiments/local_exp_config_example.json"`
     		* *Interactive Python*

         		* `import he6_cres_spec_sims`
           		* `he6_cres_spec_sims.run_local_experiment("/media/drew/T7 Shield/spec_sims_results/local_experiments/local_exp_config_example.json")`  

	* *Notes:*
		* Run experiment:
			* No clean-up is necessary. 
			* The only difference between the config files used locally and on rocks is that there is no `experiment_copies` field in the local `.json` config. The base config (`.yaml`) is identical. 

	* **Analyzing simulation results**: 
		* *Intructions:* 
			* Use the class `ExperimentResults` from the `experiment.py` module to grab the experiment results. 
			* Example code to be run in ipynb or within a script: 
				* `sys.path.append("/home/drew/He6CRES/he6-cres-spec-sims/")` (enables you to run this anywhere, change paths to be your own)
				* `from he6_cres_spec_sims.experiment import ExpResults`
				* `experiment_config_path = "/media/drew/T7 Shield/spec_sims_results/local_experiments/local_exp_config_example/local_exp_config_example_exp.yaml"`
				* `sim = ExpResults.load(experiment_config_path =experiment_config_path  )`
			* Explore all that the `sim` instance contains: `sim.__dict__`
			* All simulated tracks (as pd.DataFrame) are here: `sim.tracks`
		* *Notes:*
			* See the demo notebook for a full illustration of the above: `he6-cres-spec-sims/demo/local_sim_experiment_demo.ipynb`.

--------------------------------------------------------------------------------

## Documentation for simulation config files: 

--------------------------------------------------------------------------------

### Experiment config (.json)

#### Description:

The most basic function of the simulation package is to simulate betas in a given trap depth and main field. An "experiment" is a set of individual simulations run together with organized output. This makes it easy to simulate something closely resembling our experiments; running at a range of trap depths and main field values. The `.json` experiment config serves as the instructions for which simulations to run in the experiment.     
One can find an example of this config file here: `/he6-cres-spec-sims/config_files/rocks_exp_config_example.json`. Note that the name the experiment directory is assigned is the name of the `.json` experiment config file. 

#### An example of it's contents: 

{  
"experiment_copies": 3,    
"base_config_path": "/data/eliza4/he6_cres/simulation/sim_results/experiments/rocks_base_config_example.yaml",   
"isotope": "He6",   
"events_to_simulate": -1,   
"betas_to_simulate": 1e2,   
"rand_seeds": [4062, 3759, 3456, 3153, 2850, 2547, 2244, 1941, 1638, 1335, 1032],   
"fields_T": [3.25, 3.0, 2.75, 2.5, 2.25, 2.0, 1.75, 1.5, 1.25, 1.0, 0.75],   
"traps_A": [1.8, 1.661538, 1.523077, 1.384615, 1.246154, 1.107692, 0.969231, 0.830769, 0.692308, 0.553846, 0.4153845]  
}

#### Explanation of the required fields:

* **experiment_copies:** 
	* Number of copies to be run of the experiment. This field must be omitted if running locally. When running on the cluster the simulations are parallelized such that each field specified in each copy is sent to a different node. So for example if the `.json` config had these attributes: `{"experiment_copies": 5, "beta_num": 1000, "fields_T": [1.0, 2.0, 3.0]}`, then 5 copies x 3 fields = 15 nodes would each simulate 1000 betas. Then the clean-up phase would combine these all into a directory whose name matches the `.json` config name.  

* **base_config_path:**
	* Full path to the base_config. See instructions above for details on where to put the base config. See documentation below for the required fields of the base config. 
* **isotope:**
	* Which isotope to simulate. Currently only works for "Ne19" and "He6" but could easily be extended to other isotopes. 
* **events_to_simulate:** 
	* Number of *trapped events* to simulate. When not set to -1 the simulation will terminate when events_to_simulate betas have been trapped. If set to -1 then the simulation will be terminated once betas_to_simulate have been simulated (trapped or not). One of events_to_simulate or betas_to_simulate should be -1. 
* **betas_to_simulate:** 
	* Number of betas simulated per experiment. See description of events_to_simulate above for details. 
* **rand_seeds:** 
	* List of ints used as the random seed in a simulation; one for each simulation to be run. The lists rand_seeds, fields_T, and traps_A, must all have the same length. 
* **fields_T:**
	* List of magnetic fields to simulate in Tesla. The lists rand_seeds, fields_T, and traps_A, must all have the same length. 
* **traps_A:** 
	* List of trap currents to simulate in Amps. The lists rand_seeds, fields_T, and traps_A, must all have the same length. 

#### To build your own `.json` config: 

* Query the he6-cres database to get real experiment conditions (field, trap depth, isotope). 
* Build a python dictionary: `my_exp_dict` that contains all required fields.
* `import json`
* `json.dump(my_exp_dict, 'path_to_exp_dict.json')`

--------------------------------------------------------------------------------

### Base config (.yaml)

#### Description:

The base configuration `.yaml` file contains the specific parameters to use for each "block" of an individual simulation. When running an experiment, certain settings are over-written (like the `betas_to_simulate`) by the `.json` experiment config. But the general settings defined in the `.yaml` base config are used across all fields and trap depths defined in the `.json` experiment config. 

One can find an example of this config file here: `/he6-cres-spec-sims/config_files/rocks_base_config_example.yaml`

#### An example of it's contents:

* **Settings:**   
	* rand_seed: 234  
* **Physics:** 
	* events_to_simulate : "inf"  
	* betas_to_simulate : 100  
	* energy_spectrum:   
	    * beta_source: "Ne19"  
	* freq_acceptance_high: 18.0e+9  
	* freq_acceptance_low: 19.1e+9  
	* min_rho : 0.0  
	* max_rho : 5.78e-3  
	* min_z : -2.5e-3  
	* max_z : 2.5e-3  
	* min_theta : 89.0  
	* max_theta : 90.0  
* **EventBuilder:**  
	* main_field : 1.700468  
	* trap_current : .1  
	* decay_cell_radius : 5.78e-3  
* **SegmentBuilder:**  
	* mean_track_length : 10.0e-3  
	* jump_size_eV : 14  
	* jump_std_eV : 5  
	* pitch_angle_costheta_std : 0.0  
	* jump_num_max : 0  
* **BandBuilder:**   
	* sideband_num: 1  
	* frac_total_segment_power_cut : 0.01  
	* harmonic_sidebands: True  
	* magnetic_modulation: False   
* **TrackBuilder:**  
	* run_length: 60.0e-6  
* **DMTrackBuilder:**  
	* mixer_freq: 17.9e+9  
* **Daq:**  
	* daq_freqbw: 1.2e+9  
	* freq_bins: 32768  
	* fft_per_slice: 2  
	* band_power_override: 1.0e-16  
	* gain_override: 1.0  
* **SpecBuilder:**  
	* specfile_name: "example_spec_file"  

#### Explanation of the required fields:

 FILL IN. 
--------------------------------------------------------------------------------

### Simulation Blocks
* **Band:**
    * Base class for different types of bands (class wrapper for inter-scatter CRES signals, distinguished from the underlying betas) 
* **Config:**
    * A class used to contain the field map and configurable parameters associated with a given simulation configuration file (for example: config_example.json).
* **Physics:**
    * Creates distributions of beta kinematic parameters (position, velocity) according to the chosen distributions
* **DAQ:**
    * If called, this module  passes through list of downmixed bands through the DAQ, producing fake .spec(k) files. These can be passed through Katydid, identically to data
* **EventBuilder:**
    * Constructs a list of betas which are trapped within the detector volume (doesn’t hit waveguide walls and pitch angle is magnetically trapped)
* **SegmentBuilder:**
    * Constructs a list of betas which are trapped within the detector volume (doesn’t hit waveguide walls and pitch angle is magnetically trapped)
* **TrackBuilder:**
    * Assigns track truth parameters (e.g. start/ end times, frequencies) to created tracks
* **DMTrackBuilder:**
    * Downmixes freq\_start and freq\_stop of simulated tracks to observed frequency band out of DAQ
* **SideBandBuilder:**
    * Creates distributions of beta kinematic parameters (position, velocity) according to the chosen distributions


FILL IN. Need more complete descriptions of what each block of the simulation does/docstrings. 
--------------------------------------------------------------------------------
## Documentation for tracks output: 

Below is a description for all of the features output by the simulations. One can find instructions for generating the `tracks` pd.Dataframe above for a local or cluster based simulation.   

[
'energy':  
'gamma', 
'energy_stop', 
'initial_rho_pos', 
'initial_phi_pos',
'initial_zpos', 
'initial_theta', 
'cos_initial_theta', 
'initial_phi_dir',
'center_theta', 
'cos_center_theta', 
'initial_field', 
'initial_radius',
'center_x', 
'center_y', 
'rho_center', 
'trapped_initial_theta',
'max_radius', 
'min_radius', 
'avg_cycl_freq', 
'b_avg', 
'freq_stop',
'zmax', 
'axial_freq', 
'mod_index', 
'segment_power', 
'slope',
'segment_length', 
'band_power', 
'band_num', 
'segment_num', 
'event_num',
'beta_num', 
'fraction_of_spectrum', 
'energy_accept_high',
'energy_accept_low', 
'gamma_accept_high', 
'gamma_accept_low',
'time_start', 
'time_stop', 
'freq_start', 
'exp_copy', 
'simulation_num',
'field', 
'trap_current'] 


--------------------------------------------------------------------------------

## Important notes: 

* In case of emergency please break the glass. No actually just email me at wbyron@uw.edu. 
* The `develop` branch is tested and working on rocks and locally as of 11/16/22.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

# Development in Progress: 

## To Dos (11/18/22): 

* Make docstrings for all modules/methods! 
* Test a rocks run with a lot of stats to see what breaks.
* Work on documenting the modules (quickly) with docstrings. 
* Then move on to making sure that katydid on rocks still works without the pip install. For now I am uninstalling the package with a pip uninstall. This will break katydid on rocks!! So I need to go fix that once I'm done with this


## Tests: 

* Testing 1e4 betas for 11 fields with 10 copies at 8:49 AM 11/18/22.
	* This worked. We get some betas that fail in that the b_avg is zero, so the integral must fail. Though it resulted in over 5 million rows and 45 columns. 
