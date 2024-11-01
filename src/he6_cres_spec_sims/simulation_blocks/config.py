from pathlib import Path
import yaml

from he6_cres_spec_sims.spec_tools.trap_field_profile import TrapFieldProfile
from he6_cres_spec_sims.spec_tools.distributions.distribution_interface import DistributionInterface

class DotDict(dict):
    """Provides dot.notation access to dictionary attributes."""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class Config:
    """
    TODO: Add a default value for each of these. The dictionary gets overwritten.

    A class used to contain the field map and configurable parameters
    associated with a given simulation configuration file (for example:
    config_example.json).

    ...

    Attributes
    ----------
    simulation, physics, eventbuilder, ... : DotDict
        A dictionary containing the configurable parameters associated
        with a given simulation block. The parameters can be accessed
        with dot.notation. For example eventbuilder.main_field would
        return a field value in T.

    trap_profile: Trap_profile
        An instance of a Trap_profile that corresponds to the main_field
        and trap_strength specified in the config file. Many of the
        spec_tool.spec_calc functions take the trap_profile as a
        parameter.

    field_strength: Trap_profile instance method
        Quick access to field strength values. field_strength(rho,z)=
        field magnitude in T at position (rho,z). Note that there is no
        field variation in phi.

    Methods
    -------
    load_config_file(config_filename)
        Loads the config file.

    load_field_profile()
        Loads the field profile.
    """

    def __init__(self, config_path, load_field=True, daq_only=False):
        """
        Parameters
        ----------
        config_filename: str
            The name of the config file contained in the
            he6_cres_spec_sims/config_files directory.
        """

        # Attributes:
        self.config_path = Path(config_path)
        self.load_field = load_field
        self.daq_only = daq_only

        self.load_config_file()
        if self.load_field:
            self.load_field_profile()

    def load_config_file(self):
        """Loads the YAML config file and creates attributes associated
        with all configurable parameters.

        Parameters
        ----------
        config_filename: str
            The name of the config file contained in the
            he6_cres_spec_sims/config_files directory.

        Raises
        ------
        Exception
            If config file isn't found or can't be opened.
        """

        try:
            with open(self.config_path, "r") as read_file:
                config_dict = yaml.load(read_file, Loader=yaml.FullLoader)

                if self.daq_only:
                    self.daq = DotDict(config_dict["Daq"])

                else:
                    # Take config parameters from config_file.
                    self.settings = DotDict(config_dict["Settings"])
                    self.physics = DotDict(config_dict["Physics"])
                    self.eventbuilder = DotDict(config_dict["EventBuilder"])
                    self.trackbuilder = DotDict(config_dict["TrackBuilder"])
                    self.bandbuilder = DotDict(config_dict["BandBuilder"])
                    self.downmixer = DotDict(config_dict["DMTrackBuilder"])
                    self.daq = DotDict(config_dict["Daq"])

                print("Seed: "+str(self.settings.rand_seed))
                self.dist_interface = DistributionInterface(self.settings.rand_seed)

        except Exception as e:
            print("Config file failed to load.")
            raise e

    def load_field_profile(self):
        """Uses the he6 trap geometry (2021), specified in the
        load_he6_trap module, and the main_field and trap strength
        specified in the config file to create an instance of
        Trap_profile.
        """
        try:
            main_field = self.eventbuilder.main_field
            trap_current = self.eventbuilder.trap_current

            self.trap_profile = TrapFieldProfile(main_field, trap_current)
            self.field_strength = self.trap_profile.field_strength

        except Exception as e:
            print("Field profile failed to load.")
            raise e
