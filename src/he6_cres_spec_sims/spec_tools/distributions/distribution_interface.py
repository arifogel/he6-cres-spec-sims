from .base_distribution import BaseDistribution

from .beta_decay_distribution import BetaDecayDistribution
from .cauchy_distribution import CauchyDistribution
from .dirac_distribution import DiracDistribution
from .exponential_distribution import ExponentialDistribution
from .normal_distribution import NormalDistribution
from .rudd_distribution import RuddDistribution
from .uniform_distribution import UniformDistribution
from .uniform_annulus_distribution import UniformAnnulusDistribution

from numpy import random

class DistributionInterface:
    """ Interface that connects user config distribution choice to distribution child class
        We do not put on the @singleton decorator (for now)
        It is strongly recommended that you do not create more than 1 instance per simulation
        otherwise, you are running multiple rng's, which will not necessarily be statistically independent/ uncorrelated
    """
    def __init__(self, seed):
        # Set the random number generator (to be shared (passed by reference) between all distributions calling it)
        self.rng = random.default_rng(seed)

    def get_distribution(self, yaml_block):
        # returns child class object (which must have .generate() method) based on config input
        # python <3.10 lacks switches. Used if/else for backwards compatibility
        # Dictionary could have issues with wastefully creating unused class objects

        name = yaml_block["distribution"] # name of distribution from config file

        if name == "beta_decay":
            dist = BetaDecayDistribution()
        elif name == "cauchy" or name=="lorentz":
            dist = CauchyDistribution()
        elif name == "dirac" or name == "fixed":
            dist = DiracDistribution()
        elif name == "exponential":
            dist = ExponentialDistribution()
        elif name == "normal" or name=="gaussian":
            dist = NormalDistribution()
        elif name == "rudd":
            dist = RuddDistribution()
        elif name == "uniform":
            dist = UniformDistribution()
        elif name == "uniform_annulus":
            dist = UniformAnnulusDistribution()
        else:
            raise ValueError("Named distribution not in allowed list! See distribution_interface.py for options")

        dist.set_parameters(yaml_block)
        dist.set_random_engine(self.rng)

        return dist
