from .base_distribution import BaseDistribution

class NormalDistribution(BaseDistribution):
    """ Generator for a normal (Gaussian) probability distribution
    """
    def __init__(self, mean=0, sigma=1):
        self.mean = mean
        self.sigma = sigma

    def set_parameters(self, yaml_block):
        # if present, assign from config file
        if "mean" in yaml_block:
            self.mean = yaml_block["mean"]
        if "sigma" in yaml_block:
            self.sigma = yaml_block["sigma"]

    def generate(self, size=None):
        return self.rng.normal(self.mean, self.sigma, size)
