from .base_distribution import BaseDistribution

class CauchyDistribution(BaseDistribution):
    """ Generator for a Cauchy (Lorentz) probability distribution
    """
    def __init__(self, mean=0, gamma=1):
        self.mean = mean
        self.gamma = gamma #half-width half max

    def set_parameters(self, yaml_block):
        # if present, assign from config file
        if "mean" in yaml_block:
            self.mean = yaml_block["mean"]
        if "gamma" in yaml_block:
            self.gamma = yaml_block["gamma"]

    def generate(self, size=None):
        #found by inverse transform sampling analytically
        u = self.rng.uniform(size=size)
        return self.gamma * (self.mean  + np.tan(np.pi * (u - 0.5)))
