from .base_distribution import BaseDistribution

class RuddDistribution(BaseDistribution):
    """ Generator for a Rudd probability distribution, applicable for distribution of scattering angles
        for inelastic scatters.
        See Eq 25 for PDF: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.44.1644
        Generation done by inverse transform sampling
    """
    def __init__(self, alpha=1.):
        # Set default values
        self.alpha = alpha

    def set_parameters(self, yaml_block):
        # if present, assign from config file
        if "alpha" in yaml_block:
            self.alpha = yaml_block["alpha"]
            self.prefactor = np.sqrt(self.alpha**2 / (1. + self.alpha**2))

    def generate(self, size=None):
        u = self.rng.uniform(0,1,size)
        return np.arctan( self.prefactor * np.tan( np.pi / 2. * u));
