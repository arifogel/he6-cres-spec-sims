from abc import ABC, abstractmethod
import numpy as np
from scipy.interpolate import interp1d

class BaseDistribution(ABC):
    """ Abstract Base Class for random distribution generator. Child classes need to have generate()
    """
    def set_random_engine(self, rng):
       self.rng = rng

    def inverse_cdf_helper(self, pdf, xMin, xMax, N, **kwargs):
        # Function which takes in a PDF, an xRange over which the PDF is defined, and an interpolation factor
        # Returns a function f(u), such that the output samples according to the PDF between the set x-values
        # Ideally, do this once in setting up your random number generator and not at each f(u) call
        x = np.linspace(xMin, xMax, int(N))
        cumpdf = np.cumsum(pdf(x, **kwargs))
        cumpdf -= cumpdf[0]
        cumpdf *= 1. / cumpdf[-1]
        f = interp1d(cumpdf, x)
        return f

    @abstractmethod
    def set_parameters(self, yaml_block):
        raise NotImplementedError("set_parameters() must be implemented")

    @abstractmethod
    def generate(self, size):
        raise NotImplementedError("generate() must be implemented")
