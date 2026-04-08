from .base_distribution import BaseDistribution
import numpy as np
import scipy.stats as stats

class AseevDistribution(BaseDistribution):
    """ Generator for an Aseev-like energy loss probability distribution for inelastic scatters.
        See Eq 8 for the PDF: https://link.springer.com/article/10.1007/s100530050525
        Represents a truncated Gaussian distribution between [0, eps_c] and a truncated Cauchy distribution between [eps_c, inf) 
        Found by inverse transform sampling. Use "standard" Gaussian, Cauchy parameterizations, instead of Aseev parameterization
        Describes energy loss distributions from inelastic scattering in Katrin
        TODO: pre-compute functions that don't need to be invoked for each call (lines above p_gauss?). Allow for multiple isotopes
    """
    def __init__(self, isotope="H2"):
        # Set default values
        self.isotope = isotope

    def set_parameters(self, yaml_block):
        # if present, assign from config file
        if "isotope" in yaml_block:
            self.isotope = yaml_block["isotope"]

        self.x_c = 14.12
        self.mu_1 = 12.6
        self.sigma_1 = 1.85/2
        self.mu_2 = 14.3
        self.gamma_2 = 12.5/2.

    def generate(self, size=None):
        #Compute the probability of picking from each distribution

        pdf_gauss = lambda x: np.exp(-(x-self.mu_1)**2 / (2*self.sigma_1**2)) / np.sqrt(2*np.pi*self.sigma_1**2)
        pdf_cauchy = lambda x: 1./(np.pi * self.gamma_2 * (1. + ((x-self.mu_2)/self.gamma_2)**2))

        cdf_gauss = lambda x: stats.norm.cdf(x, loc=self.mu_1, scale=self.sigma_1)
        cdf_cauchy = lambda x: np.arctan( (x - self.mu_2) / self.gamma_2)/np.pi + 0.5

        #CDF ranges for truncated regions of gaussian, cauchy distributions
        u_range_gauss = np.array([cdf_gauss(0), cdf_gauss(self.x_c)])
        u_range_cauchy = np.array([cdf_cauchy(self.x_c), 1])

        delta_u_gauss = u_range_gauss[1] - u_range_gauss[0]
        delta_u_cauchy = u_range_cauchy[1] - u_range_cauchy[0]

        #For the piecewise PDF, probability that a sample is <x_c, and is given by the truncated Gaussian distribution
        p_gauss = delta_u_gauss / (delta_u_gauss + pdf_gauss(self.x_c) / pdf_cauchy(self.x_c) * delta_u_cauchy)
        #print("p_gauss: ",p_gauss)

        #Sample from cauchy for x > x_c (for all points, first)
        u = self.rng.uniform(u_range_cauchy[0], u_range_cauchy[1], size)
        samples = self.mu_2 + self.gamma_2 * np.tan(np.pi * (u - 0.5))

        #Randomly assign length N vector 0 or 1. 0's generate truncated Cauchy dist [x_c, inf), 1's generate truncated Gaussian [0,x_c)
        #Overwrite samples with truncated Gaussian [0, x']
        dist_choice = (self.rng.binomial(1, p_gauss, size=size)==1)
        if hasattr(dist_choice, "__len__"):
            n_gauss = sum(dist_choice)
            u = np.random.uniform(u_range_gauss[0], u_range_gauss[1], n_gauss)
            samples[dist_choice] = stats.norm.ppf(u, loc=self.mu_1, scale=self.sigma_1)
        elif dist_choice:
            u = np.random.uniform(u_range_gauss[0], u_range_gauss[1])
            samples = stats.norm.ppf(u, loc=self.mu_1, scale=self.sigma_1)

        return samples

