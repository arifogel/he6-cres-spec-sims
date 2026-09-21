import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.integrate import quad
from scipy.interpolate import interp1d

from .base_distribution import BaseDistribution

logger = logging.getLogger(__name__)

import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc
from he6_cres_spec_sims.constants import *

class BetaDecayDistribution(BaseDistribution):
    """A class used to produce and interact with the pdf of a beta
    spectrum. Note that the pickled beta spectra are for b = 0. The
    distortion due to b and the necessary renormalization is done here.
    Unnormalized beta spectrum (dNdE_unnormed_SM) based on interpolating
    the cobs_bs look-up table (df). Drew (March 2023) found that 10^4 vs
    10^3 pkled points was a 10^-6 effect and quadratic vs cubic was a
    10^-8 effect. So he went with 10^3 pts and a quadratic interpolation
    (for performance).
    N. Buzinsky after D. Byron (main author)

    Args:
        isotope (str): Isotope to build BetaSpectrum class for. Currently
            the two options are: "Ne19", "He6". More spectra can be pickled
            based on thecobs (Leendert Hayen's beta spectrum library)
            and put in the 'pickled_spectra' directory. Instructions for
            how to pickle the spectra in the correct format can be found
            ____. TODO: WHERE??? FILL IN LATER..
        b (float): Value for the Fierz Inerference coefficient (little b).
    """

    def __init__(self, isotope="He6", b=0, E_min=0, E_max=1):
        # Include all possible isotopes here.
        self.allowed_isotopes = {
            "Ne19": {"W_max": 5.337377690802349, "Z": 9, "A": 19},
            "He6": {"W_max": 7.864593361688904, "Z": 2, "A": 6},
        }

        self.isotope = isotope
        self.b = b
        self.E_min = E_min
        self.E_max = E_max

    def set_parameters(self, yaml_block):
        # if present, assign from config file
        if "isotope" in yaml_block:
            self.isotope = yaml_block["isotope"]
        if "b" in yaml_block:
            self.b = yaml_block["b"]
        if "energy_acceptance_low" in yaml_block:
            self.E_min = yaml_block["energy_acceptance_low"]
        if "energy_acceptance_high" in yaml_block:
            self.E_max = yaml_block["energy_acceptance_high"]

        self.load_beta_spectrum()

        # SM spectrum (with all corrections) from thecobs
        # 'extrapolate' will lead to bad behaviour outside of the physical gamma values (1, W0).

        self.dNdE_unnormed_SM = interp1d(
            self.cobs_bs["W"],
            self.cobs_bs["Spectrum"],
            kind="quadratic",
            bounds_error=False,
            fill_value=0,
        )

        self.W_min = 1. + self.E_min / ME
        self.W_max = 1. + self.E_max / ME
 
        nEnergyIntervals = 10**5
        # Using known PDF, use base class helper for inverse transform sampling inverse CDF
        self.beta_decay_inv_cdf = self.inverse_cdf_helper(self.dNdE, self.W_min, self.W_max, nEnergyIntervals)

    def load_beta_spectrum(self):
        # Key for relative paths from executing file is to use __file__
        pkl_dir = Path(__file__).parent.resolve() / Path("Data")
        self.pkl_spectra_path = pkl_dir / Path(f"{self.isotope}.json")

        with open(self.pkl_spectra_path) as json_file:
            isotope_info = json.load(json_file)
            #print(isotope_info)

        self.cobs_bs = pd.DataFrame.from_dict(isotope_info.pop("cobs_df"))


    def dNdE(self, W):
        """
        Add in the little b distortion. Unnormalized PDF
        """
        return np.clip(self.dNdE_unnormed_SM(W) * (1 + (self.b / W)), 0, np.inf)

    def fraction_of_spectrum(self):
        logger.debug("W_min=%s W_max=%s", self.W_min, self.W_max)
        spectrum_BW, norm_err = quad( self.dNdE, self.W_min, self.W_max,)
        spectrum_total, norm_err = quad( self.dNdE, 1, self.allowed_isotopes[self.isotope]["W_max"])
        return spectrum_BW / spectrum_total

    def generate(self, size=None):
        """Generate N random samples from dNdE(E) between E_start and E_stop
        """
        # Inverse transform sampling, get N random values between 0 to 1
        rand_cdf_values = self.rng.uniform(size=size)
        rand_Ws = self.beta_decay_inv_cdf(rand_cdf_values)
        rand_Es = ME * (rand_Ws - 1.)
        return rand_Es
