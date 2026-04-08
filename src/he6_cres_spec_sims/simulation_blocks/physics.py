import numpy as np

import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc

from he6_cres_spec_sims.constants import *


class Physics:
    """Creates distributions of beta kinematic parameters (position, velocity)
        according to the chosen distributions
    """

    def __init__(self, config):

        self.config = config

        # From the rest of the config file (e.g. freq_acceptance, magnetic_field), we pass along additional parameters to the distributions
        # This is currently only the case for energy distributions. For beta spectra we want restricted ranges, but set semi-automatically
        # These distributions are free to respect these parameters or not

        main_field = self.config.eventbuilder.main_field
        freq_acceptance_high = self.config.physics.freq_acceptance_high
        freq_acceptance_low = self.config.physics.freq_acceptance_low

        self.energy_acceptance_low = sc.freq_to_energy( freq_acceptance_high, main_field)
        self.energy_acceptance_high = sc.freq_to_energy( freq_acceptance_low, main_field)

        self.energy_acceptance_low = np.clip(self.energy_acceptance_low, 0, None)
        self.energy_acceptance_high = np.clip(self.energy_acceptance_high, 0, None)

        if self.energy_acceptance_high <= self.energy_acceptance_low:
            raise ValueError("Energy Acceptance Low must be < Energy Acceptance High! Change Frequency acceptances!")

        # Leave in if manually set
        if self.energy_acceptance_low not in self.config.physics.energy:
            self.config.physics.energy["energy_acceptance_low"] = self.energy_acceptance_low

        if self.energy_acceptance_high not in self.config.physics.energy:
            self.config.physics.energy["energy_acceptance_high"] = self.energy_acceptance_high

        # distribution of energies [eV]
        self.energy_distribution = self.config.dist_interface.get_distribution(self.config.physics.energy)

        ### What to do? Not implemented for other, non-beta distributions
        self.fraction_of_spectrum = 1.
        if hasattr(self.energy_distribution, 'fraction_of_spectrum'):
            self.fraction_of_spectrum = self.energy_distribution.fraction_of_spectrum()

        # distribution of rho positions [m]
        self.rho_distribution = self.config.dist_interface.get_distribution(self.config.physics.rho)

        # distribution of z positions [m]
        self.z_distribution = self.config.dist_interface.get_distribution(self.config.physics.z)


    def generate_beta_energy(self, size=1):
        return self.energy_distribution.generate(size)

    def generate_beta_position_direction(self, size=1):
        """
        Generates a random beta in the trap with pitch angle between
        min_theta and max_theta , and initial position (rho,0,z) between
        min_rho and max_rho and min_z and max_z.
        """

        rho_initial = self.rho_distribution.generate(size)

        # No user choice (for now)
        phi_initial = 2 * PI * self.config.dist_interface.rng.uniform(0, 1, size) * RAD_TO_DEG

        z_initial = self.z_distribution.generate(size)

        min_theta = self.config.physics["min_theta"] / RAD_TO_DEG
        max_theta = self.config.physics["max_theta"] / RAD_TO_DEG

        u_min = (1 - np.cos(min_theta)) / 2
        u_max = (1 - np.cos(max_theta)) / 2

        sphere_theta_initial = np.arccos(1 - 2 * (self.config.dist_interface.rng.uniform(u_min, u_max, size=size))) * RAD_TO_DEG
        sphere_phi_initial = 2 * PI * self.config.dist_interface.rng.uniform(0, 1, size) * RAD_TO_DEG

        # Each position, direction are vectors of vectors
        position = [rho_initial, phi_initial, z_initial]
        direction = [sphere_theta_initial, sphere_phi_initial]

        return position, direction
