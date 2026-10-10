"""Tests for Physics, using config_files/example.yaml."""

from __future__ import annotations

import unittest

import numpy as np
from python.runfiles import runfiles

from he6_cres_spec_sims.simulation_blocks.config import Config
from he6_cres_spec_sims.simulation_blocks.physics import Physics


def _example_config() -> Config:
    return Config(config_path=runfiles.Create().Rlocation("_main/config_files/example.yaml"))


class PhysicsTest(unittest.TestCase):
    def test_fixed_energy_distribution_returns_the_configured_energy(self) -> None:
        config = _example_config()
        energies = Physics(config).generate_beta_energy(size=3)
        np.testing.assert_allclose(energies, [config.physics.energy["value"]] * 3)

    def test_generated_betas_lie_within_the_configured_bounds(self) -> None:
        config = _example_config()
        (rho, _, z), (theta, _) = Physics(config).generate_beta_position_direction(size=1000)
        tol = 1e-9
        self.assertEqual(len(rho), 1000)
        self.assertTrue(np.all((rho >= config.physics.rho["rho_min"]) & (rho <= config.physics.rho["rho_max"])))
        self.assertTrue(np.all((z >= config.physics.z["low"]) & (z <= config.physics.z["high"])))
        self.assertTrue(np.all((theta >= config.physics.min_theta - tol) & (theta <= config.physics.max_theta + tol)))

    def test_inverted_frequency_acceptance_is_rejected(self) -> None:
        config = _example_config()
        config.physics.freq_acceptance_low, config.physics.freq_acceptance_high = (
            config.physics.freq_acceptance_high,
            config.physics.freq_acceptance_low,
        )
        with self.assertRaises(ValueError):
            Physics(config)


if __name__ == "__main__":
    unittest.main()
