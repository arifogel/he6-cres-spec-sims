"""Tests for Config, loading config_files/example.yaml."""

from __future__ import annotations

import unittest

from python.runfiles import runfiles

from he6_cres_spec_sims.simulation_blocks.config import Config


def _example_config() -> Config:
    return Config(config_path=runfiles.Create().Rlocation("_main/config_files/example.yaml"))


class ConfigTest(unittest.TestCase):
    def test_sections_are_read_from_the_yaml(self) -> None:
        config = _example_config()
        self.assertEqual(config.settings.rand_seed, 716)
        self.assertEqual(config.physics.events_to_simulate, 100)
        self.assertEqual(config.eventbuilder.main_field, 0.75)
        self.assertEqual(config.eventbuilder.trap_current, 0.186)

    def test_field_at_the_trap_center_is_close_to_the_main_field(self) -> None:
        config = _example_config()
        main_field = config.eventbuilder.main_field
        self.assertAlmostEqual(config.field_strength(0, 0), main_field, delta=0.01 * main_field)

    def test_field_at_the_trap_center_matches_the_known_value(self) -> None:
        # Pins the trap field model's output for the example config's main field and trap current.
        self.assertAlmostEqual(_example_config().field_strength(0, 0), 0.749255, places=6)


if __name__ == "__main__":
    unittest.main()
