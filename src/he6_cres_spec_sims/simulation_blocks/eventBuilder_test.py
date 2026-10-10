"""Tests for EventBuilder.trap_condition, using config_files/example.yaml (decay_cell_radius = 5.78e-3)."""

from __future__ import annotations

import unittest

import pandas as pd
from python.runfiles import runfiles

from he6_cres_spec_sims.simulation_blocks.config import Config
from he6_cres_spec_sims.simulation_blocks.eventBuilder import EventBuilder


def _event_builder() -> EventBuilder:
    return EventBuilder(Config(config_path=runfiles.Create().Rlocation("_main/config_files/example.yaml")))


def _track(*, initial_theta: float, rho_center: float, max_radius: float) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "initial_theta": [initial_theta],
            "trapped_initial_theta": [87.0],
            "rho_center": [rho_center],
            "max_radius": [max_radius],
            "energy": [20e3],
        }
    )


class TrapConditionTest(unittest.TestCase):
    def test_track_inside_the_cell_at_a_trapped_pitch_angle_is_trapped(self) -> None:
        track = _track(initial_theta=90.0, rho_center=0.0, max_radius=0.0)
        self.assertTrue(_event_builder().trap_condition(track))

    def test_track_below_the_trapped_pitch_angle_is_not_trapped(self) -> None:
        track = _track(initial_theta=86.0, rho_center=0.0, max_radius=0.0)
        self.assertFalse(_event_builder().trap_condition(track))

    def test_track_reaching_past_the_decay_cell_wall_is_not_trapped(self) -> None:
        track = _track(initial_theta=90.0, rho_center=0.005, max_radius=0.001)
        self.assertFalse(_event_builder().trap_condition(track))


if __name__ == "__main__":
    unittest.main()
