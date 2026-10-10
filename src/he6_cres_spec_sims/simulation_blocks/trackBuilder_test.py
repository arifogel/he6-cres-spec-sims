"""Tests for TrackBuilder, on fixed distributions so every scatter is deterministic."""

from __future__ import annotations

import unittest

import pandas as pd
from python.runfiles import runfiles

from he6_cres_spec_sims.simulation_blocks.config import Config
from he6_cres_spec_sims.simulation_blocks.eventBuilder import EventBuilder
from he6_cres_spec_sims.simulation_blocks.trackBuilder import TrackBuilder


def _build(jump_num_max: int) -> tuple[pd.DataFrame, list]:
    """Tracks and bands for three events that scatter jump_num_max times each."""
    config = Config(config_path=runfiles.Create().Rlocation("_main/config_files/example.yaml"))
    config.physics.events_to_simulate = 3
    config.physics.betas_to_simulate = 200
    config.physics.energy = {"distribution": "fixed", "value": 18.6e3}
    config.physics.rho = {"distribution": "fixed", "value": 0.001}
    config.physics.z = {"distribution": "fixed", "value": 0.0}
    config.physics.min_theta, config.physics.max_theta = 89.5, 90.0
    config.trackbuilder.update(
        jump_num_max=jump_num_max,
        frac_elastic=0.0,
        track_length={"distribution": "fixed", "value": 1e-6},
        scattering_angle={"distribution": "fixed", "value": 0.0},
        energy_loss={"distribution": "fixed", "value": 0.0},
        start_time={"distribution": "fixed", "value": 0.0},
    )
    return TrackBuilder(config).run(EventBuilder(config).run())


class TrackBuilderTest(unittest.TestCase):
    def test_n_jumps_give_n_plus_one_tracks_per_event(self) -> None:
        for jump_num_max in (0, 1, 4):
            with self.subTest(jump_num_max=jump_num_max):
                tracks, bands = _build(jump_num_max)
                self.assertEqual(set(tracks.groupby("event_num").size()), {jump_num_max + 1})
                self.assertEqual([len(event_bands) for event_bands in bands], [jump_num_max + 1] * len(bands))

    def test_each_track_starts_where_the_previous_one_ended(self) -> None:
        tracks, _ = _build(4)
        for _, event_tracks in tracks.groupby("event_num"):
            self.assertEqual(list(event_tracks["track_num"]), [0, 1, 2, 3, 4])
            self.assertEqual(list(event_tracks["start_time"][1:]), list(event_tracks["end_time"][:-1]))

    def test_every_track_has_positive_power(self) -> None:
        tracks, _ = _build(4)
        self.assertTrue((tracks["track_power"] > 0).all())


if __name__ == "__main__":
    unittest.main()
