"""Tests for SideBandBuilder, on one hand-built track and the example config's trap."""

from __future__ import annotations

import unittest

import pandas as pd
from python.runfiles import runfiles

from he6_cres_spec_sims.simulation_blocks.Band import LinearBand
from he6_cres_spec_sims.simulation_blocks.config import Config
from he6_cres_spec_sims.simulation_blocks.sideBandBuilder import SideBandBuilder

_SIDEBAND_NUM = 10
_TRACK_POWER = 1e-15
_START_FREQ = 2.0e10


def _sidebands(*, harmonic: bool, magnetic_modulation: bool, power_cut: float = 0.0) -> list[LinearBand]:
    """The sidebands SideBandBuilder makes for one track."""
    config = Config(config_path=runfiles.Create().Rlocation("_main/config_files/example.yaml"))
    config.sidebandbuilder.update(
        sideband_num=_SIDEBAND_NUM,
        frac_total_track_power_cut=power_cut,
        harmonic_sidebands=harmonic,
        magnetic_modulation=magnetic_modulation,
    )
    tracks = pd.DataFrame(
        [
            {
                "energy": 18.6e3,
                "center_theta": 89.0,
                "rho_center": 0.0,
                "start_freq": _START_FREQ,
                "axial_freq": 2.787875e7,
                "zmax": 0.008227,
                "track_power": _TRACK_POWER,
                "event_num": 0,
                "track_num": 0,
            }
        ]
    )
    main_band = LinearBand(0, 0.0, _START_FREQ, 1e-3, 1.7e10, 2.24e10, 0, 0, 0, 1e9)
    return SideBandBuilder(config).run(tracks, [[main_band]])[0]


class SideBandBuilderTest(unittest.TestCase):
    def test_every_sideband_is_kept_when_the_power_cut_is_zero(self) -> None:
        for harmonic in (False, True):
            for magnetic_modulation in (False, True):
                with self.subTest(harmonic=harmonic, magnetic_modulation=magnetic_modulation):
                    bands = _sidebands(harmonic=harmonic, magnetic_modulation=magnetic_modulation)
                    self.assertEqual([b.band for b in bands], list(range(-_SIDEBAND_NUM, _SIDEBAND_NUM + 1)))

    def test_a_higher_power_cut_keeps_fewer_sidebands(self) -> None:
        uncut = _sidebands(harmonic=True, magnetic_modulation=False, power_cut=0.0)
        cut = _sidebands(harmonic=True, magnetic_modulation=False, power_cut=0.3)
        self.assertLess(len(cut), len(uncut))

    def test_harmonic_sidebands_without_magnetic_modulation_conserve_the_track_power(self) -> None:
        bands = _sidebands(harmonic=True, magnetic_modulation=False)
        self.assertAlmostEqual(sum(b.power for b in bands) / _TRACK_POWER, 1.0, places=6)


if __name__ == "__main__":
    unittest.main()
