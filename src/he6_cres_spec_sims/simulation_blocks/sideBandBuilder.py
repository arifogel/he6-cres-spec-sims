import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc
import pandas as pd
from he6_cres_spec_sims.simulation_blocks.trackBuilder import *

class SideBandBuilder:
    """ Constructs list of sidebands and powers from main bands made in trackbuilder
    """

    def __init__(self, config):

        self.config = config

    def run(self, tracks_df, bands):

        print("~~~~~~~~~~~~SideBandBuilder Block~~~~~~~~~~~~~~\n")
        sideband_num = self.config.sidebandbuilder.sideband_num
        magnetic_modulation = self.config.sidebandbuilder.magnetic_modulation
        harmonic_sidebands = self.config.sidebandbuilder.harmonic_sidebands

        frac_total_track_power_cut = self.config.sidebandbuilder.frac_total_track_power_cut

        for tracks_index, row in tracks_df.iterrows():
            if harmonic_sidebands:
                sideband_amplitudes = sc.sideband_calc(
                    row["energy"],
                    row["rho_center"],
                    row["start_freq"],
                    row["axial_freq"],
                    row["zmax"],
                    self.config.trap_profile,
                    magnetic_modulation=magnetic_modulation,
                    num_sidebands=sideband_num,
                )[0]
            else:
                sideband_amplitudes = sc.anharmonic_sideband_calc(
                    row["energy"],
                    row["center_theta"],
                    row["rho_center"],
                    row["start_freq"],
                    row["axial_freq"],
                    row["zmax"],
                    self.config.trap_profile,
                    magnetic_modulation=magnetic_modulation,
                    num_sidebands=sideband_num,
                )[0]

            sidebands = []

            for i, band_num in enumerate(range(-sideband_num, sideband_num + 1)):
                if sideband_amplitudes[i][1] > frac_total_track_power_cut:
                    # fill in new avg_cycl_freq, band_power, band_num
                    start_freq = sideband_amplitudes[i][0]
                    # Note that the sideband amplitudes need to be squared to give power.
                    band_power = sideband_amplitudes[i][1] ** 2 * row.track_power

                    freq_shift = start_freq - row["start_freq"]
                    new_track = bands[int(row["event_num"])][int(row["track_num"])].copy()
                    new_track.shift_frequency(freq_shift)
                    new_track.set_band(band_num)
                    new_track.set_power(band_power)
                    sidebands.append(new_track)

            bands[int(row["event_num"])] = sidebands

        return bands
