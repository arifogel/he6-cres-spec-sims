import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc
import pandas as pd
from he6_cres_spec_sims.simulation_blocks.trackBuilder import * 

class BandBuilder:
    """ Constructs list of sidebands and powers for trapped "segments", between scatters
    """

    def __init__(self, config):

        self.config = config

    def run(self, tracks_df, tracks):

        print("~~~~~~~~~~~~BandBuilder Block~~~~~~~~~~~~~~\n")
        sideband_num = self.config.bandbuilder.sideband_num
        magnetic_modulation = self.config.bandbuilder.magnetic_modulation
        harmonic_sidebands = self.config.bandbuilder.harmonic_sidebands

        frac_total_track_power_cut = self.config.bandbuilder.frac_total_track_power_cut
        total_band_num = sideband_num * 2 + 1
        band_list = []

        for tracks_index, row in tracks_df.iterrows():

            if harmonic_sidebands:
                sideband_amplitudes = sc.sideband_calc(
                    row["energy"],
                    row["rho_center"],
                    row["avg_cycl_freq"],
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
                    row["avg_cycl_freq"],
                    row["axial_freq"],
                    row["zmax"],
                    self.config.trap_profile,
                    magnetic_modulation=magnetic_modulation,
                    num_sidebands=sideband_num,
                )[0]

            for i, band_num in enumerate(range(-sideband_num, sideband_num + 1)):

                if sideband_amplitudes[i][1] < frac_total_track_power_cut:
                    continue
                else:
                    # copy track in order to fill in band specific values
                    row_copy = row.copy()

                    # fill in new avg_cycl_freq, band_power, band_num
                    # TODO: properly determine band power stop.
                    row_copy["avg_cycl_freq"] = sideband_amplitudes[i][0]
                    # Note that the sideband amplitudes need to be squared to give power.
                    row_copy["band_power_start"] = sideband_amplitudes[i][1] ** 2 * row.track_power
                    row_copy["band_power_stop"] = row_copy["band_power_start"]
                    row_copy["band_num"] = band_num

                    # append to band_list, as it's better to grow a list than a df
                    band_list.append(row_copy.tolist())

                    if band_num==0:
                        for track in tracks[int(row_copy["event_num"])][int(row_copy["track_num"])][0]:
                            track.set_power(row_copy["band_power_start"])

                    else:
                        track_subset = tracks[int(row_copy["event_num"])][int(row_copy["track_num"])][0]
                        freq_shift = row_copy["avg_cycl_freq"] - row["avg_cycl_freq"]
                        new_tracks = [track.copy().shift_frequency(freq_shift).set_band(band_num) for track in track_subset]
                        tracks[int(row_copy["event_num"])][int(row_copy["track_num"])][band_num] = new_tracks

        bands_df = pd.DataFrame(band_list, columns=tracks_df.columns)

        return bands_df
