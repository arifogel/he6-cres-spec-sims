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
        total_band_num = sideband_num * 2 + 1
        band_list = []

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

            # copy track in order to fill in band specific values
            row_copy = row.copy()
            sidebands= []
            for i, band_num in enumerate(range(-sideband_num, sideband_num + 1)):
                if sideband_amplitudes[i][1] < frac_total_track_power_cut:
                    continue
                else:
                    # fill in new avg_cycl_freq, band_power, band_num
                    # TODO: properly determine band power stop.
                    row_copy["start_freq"] = sideband_amplitudes[i][0]
                    # Note that the sideband amplitudes need to be squared to give power.
                    row_copy["start_band_power"] = sideband_amplitudes[i][1] ** 2 * row.track_power
                    row_copy["end_band_power"] = row_copy["start_band_power"]
                    row_copy["band_num"] = band_num

                    # append to band_list, as it's better to grow a list than a df
                    band_list.append(row_copy.tolist())

                    if band_num==0:
                        bands[int(row_copy["event_num"])][int(row_copy["track_num"])].set_power(row_copy["start_band_power"])

                    else:
                        band_copy = bands[int(row_copy["event_num"])][int(row_copy["track_num"])].copy()
                        freq_shift = row_copy["start_freq"] - row["start_freq"]
                        new_track = band_copy.shift_frequency(freq_shift).set_band(band_num)
                        sidebands.append(new_track)
            #print(sidebands)
            bands[int(row_copy["event_num"])] += sidebands

        bands_df = pd.DataFrame(band_list, columns=tracks_df.columns)

        return bands_df, bands
