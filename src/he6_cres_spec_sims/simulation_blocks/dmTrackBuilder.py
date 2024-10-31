import numpy as np

class DMTrackBuilder:
    """ Downmixes freq_start and freq_stop of simulated tracks to observed frequency band out of DAQ
    """

    def __init__(self, config):
        self.config = config

    def run(self, tracks_df, events):
        print("~~~~~~~~~~~~DMTrackBuilder Block~~~~~~~~~~~~~~\n")
        mixer_freq = self.config.downmixer.mixer_freq
        print( "Downmixing the cyclotron frequency with a {} GHz signal".format( np.around(mixer_freq * 1e-9, 4)))
        downmixed_tracks_df = tracks_df.copy()
        downmixed_tracks_df["freq_start"] = ( downmixed_tracks_df["freq_start"] - mixer_freq)
        downmixed_tracks_df["freq_stop"] = downmixed_tracks_df["freq_stop"] - mixer_freq

        #TODO need to see how slow this is because I do not like having this quadruple nested for loop...
        for event in events:
            for track in event:
                for band in track:
                    for segment in track[band]:
                        segment.shift_frequency(-mixer_freq)


        return downmixed_tracks_df
