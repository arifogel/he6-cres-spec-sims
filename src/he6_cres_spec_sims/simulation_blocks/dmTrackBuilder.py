import numpy as np

class DMTrackBuilder:
    """ Downmixes start_freq and end_freq of simulated tracks to observed frequency band out of DAQ
    """

    def __init__(self, config):
        self.config = config

    def run(self, tracks_df, events):
        print("~~~~~~~~~~~~DMTrackBuilder Block~~~~~~~~~~~~~~\n")
        mixer_freq = self.config.downmixer.mixer_freq
        print( "Downmixing the cyclotron frequency with a {} GHz signal".format( np.around(mixer_freq * 1e-9, 4)))
        downmixed_tracks_df = tracks_df.copy()
        downmixed_tracks_df["start_freq"] -= mixer_freq
        downmixed_tracks_df["end_freq"] -= mixer_freq

        for event in events:
            for band in event:
                band.shift_frequency(-mixer_freq)

        return downmixed_tracks_df
