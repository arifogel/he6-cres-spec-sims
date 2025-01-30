import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc

class Band:
    '''
    A class for different types of bands, the most common being a linear band. All bands have a start
    and end frequency and time and power, but  have different shapes and integrals.
    '''

    def __init__(self, start_time, start_freq, event, track, band, _power=None, _end_freq=None, _end_time=None,
                  _band_type=None):
        self.start_freq = start_freq
        self.start_time = start_time
        self.event = event
        self.track = track
        self.band = band
        self.power = _power
        self.end_freq = _end_freq
        self.end_time = _end_time
        self.band_type = _band_type

    def set_power(self, power):
        self.power = power
        return self

    def set_band(self, band):
        self.band = band
        return self

    def shift_frequency(self, shift):
        self.start_freq += shift
        self.end_freq += shift
        return self

    def copy(self):
        return Band(self.start_freq, self.start_time,self.event,self.track,self.band,self.power,self.end_freq,
                          self.end_time, self.band_type)

    def __repr__(self):
        return f"{self.band_type} Track"

    def __str__(self):
        return f"{self.band_type} Track \n Event: {self.event} \n Track: {self.track} \n Band: {self.band}"


class LinearBand(Band):
     def __init__(self, start_time, start_freq, total_power, event, track, band, max_time, max_freq, field):
        super().__init__(start_time, start_freq, event, track, band)
        self.track_type = "Linear"

        start_energy = sc.freq_to_energy(start_freq, field)
        self.slope = sc.df_dt( start_energy, field, total_power)

        if self.slope*max_time < max_freq:
            self.end_freq = self.slope*max_time + start_freq
            self.end_time = max_time
        else:
            self.end_freq = max_freq
            self.end_time = (max_freq-start_freq)/self.slope

        #print(f"Slope: {self.slope} \n t0: {start_time}, tf: {self.end_time} \n f0: {start_freq}, ff: {self.end_freq}")
