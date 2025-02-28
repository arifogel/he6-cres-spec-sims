from he6_cres_spec_sims.constants import *

class Band:
    '''
    Base class for different types of bands (class wrapper for inter-scatter CRES signals, distinguished from the underlying betas)
    The linear band is the simplest band model
    All bands have a start and end frequency and time and power, but have different shapes and integrals.
    Child classes need to implement f(t), t(f), Phi(t)
    NOTE: in the base class constructor, start_time, start_freq, end_time ARE modified to fit in the BW (end_freq is computed)
    '''

    def __init__(self, start_time, start_freq, end_time, min_freq, max_freq, event, track, band):
        self.start_freq = start_freq
        self.start_time = start_time
        self.end_time = end_time
        self.end_freq = self.f(end_time)

        self.min_freq = min_freq
        self.max_freq = max_freq

        #IDs for event, track, band (i.e. beta, scatter-free time interval, sideband/mainband order f_c + k * f_a)
        self.event = event
        self.track = track
        self.band = band

        self.power = None
        self.end_freq = None
        self.band_type = None

        # Returns true if no part of the band is within the bandwidth (so no need to write or store this band)
        # Necessary to create band before we know if it is in bandwidth, for stuff "from below". May be born outside BW, enter BW later
        self.outside_BW = ((self.start_freq > self.max_freq) or (self.end_freq < self.min_freq))

    def shrink_to_BW(self):
        # Modify start_time, start_freq, end_time for entering/exiting bandwidth. Necessary to prevent aliasing
        # Note the two separate if statements: both are necessary (not an elif) for an event passing all the way through the BW
        # If you are copying a mainband to create sidebands, create sidebands THEN apply this function
        # otherwise, a sideband will have a modified start/end time based on when mainband enters/leaves BW (incorrect)
        if self.end_freq > self.max_freq:
            self.end_freq = self.max_freq
            self.end_time = self.t(self.max_freq)
        if self.start_freq < self.min_freq:
            self.start_freq = self.min_freq
            self.start_time = self.t(self.min_freq)

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

    def f(self, t):
        #frequency of time for band, for time series signal generation
        raise Exception("function f(t) not implemented by child class. Do so!")
        return None

    def t(self, f):
        #time of frequency for band, for alias prevention.
        #It is on the child class implementation for f(t), t(f) to be consistent!
        raise Exception("function t(f) not implemented by child class. Do so!")
        return None

    def Phi(self, t):
        #Cyclotron phase vs time for band. Integral \int ω(t') from 0 to t. Recommended to have Phi(start_time) = 0.
        #It is on the class implementation for f(t), Phi(t) to be consistent!
        raise Exception("function Phi(t) not implemented by child class. Do so!")
        return None

    def __repr__(self):
        return f"{self.band_type} Track"

    def __str__(self):
        return f"{self.band_type} Track \n Event: {self.event} \n Track: {self.track} \n Band: {self.band}"


class LinearBand(Band):
     def __init__(self, start_time, start_freq, end_time, min_freq, max_freq, event, track, band, slope):
        super().__init__(start_time, start_freq, end_time, min_freq, max_freq, event, track, band)
        self.track_type = "linear"

        self.slope = slope

    def f(self, t):
        return self.start_freq + self.slope * (t - self.start_time)

    def t(self, f):
        return self.start_time + (f - self.start_freq) / self.slope

    def Phi(self, t):
        #Phi(start_time) = 0.
        return 2*PI*(self.start_freq * (t-self.start_time)  + self.slope / 2. * (t - self.start_time)**2)

