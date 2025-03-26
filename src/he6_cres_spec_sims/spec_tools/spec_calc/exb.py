import numpy as np

class ExB:
    """
        Helper functions for easily computing quantities related to the ExB/ vaunix
        e.g.) given vaunix voltage ON/OFF times & phase, given a time, when is the next ExB pulse?
    """
    def __init__(self, voltage_off_time=np.inf, voltage_on_time=0, voltage_fractional_offset=0):

        self.tVoltageOFF = voltage_off_time
        self.tVoltageON = voltage_on_time
        self.tVoltagePeriod = self.tVoltageON + self.tVoltageOFF
        self.tShift =  voltage_fractional_offset * self.tVoltagePeriod

    def time_in_trap_acq(self, t):
        # given loaded vaunix parameters and a time, returns the time in the trap acquisition
        # t_trap_in_acq defined as 0 at instant the ExB voltage turns off
        return (t - self.tVoltageON - self.tShift) % self.tVoltagePeriod

    def trap_cycle_index(self, t):
        # given loaded vaunix parameters and a time, returns the trap acquisition number
        # same as time_in_trap_acq, but returns quotient instead of remainder
        # return integer values >= -1
        return (t - self.tVoltageON - self.tShift) // self.tVoltagePeriod

    def next_empty(self, t):
        # given loaded vaunix parameters and a time, returns the onset of the next ExB pulse
        t_in_trap_acq = self.time_in_trap_acq(t)
        return  t - t_in_trap_acq + self.tVoltageOFF + (t_in_trap_acq > self.tVoltageOFF) * self.tVoltagePeriod

    def vaunix_time_series(self, t, fSignal):
        # given loaded vaunix parameters and a time series, return the vaunix signal (in data)
        # assumes envelope x high-frequency pulse. We don't know how phases are correlated across pulses
        # modulus creates periodic vaunix pulse. "%" operator does work for floats

        #TODO: it is unclear how the high-frequency phase of the vaunix is correlated across pulses.
        # Should/ could set to random [0,2 pi]. Unobservable without time-domain or complex data. Punt for now

        frac_time = ((t - self.tShift) / self.tVoltagePeriod)
        mask = ((frac_time - np.floor(frac_time)) < (self.tVoltageON / self.tVoltagePeriod))
        #mask = ((t - self.tShift) % self.tVoltagePeriod)<self.tVoltageON
        result = np.zeros_like(t)
        result[mask] = np.sin(2 * np.pi * fSignal * t[mask])
        return result
