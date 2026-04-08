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

    def get_voltage_on_slices(self, tMin, tMax, N):
        #given t = np.linspace(tMin, tMax, N), return the index range(s) in which the vaunix is on
        # we don't actually want to construct the linspace array, as it is very slow to repeatedly allocate/deallocate this
        dt = (tMax - tMin) / (N-1)
        kRange = np.array([(tMin - self.tShift) / self.tVoltagePeriod, (tMax-self.tShift) / self.tVoltagePeriod])
        kRange = kRange.astype(int)
        slices = []
        for k in range(kRange[0], kRange[1]+1):
            nRising = max(int((self.tShift + k*self.tVoltagePeriod - tMin) / dt), 0)
            nFalling = min(int((self.tShift + k*self.tVoltagePeriod + self.tVoltageON - tMin) / dt), N-1)
            if (nRising < N) and (nFalling >= 0):
                slices.append(slice(nRising, nFalling))

        return slices
