from scipy import interpolate
import pandas as pd
import numpy as np
from time import process_time
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

import he6_cres_spec_sims.spec_tools.spec_calc.exb as exb
from he6_cres_spec_sims.spec_tools.spec_calc.spec_calc import waveguide_beta
from he6_cres_spec_sims.constants import *

class DAQ:
    """  If called, this module  passes through list of downmixed bands through the DAQ, producing fake .spec(k) files
         These can be passed through Katydid, identically to data
         Converts bands to time-domain signals, s(t). FFTs give S(f), for each slice
         Data = |S(f) + N(f)|**2, converted to uint8, and written to .spec(k) files
    """

    def __init__(self, config):

        self.config = config

        # DAQ parameters derived from the config parameters.
        self.delta_f = config.daq.freq_bw / config.daq.freq_bins # frequency bin size
        self.delta_t = 1 / self.delta_f # time bin size (before averaging/ tossing)
        self.slice_time = self.delta_t * self.config.daq.roach_avg # time bin size (after averaging/ tossing)
        self.pts_per_fft = config.daq.freq_bins * 2
        self.freq_axis = np.linspace( 0, self.config.daq.freq_bw, self.config.daq.freq_bins)

        self.antenna_z = 50  # Ohms

        self.slices_in_roach = int( config.daq.acq_length/ self.delta_t) # num slices in "data" (before averaging/tossing)
        self.slices_in_spec = int( config.daq.acq_length/ self.delta_t / self.config.daq.roach_avg) # num slices in file (after averaging/tossing)

        self.n_acquisitions = self.config.daq.n_acquisitions #number of seconds of data to write
        self.n_channels = self.config.daq.n_channels #1 or 2, for lower/ upper halves of spectrogram
        self.bins = [slice(0,4096), slice(4096,8192)] #which frequency bins to write for each channel

        # This block size is used to create chunks of spec file that don't overwhelm the ram.
        self.slice_block = int(250 * 32768 / config.daq.freq_bins) * self.config.daq.roach_avg

        # Read in .spec files for noise-only reference
        self.noise_mean = np.ones(config.daq.freq_bins, dtype=float)
        try:
            for n in range(self.n_channels):
                #Should we do more than 10k slices read in? Perhaps...
                self.noise_mean[self.bins[n]] = self.spec_to_array(self.config.daq.noise_paths[n]).mean(axis=0)
        except Exception as e:
            print("Noise loading failed!")
            print(str(e))

        self.noise_tau = np.nan_to_num(1./np.log(1 + 1./self.noise_mean), 0)

        #amplitude gain g_overall(f) experienced by both signal and noise. Class object is interpolation function g(f)
        #If frequency outside of bandwidth, automatically returns g(f) = 0 (aka, alias prevention)
        self.gain_overall_array = self.estimate_gain()

        #Save signal gain as arrays to reference for later (much faster than generating on the fly for each signal)
        self.signal_gain_array = self.signal_gains(self.config.sidebandbuilder.sideband_num, self.freq_axis)

        # Fast estimation of zero-suppression thresholds
        # We call it for both spec and speck, so that the rng for fake bands are the same
        self.thresholds = self.set_thresholds()

        self.ExB = exb.ExB(self.config.trackbuilder.voltage_off_time_ms/1000., self.config.trackbuilder.voltage_on_time_ms/1000., self.config.trackbuilder.voltage_fractional_offset)

    def estimate_gain(self):
        G = self.noise_tau * 1. #multiply for copy by reference

        #Get list of all peaks. Decrease prominence for more sensitivity to small peaks, though prone to picking up erroneous peaks
        noise_peaks, properties = find_peaks(self.noise_tau, width=3, prominence=0.5)
        #Pick out only the noise resonances (narrow) instead of gain/SNR oscillations (wide)
        noise_peaks = noise_peaks[properties["widths"] < 40]

        bins = np.arange(self.config.daq.freq_bins)
        fit_width = 13
        for peak in noise_peaks:
            #region Lorentz + linear fit is computed over
            binRange = [peak - fit_width, peak + fit_width]
            mask = ( binRange[0] < bins)*(bins< binRange[1])
            #Fitting model (linear background + lorentzian peak) (linear part defined with respect to bin_min)
            fLinearLorentzian = lambda f, a, b, c, f0, HWHM: a + b*(f-binRange[0]) + c / (1 + ((f-f0)/HWHM)**2)
            x = bins[mask]
            y = G[mask]
            popt, pcov = curve_fit(fLinearLorentzian, x, y, bounds=([0, -1,0,min(x), 1],[max(y), 1, 20., max(x), 20]))
            popt[0:2] = 0 #Set linear parameters to 0, leaving only Lorentzian component. Subtract from tau to get gain w/o resonances
            G -= fLinearLorentzian(bins,*popt)

        #prevent negative gains after subtracting out the Lorentzians
        G = np.clip(G,a_min=0,a_max=None)
        #Uncomment and plot if debugging? Could become permanent feature. Should agree except at resonances, where G looks like resonances are subtracted out
        #np.savetxt("gains.txt", G)
        #np.savetxt("tau.txt", self.noise_tau)

        return np.sqrt(G/2.)

    def signal_gains(self, max_sideband_order, f, r=[0.10,0.10], L = [0.15,0.92]):
        gSignal = np.ones(shape=(max_sideband_order+1, f.size)) + 0j
        #perhaps there is a cleaner/clearer/more clever way to do this. Sum over reflective surfaces, compute for each sideband
        #(-1)^s only valid exactly for harmonic traps
        for i in range(len(r)):
            vTmp = r[i] * np.exp(2 * 1j * waveguide_beta(2*np.pi*(f + self.config.downmixer.mixer_freq)) * L[i])
            for s in range(max_sideband_order + 1):
                gSignal[s,:] += (-1)**s * vTmp

        return np.abs(gSignal)

    def run(self, bands):
        """
        This function is responsible for building out the spec files and calling the below methods.
        """
        # Flatten into a 1D NumPy array
        self.bands = np.hstack(bands)

        # Define a random phase for each band. Need to be associated per track (lasting multiple chunks)
        # TODO: This is technically (actually) incorrect, there is an overall random phase that arises from
        # initial particle position. Different bands in the same event have correlated phases depending on z0
        # Should be done earlier, probably
        #self.bands["phi_0"] = self.config.dist_interface.rng.uniform(0,2 * PI, size=len(self.bands))

        self.create_results_dir()
        self.spec_file_paths = self.build_file_paths(self.n_acquisitions, self.n_channels, self.spec_files_dir)
        self.write_empty_files(self.spec_file_paths)

        spec_array = np.zeros(shape=(self.slice_block, self.config.daq.freq_bins))
        initial_packet = 0

        for acq in range(self.n_acquisitions):
            print( f"Building spec acquistion {acq}. {self.config.daq.acq_length} s, {self.slices_in_spec} slices.")
            build_file_start = process_time()
            # Iterate by the slice_block until you hit the end of the spec file.
            for start_slice in np.arange(0, self.slices_in_roach, self.slice_block):
                stop_slice = min(start_slice + self.slice_block, self.slices_in_roach)

                num_slices = stop_slice - start_slice
                requant_gain_scaling = 2**self.config.daq.requant_gain

                #dimensions: slices x FFT bins
                spec_array = self.get_signal_array( acq, start_slice, stop_slice)
                spec_array *= np.sqrt( requant_gain_scaling)

                #shape[1] is the number of slices, though by doing it like this, we handle automatically if
                # roach_inverted_flag=True (so this is number of slices either before or after summing/tossing)
                spec_array += self.get_noise_array(spec_array.shape[0])

                # Computer Fourier power (magnitude_squared)
                spec_array = np.abs(spec_array)**2

                if not self.config.daq.roach_inverted_flag:
                    spec_array = self.roach_slice_sum(spec_array)

                spec_array = np.clip(spec_array, a_min=0, a_max=255)

                # Write chunk to spec file.
                for channel in range(self.n_channels):
                    # self.bins[channel] tells to write frequency bins [0-4095], [4096,8191]
                    # in _0.spec(k) and _1.spec(k) respectively. First ":" indicates write all time slices
                    if self.config.daq.spec_suffix == "spec":
                        self.write_to_spec(spec_array[:,self.bins[channel]], self.spec_file_paths[acq][channel], initial_packet)
                    elif self.config.daq.spec_suffix == "speck":
                        self.write_to_speck(spec_array[:,self.bins[channel]], self.spec_file_paths[acq][channel], initial_packet, channel)
                    else:
                        raise ValueError('Invalid spec_suffix: spec || speck')

                initial_packet += spec_array.shape[0]
                initial_packet = initial_packet % 2**20

            build_file_stop = process_time()
            print( f"Time to build acq {acq}: {build_file_stop- build_file_start:.3f} s \n")

        print("Done building {} files. ".format(self.config.daq.spec_suffix))

    def get_signal_time_series(self, acq, start_slice, stop_slice):
        """
        Build a time-domain array of signal (Dimensions = N_FFT Bins x num_slices)
        Later, this will be converted to the frequency domain S(f) via FFT, with the same dimensions
        """
        print(f"acq = {acq}, slices = [{start_slice}:{stop_slice}]")
        slice_start_time = start_slice * self.delta_t
        slice_stop_time = stop_slice * self.delta_t
        num_slices = stop_slice - start_slice

        t = np.linspace(slice_start_time, slice_stop_time, self.pts_per_fft * num_slices )
        dt = t[1] - t[0]
        signal_time_series = np.zeros(shape=self.pts_per_fft * num_slices)

        # shape of signal_alive_condition: num_bands
        signal_alive_condition = np.where([(
            (b.outside_BW == False)
            & (b.acquisition == acq)
            & (b.start_time <= slice_stop_time)
            & (b.end_time >= slice_start_time))
            for b in self.bands])[0]

        eligible_bands = self.bands[signal_alive_condition]

        # Sum all signals in bandwidth to get total (CRES) time-series, to be FFT'ed
        # The factor of 2 is needed because the instantaneous frequency is the derivative of the phase
        # The band_phase is a random phase assigned to each band.
        time_to_index = lambda tTime: int( (tTime - t[0]) / dt)

        for band in eligible_bands:
            # Slice object - selects the time indices in which the band is active
            band_mask = slice(max(time_to_index(band.start_time), 0), time_to_index(band.end_time), 1)

            # TODO: Put back in time-dependence of amplitudes. Want to add in frequency-dependence too
            inds = (band.f(t[band_mask]) / self.delta_f).astype(int)
            #The abs is so you pick out the ORDER of the sidebands for the signal gain profile. Negative indices don't work here
            voltage = np.sqrt(band.power * self.antenna_z) * self.signal_gain_array[abs(band.band)][inds]
            signal_time_series[band_mask] += voltage * np.sin( band.Phi(t[band_mask]))

        return signal_time_series.reshape((num_slices, self.pts_per_fft)).transpose()

    def get_vaunix_time_series(self, acq, start_slice, stop_slice):
        """
        Build a time-domain array of vaunix signal (Dimensions = N_FFT Bins x num_slices)
        Later, this will be converted to the frequency domain S(f) via FFT, with the same dimensions
        """
        slice_start_time = start_slice * self.delta_t
        slice_stop_time = stop_slice * self.delta_t
        num_slices = stop_slice - start_slice

        t = np.linspace(slice_start_time, slice_stop_time, self.pts_per_fft * num_slices )
        fVaunix = self.config.daq.vaunix_bin * self.config.daq.freq_bw / self.config.daq.freq_bins

        # vaunix power scaling given the axolotl controls (power in dB)
        # Nick found this by setting voltage_off_time = 0, producing a spectrogram.
        # From https://drive.google.com/file/d/197czYZ2x9wSeNMNPTpIlJNUG2gTmIewV/view?usp=sharing (slide 11)
        # we should get an average power in the spectrogram of 98.3 at input of -1 dB. Scale reference power accordingly
        reference_power = 4.747e-6
        voltage = np.sqrt(reference_power * self.antenna_z) * 10**(self.config.daq.vaunix_power_db / 20.)
        # modulus creates periodic vaunix pulse. "%" operator does work for floats
        vaunix_time_series =  voltage * self.ExB.vaunix_time_series(t, fVaunix)

        return vaunix_time_series.reshape((num_slices, self.pts_per_fft)).transpose()

    def get_signal_array(self, acq, start_slice, stop_slice):
        """
        Build a frequency-domain array of signal (Dimensions = some number of slices x FFT bins)
        The number of slices is usually slice_block. If roach_inverted_flag, it is usually slice_block / roach_avg, except near the end of the spec
        Given signal time-series s(t), convert to frequency domain S(f) via FFT
        Returns frequency-domain (with phase)
        """
        signal_time_series = self.get_signal_time_series(acq, start_slice, stop_slice)
        # LNA gain of 67dB (this should be a user parameter)
        signal_time_series *= np.sqrt(1e9)
        signal_time_series += self.get_vaunix_time_series(acq, start_slice, stop_slice)

        #shape of signal_time_series: (pts_per_fft, num_slices). Conduct a 1d FFT along axis = 0 (the time axis).
        #Avoid taking FFT's of time slices were are going to toss anyways!
        if self.config.daq.roach_inverted_flag:
            #Note: We want to keep every daq_roach_avg'th slice in the full data, not the data block
            slices_to_keep = np.arange(start_slice, stop_slice) % self.config.daq.roach_avg == 0
            signal_time_series = signal_time_series[:,slices_to_keep]

        Y_fft = np.fft.fft(signal_time_series, axis=0, norm="ortho")[:self.pts_per_fft // 2]

        #Multiply the Fourier transform of all signals by the overall amplifier gain profile
        #Y_fft shape (fourier bins x slices), gain shape (fourier bins x 1). Product multiplies each slice by gain
        Y_fft *= self.gain_overall_array[:,np.newaxis]

        return Y_fft.T


    def get_noise_array(self, num_slices):
        """
        Build a frequency-domain array of noise (Dimensions = N_FFT Bins x self.slice_block slices)
        For additive white Gaussian noise n, N(f) = FFT(n) is also Gaussian distributed
        Note that noise is scaled according to noise-only data (tau(f)).
        For the signal only, we will multiply by gain_tot(f), it is already included in the noise here
        Returns frequency-domain (with phase)
        """

        delta_f_12 = 2.4e9 / 2**13

        noise_power_scaling = self.delta_f / delta_f_12
        requant_gain_scaling = (2**self.config.daq.requant_gain) / (2**self.config.daq.noise_file_gain)
        noise_scaling = noise_power_scaling * requant_gain_scaling

        array_size = (num_slices, self.config.daq.freq_bins)

        # Additive white Gaussian noise has FFT which is complex Gaussian
        # Real time-series only imply symmetry between positive/negative frequencies. FFT still complex
        # Imaginary first gets array typing to complex128, avoid recast
        noise_array = 1j * self.config.dist_interface.rng.normal(size=array_size)
        noise_array += self.config.dist_interface.rng.normal(size=array_size)

        # Want to scale so that mean power agrees with config (based on Chi-Squared k=2 for unsummed bins)
        tau_noise = 1./np.log(1 + 1./ self.noise_mean)
        noise_array *= np.sqrt(tau_noise /  2.)
        # Scale by noise power
        noise_array *= noise_scaling

        return noise_array

    def roach_slice_sum(self, signal_array):
        #input array dimensions: nSlices x nFFT
        #WARNING: This breaks if nSlices in the block size is not divisible by roach_avg (i.e. we do incomplete sums over slices)
        num_slices = signal_array.shape[0]
        #need to do sum over slices that is divisible by roach_avg
        if num_slices %  self.config.daq.roach_avg:
            print("Num slices really should be divisible by roach_avg! Why is it not!? Trimming")
            num_slices_divisible = num_slices - (num_slices % self.config.daq.roach_avg)
            signal_array = signal_array[:num_slices_divisible,:] # Trim remainder rows

        #This command takes our nSlices x nFFT array, and groups the nSlice columns into groups of roach_avg and sums over them
        #the -1 is a filler telling numpy to automatically compute the number of "groups", which would be the slices post-summing
        return signal_array.reshape(-1, self.config.daq.roach_avg, signal_array.shape[1]).sum(axis=1)

    #################### File Writing Utilities ####################

    def build_file_paths(self, n_acqs, n_channels, files_dir):
        """
            Creates list of all filenames.
            Returns 2D array like: file_paths[file_id][channel_id]
            where file_id is like the acquisition or the second of data in the run
            and channel_id is 0,1 for 0-1200 MHz, or 1200-2400 MHz
        """
        file_paths = []
        for acq_id in range(n_acqs):
            acq_n_paths = []
            for channel_id in range(n_channels):
                file_path = files_dir / "{}_{}_{}.{}".format( self.config.daq.spec_prefix, acq_id, channel_id, self.config.daq.spec_suffix)
                acq_n_paths.append(file_path)
            file_paths.append(acq_n_paths)
        return file_paths

    def safe_mkdir(self, new_dir):
        # If new_dir doesn't exist, then create it.
        if not new_dir.is_dir():
            new_dir.mkdir()
            print("created directory : ", new_dir)

    def create_results_dir(self):
        # First make a results_dir with the same name as the config.
        config_name = self.config.config_path.stem
        parent_dir = self.config.config_path.parents[0]

        self.results_dir = parent_dir / config_name
        self.safe_mkdir(self.results_dir)

        self.spec_files_dir = self.results_dir / "spec_files"
        self.safe_mkdir(self.spec_files_dir)

    def write_empty_files(self, files):
        """
        Create empty files to be filled with data (to be appended later)
        files is list of paths of output spec(k) files organized
        [[acq0_0.spec(k), _1.spec(k)],[acq1_0.spec(k), _1.spec(k)]...]
        """
        for acq in files:
            for file_path in acq:
                open(file_path, "wb")

    def spec_to_array(self, spec_path, slices=10000, start_packet=0):
        """
            Stolen (though modified) from He6DAQ: Data_Quality_Control.py, with 2^15 bitcode deprecated (as in data)
            Needed to read in noise data (.spec files) to set the noise levels vs frequency
            Returns object with dimensions spec_array[slice][freq], unlike in Data_Quality_Control.py
        """
        BYTES_IN_HEADER = 32
        BYTES_IN_PAYLOAD = 4096
        BYTES_IN_PACKET = BYTES_IN_PAYLOAD + BYTES_IN_HEADER

        if slices == -1:
            spec_array = np.fromfile(spec_path, dtype="uint8", count=-1).reshape(
                (-1, BYTES_IN_PACKET)
            )[:, BYTES_IN_HEADER:]
        else:
            spec_array = np.fromfile(
                spec_path, dtype="uint8", count=BYTES_IN_PACKET * slices
            ).reshape((-1, BYTES_IN_PACKET))[:, BYTES_IN_HEADER:]

        return spec_array

    def packet_num_base_256(self, packet_num):
        # We want to write (packet number) to spec(k) file
        # packet_number (0 - 2^20-1)
        # Fit in 3 bytes via: 2^16 a[0] + 2^8 a[1] + a[2]

        aOnes = packet_num % 256
        intermediate = packet_num // 256 #2^8 a[0] + a[1]
        aTens = intermediate % 256
        aHunds = intermediate // 256
        return np.array([aHunds, aTens, aOnes])


    def write_to_spec(self, spec_array, spec_file_path, initial_packet):
        """
        Append to an existing spec file. This is necessary because the spec arrays get too large for 1s
        worth of data.
        """
        # Make spec file:
        slices_in_spec, freq_bins_in_spec = spec_array.shape

        #32 is the hardcoded header size, in bytes
        zero_hdrs = np.zeros((slices_in_spec, 32),dtype=int) #shape = slices x 32

        packets = np.arange(initial_packet,initial_packet + slices_in_spec)
        packets_base256 = self.packet_num_base_256(packets).transpose() #shape = slices x 3
        zero_hdrs[:,9:12] = packets_base256

        # Append mostly empty headers to the spec array.
        spec_array_hdrs = np.hstack((zero_hdrs, spec_array))

        data = spec_array_hdrs.flatten().astype("uint8")

        # Pass "ab" to append to a binary file
        with open(spec_file_path, "ab") as spec_file:
            # Write data to spec_file.
            data.tofile(spec_file)

        return None

    def set_thresholds(self):
        # For zero-suppression, need to set zero-suppression thresholds based on noise
        # Instead of generating 1s of noise for each frequency bin, use central limit theorem for mean
        # power in each bin. Sum of N Chi-squared (k=4,2) depending on if summed or not

        means = self.noise_mean
        # DOF = 2 when summing off (inverted_flag == True), 4 when summing on (inverted_flag == False)
        kDOF = 2
        if not self.config.daq.roach_inverted_flag:
            kDOF *= self.config.daq.roach_avg

        #number of slices used in average for zero-suppression thresholding. Tends to be slow to use 146,484. 10k good enough in practice
        nSlicesThresholding = 10000
        # CLT: sigma of sum = sigma(Chi-squared) / sqrt(N)
        sigma_thresholds = np.sqrt(2. / (kDOF * nSlicesThresholding))
        #Add some noise to the mean based on the number of samples (average over previous second)
        means *= self.config.dist_interface.rng.normal(1, sigma_thresholds, size=self.config.daq.freq_bins)

        # use tau from Non-Exponential noise doc: https://drive.google.com/file/d/10EGOZGXkmiXHXLeyHQ1qnNcc_HK8FPxj/view
        #thresholds = means
        thresholds = 1. / np.log(1. + 1./means)
        thresholds *= self.config.daq.threshold_factor
        thresholds = np.clip(thresholds, 1,None)
        return thresholds

    def add_high_power_point(self, frequency_bin):
        # We want to write (bin number, power) to zero-suppressed file
        # aIndex (0 - 4095) is a 12-bit number, does not fit in a byte.
        # Fit in 2 bytes via: index = 2^8 a[0] + a[1]

        aOnes = frequency_bin % 256
        aTens = (frequency_bin - aOnes) // 256

        return [aTens, aOnes]

    def write_to_speck(self, spec_array, speck_file_path, initial_packet, channel):
        """
        Append to an existing speck file. This is necessary because the raw spec arrays get too large for 1s
        worth of data.
        """
        slices_in_spec, freq_bins_in_spec = spec_array.shape

        # Append mostly empty packet header to data
        header = np.zeros(32)

        # Append empty (zero) footer. 3 zeros signals end of spectrogram slice
        footer = np.zeros(3)

        if self.config.daq.threshold_factor is None or self.config.daq.threshold_factor < 0:
                raise ValueError('Invalid DAQ::threshold_factor. Set to non-negative real value!')

        data = np.array([])

        #initial index (e.g. 0 or 4096 for channels 0,1) in thresholds to compare to
        jThreshold0 = channel * freq_bins_in_spec
        thresholds = self.thresholds[jThreshold0:jThreshold0+freq_bins_in_spec]

        # Pass "ab" to append to a binary file
        with open(speck_file_path, "ab") as speck_file:
            for s in range(slices_in_spec):
                header[9:12] = self.packet_num_base_256(initial_packet + s)
                data = np.append(data, header)
                # select indices of spectrogram [0-4096] above threshold. Loop is slow!
                indices = np.where(spec_array[s] >  thresholds)[0]
                for j in indices:
                    data = np.append(data, self.add_high_power_point(j))
                    data = np.append(data, spec_array[s][j])
                data = np.append(data, footer)

            data = data.flatten().astype("uint8")
            data.tofile(speck_file)

        #fractionHighPowerPoints =  (len(data) - (len(header) + len(footer))  * slices_in_spec)  / (slices_in_spec * freq_bins_in_spec)
        #print("Fraction passing 0-supp: ",fractionHighPowerPoints)

        return None
