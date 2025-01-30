from .eventBuilder import *
import he6_cres_spec_sims.spec_tools.spec_calc.power_calc as pc
from .Band import *

class TrackBuilder:
    """ Constructs a list of tracks (interrupted by scatters) making up the trapped event
    """

    def __init__(self, config):

        self.config = config
        self.eventbuilder = EventBuilder(config)

        # distribution of energy losses [eV]
        self.jump_distribution = config.dist_interface.get_distribution(self.config.trackbuilder.energy_loss)

        # distribution of track durations [s]
        self.track_length_distribution = config.dist_interface.get_distribution(self.config.trackbuilder.track_length)

        # distribution of scattering angles [degrees]
        self.scattering_angle_distribution = config.dist_interface.get_distribution(self.config.trackbuilder.scattering_angle)

        # distribution of start times [s]
        self.start_time_distribution = config.dist_interface.get_distribution(self.config.trackbuilder.start_time)

        self.verbosity = self.config.trackbuilder.verbose
        print(self.config.trackbuilder.verbose)

    def run(self, trapped_event_df):
        """
        Builds scattered tracks for each event.
        """
        print("~~~~~~~~~~~~TrackBuilder Block~~~~~~~~~~~~~~\n")
        # Empty list to be filled with tracks.
        tracks_list = []
        bands = []
        #create tracks for every event
        for event_index, event in trapped_event_df.iterrows():
            if event_index % 25 == 0:
                print("\nBuilding Event :", event_index)

            # Fill the event with computationally intensive properties.
            event = self.fill_in_properties(event)

            event["time_start"] = self.start_time_distribution.generate()
            event["freq_start"] = event["avg_cycl_freq"]

            # Assign track 0 of event with a scatter time.
            track_duration = self.track_length_distribution.generate()
            scatter_time = event["time_start"] + track_duration

            # Begin with trapped beta (track 0 of event).
            tracks = [event]
            is_trapped = True
            jump_num = 0
            track_num = 0

            #TODO this may need to be more nuanced to account for lower sidebands
            max_freq = self.config.physics.freq_acceptance_high

            #TODO maybe this should be a different variable in the config... or maybe just renamed
            trap_on_time = np.inf
            if self.config.trackbuilder.voltage_off_time_ms:
                trap_on_time = self.config.trackbuilder.voltage_off_time_ms * 1e-3

            ###XXX this is bah-roken
            #end_time = min(trap_on_time, scatter_time)
            end_time = event["time_start"] + 30.2e-3

            #list of band objects (to be added to bands list)
            event_main_bands = []

            '''
            Track building loop:
            Basic idea is we have 2 while loops: one for scattering, and one for frequency/time dependent bands
            The reason for this is the physical attributes of an event change after a scatter but dont for the other
            features handled by freq/time dependent features (this is not strictly true but close enough for us)
            '''
            # TODO: we need to add something to better handle ExB turning on. not sure if that should be treated like a
            # scatter or not.
            while is_trapped and jump_num<=self.config.trackbuilder.jump_num_max:
                if self.verbosity == True: print(f"Event {event_index}, Jump {jump_num}")
                t, freq, field = tracks[-1]["time_start"], tracks[-1]["freq_start"], tracks[-1]["b_avg"]
                track_radiated_power_tot = sc.power_larmor(field, freq)
                while (t < end_time) and (freq < max_freq):

                    band = self.create_band(t, freq, track_radiated_power_tot, end_time, max_freq, event_index,
                                                   jump_num, track_num, field)
                    t, freq = band.end_time, band.end_freq
                    event_main_bands.append(band)

                    tracks[-1]["freq_stop"] = freq
                    tracks[-1]["time_stop"] = t + 3e-3
                    tracks[-1]["energy_stop"] = sc.freq_to_energy(freq, tracks[-1]["b_avg"])
                    tracks[-1]["track_num"] = track_num

                    new_track = tracks[-1].copy()
                    new_track["time_start"] = t
                    new_track["freq_start"] = freq
                    tracks.append(new_track)
                    track_num += 1
                    tracks_list.append(tracks[-1].values.tolist())

                # break out of loop if this track reached end of trap on time
                if t >= trap_on_time: break

                new_track = self.scatter(tracks[-1])

                if self.eventbuilder.trap_condition(new_track) == True:
                    #TODO there is almost certainly a better way to pass/grab the new track info
                    new_track = next(self.fill_in_properties(new_track).iterrows())[1]
                    new_track["time_start"] = t
                    new_track["freq_start"] = new_track["avg_cycl_freq"]
                    tracks.append(new_track)
                    jump_num += 1
                    scatter_time = new_track["time_start"] + self.track_length_distribution.generate()
                    end_time = (trap_on_time, scatter_time) [scatter_time<trap_on_time]
                else:
                    is_trapped=False

            bands.append(event_main_bands)

        # TODO there may be a more elegant way to update the columns... but this works for now
        columns = np.append(trapped_event_df.columns.to_numpy(), ["time_start","freq_start","time_stop"])
        tracks_df = pd.DataFrame(tracks_list, columns=columns)


        return tracks_df, bands

    def scatter(self, event):
        """Creates Scattered track from initial event conditions.
        """

        center_x, center_y = event["center_x"], event["center_y"]
        rho_pos = event["initial_rho_pos"]
        phi_pos = event["initial_phi_pos"]
        zpos = 0
        center_theta = event["center_theta"]
        phi_dir = event["initial_phi_dir"]
        energy_stop = event["energy_stop"]
        event_num = event["event_num"]
        beta_num = event["beta_num"]

        # Jump Size
        jump_size_eV = self.jump_distribution.generate()

        # Delta Pitch Angle: Sampled from normal dist.
        scattering_angle = self.scattering_angle_distribution.generate()

        # Original beta direction vector in cartesian coordinates
        theta_dir = center_theta # tmp see TODO below
        v_vec_old = np.array([np.sin(theta_dir/RAD_TO_DEG) *np.cos(phi_dir/RAD_TO_DEG),
                    np.sin(theta_dir/RAD_TO_DEG)*np.sin(phi_dir/RAD_TO_DEG), np.cos(theta_dir/RAD_TO_DEG)])

        # This always produces a vector lying on cone with angle = scattering_angle with respect to initial velocity vector
        tmp_theta = theta_dir + scattering_angle
        tmp_v_vec = np.array([np.sin(tmp_theta/RAD_TO_DEG) *np.cos(phi_dir/RAD_TO_DEG),
                    np.sin(tmp_theta/RAD_TO_DEG)*np.sin(phi_dir/RAD_TO_DEG), np.cos(tmp_theta/RAD_TO_DEG)])

        # Using Rodrigues' Rotation Formula (rotate tmp_velocity_vec around velocity_vec_old)
        # Depending on dPhi, get random vector along that cone
        dPhi = 2*PI*self.config.dist_interface.rng.uniform()
        vNew = tmp_v_vec * np.cos(dPhi) + np.cross(v_vec_old,tmp_v_vec) * np.sin(dPhi) + v_vec_old * np.dot(v_vec_old,tmp_v_vec) * (1-np.cos(dPhi))
        vNew /= np.sqrt(np.dot(vNew, vNew)) # Probably unnecessary, don't want this to drift too much from floating point errors

        # Second, calculate new pitch angle and energy.
        theta_new = np.arccos( vNew[2] ) * RAD_TO_DEG
        phi_dir = np.arctan2(vNew[1],vNew[0]) * RAD_TO_DEG

        # TODO: We should 1) randomly choose z's to scatter at (PDF proportional to v_z^-1, which means saving z-motion, inverse transform sampling)
        # Then 2, convert local scattered pitch angle to the new center_theta
        # Since we scatter only at z=0, don't bother propagating scattering change in instantaneous theta_dir to center_theta
        center_theta = theta_new

        # Ensure that pitch angle is defined to be smaller than 90.
        if center_theta > 90:
            center_theta = 180 - center_theta

        # New energy:
        energy = energy_stop - jump_size_eV

        # New position and direction. Only center_theta is changing right now.
        beta_position, beta_direction = (
            [rho_pos, phi_pos, zpos],
            [center_theta, phi_dir],
        )

        # Third, construct a scattered, meaning potentially not-trapped, segment df
        return self.eventbuilder.construct_untrapped_track_df(beta_position, beta_direction, energy, event_num, beta_num)

    def fill_in_properties(self, incomplete_scattered_events_df):
        """ Assigns calculated properties (e.g. axial frequency, power, slope, etc.)
            to beta with given (E, theta, rho) in the magnetic field profile
        """

        df = incomplete_scattered_events_df.copy()
        trap_profile = self.config.trap_profile
        main_field = self.config.eventbuilder.main_field
        decay_cell_radius = self.config.eventbuilder.decay_cell_radius

        # Calculate all relevant track parameters. Order matters here.
        axial_freq = sc.axial_freq( df["energy"], df["center_theta"], df["rho_center"], trap_profile)

        # TODO: Make this more accurate as per discussion with RJ.
        b_avg = sc.b_avg( df["energy"], df["center_theta"], df["rho_center"], trap_profile, axial_freq)
        avg_cycl_freq = sc.energy_to_freq(df["energy"], b_avg)
        grad_b_freq = sc.grad_b_freq( df["energy"], df["center_theta"], df["rho_center"], trap_profile, axial_freq)
        zmax = sc.max_zpos( df["energy"], df["center_theta"], df["rho_center"], trap_profile)
        mod_index = sc.mod_index(avg_cycl_freq, zmax)


        track_radiated_power_te11 = (
            pc.power_calc(
                df["center_x"],
                df["center_y"],
                avg_cycl_freq,
                main_field,
                decay_cell_radius,
            )
            * 2
        )

        track_radiated_power_tot = sc.power_larmor(main_field, avg_cycl_freq)

        # slope = sc.df_dt( df["energy"], self.config.eventbuilder.main_field, track_radiated_power)

        energy_stop = ( df["energy"] - track_radiated_power_tot * df["track_length"] * J_TO_EV)

        # Replace negative energies if energy_stop is a float or pandas series
        if isinstance(energy_stop, pd.core.series.Series):
            energy_stop[energy_stop < 0]  = 1e-10
        elif energy_stop < 0:
            energy_stop = 1e-10

        freq_stop = sc.avg_cycl_freq( energy_stop, df["center_theta"], df["rho_center"], trap_profile)
        slope = (freq_stop - avg_cycl_freq) / df["track_length"]

        track_power = track_radiated_power_te11 / 2

        df["axial_freq"] = axial_freq
        df["avg_cycl_freq"] = avg_cycl_freq
        df["b_avg"] = b_avg
        df["grad_b_freq"] = grad_b_freq
        df["freq_stop"] = freq_stop
        df["energy_stop"] = energy_stop
        df["zmax"] = zmax
        df["mod_index"] = mod_index
        df["slope"] = slope
        df["track_power"] = track_power

        return df

    def create_band(self, time, freq, power, max_time, max_freq, event_num, track_num, band_num, b_avg):
        '''
        This function creates a band object for a given time and frequency. Currently only check that we are within
        track length.
        TODO add more band options
        '''

        # set different ranges of frequencies where different things can happen, the largest range is normal linear
        # tracks, but there can be cutoff regions, field slewing regions, etc where shape and slope changes

        # TODO currently code is creating events outside of physics frequency range... this buffer is a bandaid on
        # this problem which I believe is from the physics part of the code, bandaiding so I can continue debugging
        # this code...
        linear_range = [self.config.physics.freq_acceptance_low-0.1e9, self.config.physics.freq_acceptance_high]

        if linear_range[0] <= freq < linear_range[1]:
            band = LinearBand(time, freq, power, event_num, track_num, band_num, max_time, max_freq, b_avg)

        return band
