import numpy as np
import pandas as pd

from .physics import *
from he6_cres_spec_sims.constants import *

class EventBuilder:
    """  Constructs a list of betas which are trapped within the detector volume
         (Doesn't hit waveguide walls && pitch angle is magnetically trapped)
    """
    def __init__(self, config):

        self.config = config
        self.physics = Physics(config)

    def run(self):

        print("~~~~~~~~~~~~EventBuilder Block~~~~~~~~~~~~~~\n")
        print("Constructing a set of trapped events:")
        # event_num denotes the number of trapped electrons simulated.
        event_num = 0
        # beta_num denotes the total number of betas produced in the trap.
        beta_num = 0

        # if simulating full daq we instead use the beta monitor rate to determine the number of events we should be seeing
        if self.config.settings.sim_daq==True:
            events_to_simulate = self.physics.number_of_events()
            betas_to_simulate = np.inf
        else:
            events_to_simulate = self.config.physics.events_to_simulate
            betas_to_simulate = self.config.physics.betas_to_simulate

            if events_to_simulate == -1:
                events_to_simulate = np.inf
            if betas_to_simulate == -1:
                betas_to_simulate = np.inf

        print( f"Simulating: num_events:{events_to_simulate}, num_betas:{betas_to_simulate}")

        while (event_num < events_to_simulate) and (beta_num < betas_to_simulate):
            # generate trapped beta
            is_trapped = False

            while not is_trapped and beta_num < betas_to_simulate:
                if beta_num % 2500 == 0:
                    print( f"\nBetas: {beta_num}/{betas_to_simulate - 1} simulated betas.")
                    print( f"\nEvents: {event_num}/{events_to_simulate-1} trapped events.")

                initial_position, initial_direction  = self.physics.generate_beta_position_direction()
                energy = self.physics.generate_beta_energy()
                beta_num += 1

                single_event_df = self.construct_untrapped_track_df(initial_position, initial_direction, energy, event_num, beta_num)

                is_trapped = self.trap_condition(single_event_df)

            if event_num == 0:
                trapped_event_df = single_event_df

            elif beta_num == betas_to_simulate:
                break

            else:
                trapped_event_df = pd.concat([trapped_event_df, single_event_df], ignore_index=True)

            event_num += 1
        return trapped_event_df

    def construct_untrapped_track_df( self, beta_position, beta_direction, beta_energy, event_num, beta_num):
        """ Computes e.g. guiding center position, range of cyclotron radii from beta parameters
        """
        # Initial beta position and direction.
        initial_rho_pos = beta_position[0]
        initial_phi_pos = beta_position[1]
        initial_zpos = beta_position[2]

        initial_theta = beta_direction[0]
        initial_phi_dir = beta_direction[1]

        initial_field = self.config.field_strength(initial_rho_pos, initial_zpos)
        initial_radius = sc.cyc_radius(beta_energy, initial_field, initial_theta)

        # Given initial position, velocity vectors, compute guiding center position (x,y)
        # Note initial velocity vector (in x-y plane) is orthogonal to vector connecting guiding center to beta
        # \vec{r}_{GC} = \vec{r}_{init} - Rc \vec{n}_\perp, where \vec{v}_{init} \cdot \vec{n}_\perp = 0 with both unit length
        # Slightly inaccurate using Rc at beta position, and not at the guiding center (root-finding problem)
        center_x = initial_rho_pos * np.cos( initial_phi_pos / RAD_TO_DEG) - initial_radius * np.sin( initial_phi_dir / RAD_TO_DEG)
        center_y = initial_rho_pos * np.sin( initial_phi_pos / RAD_TO_DEG) + initial_radius * np.cos( initial_phi_dir / RAD_TO_DEG)

        rho_center = np.sqrt(center_x**2 + center_y**2)

        center_theta = sc.theta_center( initial_zpos, rho_center, initial_theta, self.config.trap_profile)

        # Use trapped_initial_theta to determine if trapped.
        trapped_initial_theta = sc.min_theta( rho_center, initial_zpos, self.config.trap_profile)
        max_radius = sc.max_radius( beta_energy, center_theta, rho_center, self.config.trap_profile)
        min_radius = sc.min_radius( beta_energy, center_theta, rho_center, self.config.trap_profile)

        track_properties = {
            "energy": beta_energy, #start energy
            "gamma": sc.gamma(beta_energy),
            "end_energy": 0.0,
            "initial_rho_pos": initial_rho_pos,
            "initial_phi_pos": initial_phi_pos,
            "initial_zpos": initial_zpos,
            "initial_theta": initial_theta,
            "cos_initial_theta": np.cos(initial_theta / RAD_TO_DEG),
            "initial_phi_dir": initial_phi_dir,
            "center_theta": center_theta,
            "cos_center_theta": np.cos(center_theta / RAD_TO_DEG),
            "initial_field": initial_field,
            "initial_radius": initial_radius,
            "center_x": center_x,
            "center_y": center_y,
            "rho_center": rho_center,
            "trapped_initial_theta": trapped_initial_theta,
            "max_radius": max_radius,
            "min_radius": min_radius,
            "b_avg": 0.0,
            "start_freq": 0.0,
            "end_freq": 0.0,
            "start_time": np.nan,
            "end_time": np.nan,
            "start_time_in_trap_acq": np.nan,
            "zmax": 0.0,
            "axial_freq": 0.0,
            "grad_b_freq": 0.0,
            "mod_index": 0.0,
            "track_power": 0.0,
            "slope": 0.0,
            "track_length": 0.0,
            "track_num": 0,
            "event_num": event_num,
            "beta_num": beta_num,
            "acq_num": np.nan,
            "trap_acq_num": np.nan,
        }

        event_df = pd.DataFrame(track_properties, index=[event_num])

        return event_df

    def trap_condition(self, track_df):
        """ Returns whether beta (described by track_df column row) is trapped or not
        """
        track_df = track_df.reset_index(drop=True)

        if track_df.shape[0] != 1:
            raise ValueError("trap_condition(): Input track not a single row.")

        initial_theta = track_df["initial_theta"][0]
        trapped_initial_theta = track_df["trapped_initial_theta"][0]
        rho_center = track_df["rho_center"][0]
        max_radius = track_df["max_radius"][0]
        energy = track_df["energy"][0]

        if initial_theta < trapped_initial_theta:
            # print("Not Trapped: Pitch angle too small.")
            return False

        if rho_center + max_radius > self.config.eventbuilder.decay_cell_radius:
            # print("Not Trapped: Collided with guide wall.")
            return False

        return True
