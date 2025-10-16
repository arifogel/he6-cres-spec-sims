"""spec_calc

This set of calculators takes cres electron properties (such as energy, pitch angle), as
well as information about the trap (instance of a trap_profile object) and outputs other
cres electron properties such as axial frequency, z_max, B_avg.

---
Units for all module functions' inputs/outputs:

Energy   : eV
B-field  : T
Time     : s
Angle    : degrees
Distance : m
Frequency: Hz
Power    : W
---
"""
import numpy as np

import scipy.integrate as integrate
from scipy.fft import fft
from scipy.optimize import root_scalar
from scipy.misc import derivative
from scipy.special import jv

from he6_cres_spec_sims.constants import *

# Simple special relativity functions.


def gamma(energy):
    gamma = (energy + ME) / ME
    return gamma

def energy(gamma):
    energy = ME*(gamma - 1)
    return energy

def momentum(energy):
    momentum = (((energy + ME) ** 2 - ME**2) / C**2) ** 0.5 / J_TO_EV
    return momentum


def velocity(energy):
    velocity = momentum(energy) / (gamma(energy) * M)
    return velocity


# CRES functions.
def energy_to_freq(energy, field):

    """Converts kinetic energy to cyclotron frequency."""

    cycl_freq = Q * field / (2 * PI * gamma(energy) * M)

    return cycl_freq


def freq_to_energy(frequency, field):

    """Calculates energy of beta particle in eV given cyclotron
    frequency in Hz, magnetic field in Tesla, and pitch angle
    at 90 degrees.
    """

    gamma = Q * field / (2 * PI * frequency * M)
    if np.any(gamma < 1):
        gamma = 1
        max_freq = Q * field / (2 * PI * M)
        warning = "Warning: {} higher than maximum cyclotron frequency {}".format(
            frequency, max_freq
        )
        print(warning)
    return gamma * ME - ME


def energy_and_freq_to_field(energy, freq):

    """Converts kinetic energy to cyclotron frequency."""

    field = (2 * PI * gamma(energy) * M * freq) / Q

    return field


def power_from_slope(energy, slope, field):

    """Converts slope, energy, field into the associated cres power."""

    energy_Joules = (energy + ME) / J_TO_EV

    power = slope * (2 * PI) * ((energy_Joules) ** 2) / (Q * field * C**2)

    return power

def theta_center(zpos, rho, pitch_angle, trap_profile):

    """Calculates the pitch angle an electron with current z-coordinate
    zpos, rho, and current pitch angle pitch_angle takes at the center
    of given trap.
    """

    if trap_profile.is_trap:

        Bmin = trap_profile.field_strength(rho, 0.0)
        Bcurr = trap_profile.field_strength(rho, zpos)

        theta_center_calc =  np.arcsin((np.sqrt(Bmin / Bcurr)) * np.sin(pitch_angle / RAD_TO_DEG)) * RAD_TO_DEG

        return theta_center_calc

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def cyc_radius(energy, field, pitch_angle):

    """Calculates the instantaneous cyclotron radius of a beta electron
    given the energy, magnetic field, and current pitch angle.
    """

    vel_perp = velocity(energy) * np.sin(pitch_angle / RAD_TO_DEG)

    cyc_radius = (gamma(energy) * M * vel_perp) / (Q * field)

    return cyc_radius


def max_radius(energy, center_pitch_angle, rho, trap_profile):

    """Calculates the maximum cyclotron radius of a beta electron given
    the kinetic energy, trap_profile, and center pitch angle (pitch angle
    at center of trap).
    """

    if trap_profile.is_trap:
        min_field = trap_profile.field_strength(rho, 0)
        max_field = min_field / (np.sin(center_pitch_angle / RAD_TO_DEG)) ** 2

        center_radius = cyc_radius(energy, min_field, center_pitch_angle)
        end_radius = cyc_radius(energy, max_field, pitch_angle=90)

        if np.all(center_radius >= end_radius):
            return center_radius
        else:
            print(
                "Warning: max_radius is occuring at end of trap (theta=90). \
                Something odd may be going on."
            )
            return False

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def min_radius(energy, center_pitch_angle, rho, trap_profile):

    """Calculates the minimum cyclotron radius of a beta electron given
    the kinetic energy, trap_profile, and center pitch angle (pitch angle
    at center of trap).
    """

    if trap_profile.is_trap:

        min_field = trap_profile.field_strength(rho, 0)
        max_field = min_field / (np.sin(center_pitch_angle / RAD_TO_DEG)) ** 2

        center_radius = cyc_radius(energy, min_field, center_pitch_angle)
        end_radius = cyc_radius(energy, max_field, pitch_angle=90)

        if np.all(center_radius >= end_radius):
            return end_radius
        else:
            print(
                "Warning: min_radius is occuring at center of trap, something \
                odd may be going on."
            )
            return False

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def min_theta(rho, zpos, trap_profile):

    """Calculates the minimum pitch angle theta at which an electron at
    zpos is trapped given trap_profile.
    """
    if trap_profile.is_trap:
        # Be careful here. Technically the Bmax doesn't occur at a constant z.
        Bmax = trap_profile.field_strength(rho, trap_profile.trap_width[1])
        Bz = trap_profile.field_strength(rho, zpos)

        theta = np.arcsin((Bz / Bmax) ** 0.5) * RAD_TO_DEG
        return theta

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


@np.vectorize
def max_zpos(energy, center_pitch_angle, rho, trap_profile, debug=False):

    """Calculates the maximum axial length from center of trap as a function of center_pitch_angle and rho.
       One would ideally prefer minimization which was not looped (with np.vectorize)
       In practice, library multi-dimensional minimizers are slower, even if Jacobian is diagonal
    """

    if trap_profile.is_trap:

        if center_pitch_angle < min_theta(rho, 0, trap_profile):
            print("WARNING: Electron not trapped (zpos)")
            return False

        else:
            # Ok, so does this mean we now have an energy dependence on zmax? Yes.
            c_r = cyc_radius( energy, trap_profile.field_strength(rho, 0), center_pitch_angle)
            rho_p = np.sqrt(rho**2 + c_r**2 / 2)

            min_field = trap_profile.field_strength(rho_p, 0)
            max_field = trap_profile.field_strength(rho_p, trap_profile.trap_width[1])

            max_reached_field = min_field / pow( math.sin(center_pitch_angle / RAD_TO_DEG), 2)

            func = lambda z: trap_profile.field_strength(rho_p, z) - max_reached_field

            # root-finding generally easier than minimization (trivial to step z +- dz by sign of func(z))
            solution = root_scalar(func, bracket=[0, trap_profile.trap_width[1]], method='brentq', xtol=1e-14)
            max_z = solution.root
            curr_field = trap_profile.field_strength(rho_p, max_z)

            if debug and (curr_field > max_reached_field):
                print( "Final field greater than max allowed field by: ", curr_field - max_reached_field)
                print("Bmax reached: ", curr_field)

            if debug == True:
                print("zlength: ", max_z)

            if max_z > trap_profile.trap_width[1]:
                print("Error Rate: ", max_z - trap_profile.trap_width[1])

            return max_z

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False

def mod_index(avg_cycl_freq, zmax):

    """Calculates modulation index from average cyclotron frequency
    (avg_cycl_freq) and maximum axial amplitude (zmax).
    """
    # calculated parameters
    omega = 2 * PI * avg_cycl_freq
    beta = waveguide_beta(omega)

    mod_index = zmax * beta

    return mod_index


def df_dt(energy, field, power):

    """Calculates cyclotron frequency rate of change of electron with
    given kinetic energy at field in T radiating energy at rate power.
    """

    energy_Joules = (energy + ME) / J_TO_EV

    slope = (Q * field * C**2) / (2 * PI) * (power) / (energy_Joules) ** 2

    return slope


def curr_pitch_angle(rho, zpos, center_pitch_angle, trap_profile):

    """Calculates the current pitch angle of an electron at zpos given
    center pitch angle and main field strength.
    """

    if trap_profile.is_trap:

        min_field = trap_profile.field_strength(rho, 0)
        max_z = max_zpos(center_pitch_angle, rho, trap_profile)
        max_reached_field = trap_profile.field_strength(rho, max_z)

        if np.any(abs(zpos) > max_z):
            print("Electron does not reach given zpos")
            curr_pitch = "FAIL"
        else:
            curr_field = trap_profile.field_strength(rho, zpos)
            curr_pitch = np.arcsin(np.sqrt(curr_field / max_reached_field)) * RAD_TO_DEG

        return curr_pitch

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False

def semiopen_simpson(v):
    """Semi-Open Composite Simpson's Rule for fast/ vectorized integral evaluation
       Numerical Recipes, pg 162. Requires you to set v[0] = v[1] previously
    """
    return np.sum(v,axis=0)  - v[1] * 1./12 - v[2] * 5./12 - v[-1] * 7./12 + v[-2]  * 1./12

def axial_freq(energy, center_pitch_angle, rho, trap_profile, nIntegralPoints=200):
    """Calculates the axial frequency of trapped electrons."""
    if trap_profile.is_trap:
        # axial_freq of 90 deg otherwise returns 1./0. Would prefer to get correct limit
        # Sets allowed range, "clipping" the pitches at 90 deg.
        center_pitch_angle = np.clip(center_pitch_angle, 0, 89.9999)

        # Get cyclotron-radius modified "effective" guiding center rho
        c_r = cyc_radius( energy, trap_profile.field_strength(rho, 0), center_pitch_angle)
        rho_p = np.sqrt(rho**2 + c_r**2 / 2)

        # Field at center of trap
        B0 = trap_profile.field_strength(rho_p, 0)
        # Field at turning point
        Bmax = B0 / np.sin(center_pitch_angle / RAD_TO_DEG)**2

        B = lambda z: trap_profile.field_strength(rho_p, z)

        # Should optionally pass these in as argument to reuse calculations!
        zmax = max_zpos(energy, center_pitch_angle, rho, trap_profile)

        # See write-ups XXX for more information on this integral
        u = np.linspace(0,1., nIntegralPoints)
        du = u[1]
        # Semi-open simpsons rule avoids evaluation at t=0. Just replace with next entry (semi-open)
        u[0] = u[1]

        zmax_arr = np.atleast_1d(np.array(zmax))
        Bmax_arr = np.atleast_1d(np.array(Bmax))

        u = u[:,np.newaxis]
        zmax_arr = zmax_arr[np.newaxis,:]
        Bmax_arr = Bmax_arr[np.newaxis,:]

        integrand = u / np.sqrt(1. - B(zmax_arr *(1.-u**2)) / Bmax_arr)

        T_a = 8. * zmax / velocity(energy) * semiopen_simpson(integrand) * du

        axial_frequency = 1. / T_a

        return axial_frequency

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False

def avg_cycl_freq(energy, center_pitch_angle, rho, trap_profile):
    field = b_avg(energy, center_pitch_angle, rho, trap_profile)
    return energy_to_freq(energy, field)

def b_avg(energy, center_pitch_angle, rho, trap_profile, ax_freq=None, nIntegralPoints=200):

    """Calculates the average magnetic field experienced by an array of electrons
    with given kinetic energies, main fields, and center pitch angles.
    Returns 0 if not trapped.
    """
    if trap_profile.is_trap:

        # For convenience, to avoid if statements for 90 deg. (stationary in z)
        center_pitch_angle = np.clip(center_pitch_angle, 0, 89.9999)

        # we should not be redoing integrals which are the bottleneck of the Monte Carlo, by default
        if ax_freq is None:
            ax_freq = axial_freq(energy, center_pitch_angle, rho, trap_profile)

        c_r = cyc_radius(energy, trap_profile.field_strength(rho, 0), center_pitch_angle)

        rho_p = np.sqrt(rho**2 + c_r**2 / 2)
        rho_pp = np.sqrt(rho**2 + c_r**2)

        # Field at center of trap
        B0 = trap_profile.field_strength(rho_p, 0)
        # Field at turning point
        Bmax = B0 / np.sin(center_pitch_angle / RAD_TO_DEG)**2
        # Field in between (vs. instantaneous pitch angle)
        Bp = lambda z: trap_profile.field_strength(rho_p, z)
        Bpp = lambda z: trap_profile.field_strength(rho_pp, z)

        # Should optionally pass these in as argument to reuse calculations!
        zmax = max_zpos(energy, center_pitch_angle, rho, trap_profile)
        # Should optionally pass these in as argument to reuse calculations!
        f_a = axial_freq(energy, center_pitch_angle, rho, trap_profile, nIntegralPoints)

        # See write-ups XXX for more information on this integral
        u = np.linspace(0,1., nIntegralPoints)
        du = u[1]
        # Semi-open simpsons rule avoids evaluation at t=0. Just replace with next entry (semi-open)
        u[0] = u[1]

        zmax_arr = np.atleast_1d(np.array(zmax))
        Bmax_arr = np.atleast_1d(np.array(Bmax))

        u = u[:,np.newaxis]
        zmax_arr = zmax_arr[np.newaxis,:]
        Bmax_arr = Bmax_arr[np.newaxis,:]

        integrand = u * Bpp(zmax_arr*(1-u**2)) / np.sqrt(1. - Bp(zmax_arr*(1.-u**2)) / Bmax_arr)

        b_avg = 8. * zmax * f_a / velocity(energy) * semiopen_simpson(integrand) * du

        return b_avg

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False


def grad_b_freq(energy, center_pitch_angle, rho, trap_profile, ax_freq=None, nIntegralPoints=200):
    """Calculates the axial frequency of trapped electrons."""
    if trap_profile.is_trap:
        # axial_freq of 90 deg otherwise returns 1./0. Would prefer to get correct limit
        # Sets allowed range, "clipping" the pitches at 90 deg.
        center_pitch_angle = np.clip(center_pitch_angle, 0, 89.9999)

        # Get cyclotron-radius modified "effective" guiding center rho
        #c_r = cyc_radius( energy, trap_profile.field_strength(rho, 0), center_pitch_angle)
        #rho_p = np.sqrt(rho**2 + c_r**2 / 2)

        if ax_freq is None:
            ax_freq = axial_freq(energy, center_pitch_angle, rho, trap_profile)

        # Field at center of trap
        B0 = trap_profile.field_strength(rho, 0)
        # Field at turning point
        Bmax = B0 / np.sin(center_pitch_angle / RAD_TO_DEG)**2

        B = lambda z: trap_profile.field_strength(rho, z)

        rho_delta = 1e-7
        dBdRho = lambda z: (trap_profile.field_strength(rho + rho_delta, z) - trap_profile.field_strength(rho - rho_delta, z)) / (2. * rho_delta)

        # Should optionally pass these in as argument to reuse calculations!
        zmax = max_zpos(energy, center_pitch_angle, rho, trap_profile)

        # See write-ups XXX for more information on this integral
        u = np.linspace(0,1., nIntegralPoints)
        du = u[1]
        # Semi-open simpsons rule avoids evaluation at t=0. Just replace with next entry (semi-open)
        u[0] = u[1]

        zmax_arr = np.atleast_1d(np.array(zmax))
        Bmax_arr = np.atleast_1d(np.array(Bmax))

        u = u[:,np.newaxis]
        zmax_arr = zmax_arr[np.newaxis,:]
        Bmax_arr = Bmax_arr[np.newaxis,:]

        z_arg = zmax_arr *(1.-u**2)
        integrand = u * (2 - B(z_arg) / Bmax_arr) * dBdRho(z_arg)  / (np.sqrt(1. - B(z_arg) / Bmax_arr) * B(z_arg)**2)

        ### Energy == KINETIC ENERGY (γ-1) m c**2. Multiply by Q because energy is in eV, not Joules
        grad_B_frequency = 4 / PI * (energy * zmax * ax_freq) / (rho  * velocity(energy)) * semiopen_simpson(integrand) * du

        #We don't really care which direction it goes: just want to report a frequency
        grad_B_frequency = np.abs(grad_B_frequency)

        return grad_B_frequency

    else:
        print("ERROR: Given trap profile is not a valid trap")
        return False



def waveguide_beta(omega):
    """  Computes the (waveguide definition) of beta (propagation constant for TE11 mode
    """
    # fixed experiment parameters
    waveguide_radius = 0.578e-2
    kc = P11_PRIME / waveguide_radius

    # calculated parameters
    k_wave = omega / C
    beta = np.sqrt(k_wave**2 - kc**2)
    return beta

def sideband_calc(energy, rho, avg_cycl_freq, axial_freq, zmax, trap_profile, magnetic_modulation=True, num_sidebands=7, nHarmonics=128):

    """Calculates relative magnitudes of num_sidebands sidebands from
    average cyclotron frequency (avg_cycl_freq), axial frequency
    (axial_freq), and maximum axial amplitude (zmax).
    """
    if magnetic_modulation:
        T = 1./ axial_freq
        dt = T/ nHarmonics
        t = np.arange(0,T,dt)
        omegaA = 2. * np.pi * axial_freq

        z = zmax*np.sin(omegaA * t)
        vz = zmax*omegaA*np.cos(omegaA * t)

        sidebands =  FFT_sideband_amplitudes(energy, rho, avg_cycl_freq, axial_freq, vz, z, trap_profile, True, nHarmonics)
        return format_sideband_array(sidebands, avg_cycl_freq, axial_freq, np.nan, num_sidebands)

    else:
        h =  mod_index(avg_cycl_freq, zmax)
        sidebands = [abs(jv(k, h)) for k in range(num_sidebands+1)]
        return format_sideband_array(sidebands, avg_cycl_freq, axial_freq, h, num_sidebands)



def anharmonic_sideband_calc(energy, center_pitch_angle, rho, avg_cycl_freq, axial_freq, zmax, trap_profile, magnetic_modulation=True, num_sidebands=7, nHarmonics=128):

    #c_r = cyc_radius(energy, trap_profile.field_strength(rho, 0), center_pitch_angle)
    #rho= np.sqrt(rho**2 + c_r**2 / 2)

    ### Compute particle trajectory over single period
    sol = anharmonic_axial_trajectory(energy, center_pitch_angle, rho, axial_freq, zmax, trap_profile, nHarmonics)
    z = sol[0]
    vz = sol[1]

    sidebands =  FFT_sideband_amplitudes(energy, rho, avg_cycl_freq, axial_freq, vz, z, trap_profile, magnetic_modulation, nHarmonics)
    return format_sideband_array(sidebands, avg_cycl_freq, axial_freq, np.nan, num_sidebands)


def anharmonic_axial_trajectory(energy, center_pitch_angle, rho, axial_freq, zmax, trap_profile, nHarmonics):
    """ Computes the time series of the beta axial motion over a single
    found by integrating the relevant ODE. Returns [z(t), vz(t)].
    """
    if not trap_profile.is_trap:
        print("ERROR: Given trap profile is not a valid trap")
        return False

    T = 1./ axial_freq
    dt = T/ nHarmonics
    t = np.arange(0,T,dt)

    p0 = M * velocity(energy)
    ### Note: This is the non-relativistic magnetic moment. One gets the same ODE if M -> gamma M in the lambda ode.
    Bmin = trap_profile.field_strength(rho, 0)
    mu = p0**2 * np.sin(center_pitch_angle / RAD_TO_DEG )**2 / (2. * M * Bmin)
    Bz = lambda z: trap_profile.field_strength(rho,z)
    dBdz = lambda z: derivative(Bz, z, dx=1e-6)
    ### Coupled ODE for z-motion: z = y[0], vz = y[1]. z'=vz. vz' = -mu * B'(z) / m
    ode = lambda t, y: [y[1], - mu / M * dBdz(y[0])]
    result = integrate.solve_ivp(ode, [t[0], t[-1]], (zmax, 0), t_eval=t,rtol=1e-7)

    ####### [0] is z array, [1] is vz array ######
    return result.y

def instantaneous_frequency(energy, rho, avg_cycl_freq, vz, z=None, trap_profile=None, magnetic_modulation=False):
    """ Computes the instantaneous (angular) frequency as a function of time. Magnetic modulation controls whether
        instantaneous changes in magnetic field are included, or just the Doppler effect. Defaults chosen so one could pass
        fewer arguments if not using magnetic_modulation
    """
    omega = 2 * PI * avg_cycl_freq
    beta = waveguide_beta(omega)
    phase_vel = omega / beta

    if magnetic_modulation:
        if not trap_profile.is_trap:
            print("ERROR: Given trap profile is not a valid trap")
            return False
        Bz = lambda z: trap_profile.field_strength(rho, z)
        return Q * Bz(z) / (M * gamma(energy)) * ( 1. + vz / phase_vel)
    else:
        ### Instead of evaluating at B(0), use average_cyc_freq to say B field contribution is const. equal to avg
        return omega * ( 1. + vz / phase_vel)


def FFT_sideband_amplitudes(energy, rho, avg_cycl_freq, axial_freq, vz, z, trap_profile, magnetic_modulation, nHarmonics=128):
    """  Computes sideband amplitudes as a function of axial trajectory, magnetic field profile. Returns list with sidebands
    """
    ### Convert particle trajectory to instantaneous radiated frequency (Doppler + B-field)
    dt = 1./ (nHarmonics * axial_freq)

    omega_c = instantaneous_frequency(energy, rho, avg_cycl_freq, vz, z, trap_profile, magnetic_modulation)
    omega_c -= np.mean(omega_c)
    Phi = np.cumsum(omega_c) * dt
    expPhi = np.exp(1j * Phi)
    yf = np.abs(fft(expPhi,norm="forward"))
    yf = yf[:nHarmonics//2]
    return yf

def format_sideband_array(sidebands_one, avg_cyc_freq, axial_freq, mod_index=np.nan, num_sidebands = 7):
    """ Does formatting for array with list of sideband magnitudes (normalized), and their start frequencies.
        Takes in 1-sided list of sideband magnitudes
    """
    # Calculate (2-sided) list of (frequency, amplitude) of sidebands
    sidebands = []

    for k in range(-num_sidebands, num_sidebands + 1):
        freq = avg_cyc_freq + k * axial_freq
        magnitude = sidebands_one[abs(k)]
        pair = (freq, magnitude)
        sidebands.append(pair)

    ### Intentionally returns modulation index of nan as it is only (meaningfully) defined for harmonic traps
    return sidebands, mod_index

def power_larmor(field, frequency):

    omega = 2 * PI * frequency
    energy = freq_to_energy(frequency, field)
    r_c = cyc_radius(energy, field, 90)
    beta = velocity(energy) / C
    p = gamma(energy) * M * velocity(energy)

    power_larmor = (2 / 3 * Q**2 * C * beta**4 * gamma(energy) ** 4) / (
        4 * PI * EPS_0 * r_c**2
    )

    return power_larmor

def power_larmor_e(field, energy):
    """Takes energy instead of frequency as input. """

    r_c = cyc_radius(energy, field, 90)
    beta = velocity(energy) / C
    p = gamma(energy) * M * velocity(energy)

    power_larmor = (2 / 3 * Q**2 * C * beta**4 * gamma(energy) ** 4) / (
        4 * PI * EPS_0 * r_c**2
    )

    return power_larmor
