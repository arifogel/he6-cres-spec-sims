import numpy as np
import scipy.special as sp

import he6_cres_spec_sims.spec_tools.spec_calc.spec_calc as sc

from he6_cres_spec_sims.constants import *

def power_calc(center_x, center_y, frequency, field, trap_radius):

    """Calculates the average cyclotron radiation power (in one direction) in Watts in the
    TE11 mode of an electron undergoing cyclotron motion in the
    cylindrical waveguide around the point (center_x,center_y) with
    frequency in Hz and field in Tesla.
    See https://arxiv.org/abs/2405.06847 Eqs 18-20 (n,m,h=1)
    """

    center_rho = np.sqrt(center_x**2 + center_y**2)

    kc = P11_PRIME / trap_radius
    Rcycl = sc.cyc_radius(sc.freq_to_energy(frequency, field), field, 90)

    # values in power equation
    omega = 2 * PI * frequency
    k = omega / C
    v_perp = Rcycl * omega

    beta = np.sqrt(pow(k, 2) - pow(kc, 2))

    P_lambda = PI * beta / (2*kc**2 *MU_0 * omega) * (P11_PRIME**2 - 1) * sp.jv(1,P11_PRIME)**2
    power = (Q*v_perp/2.) **2 / P_lambda * sp.jv(0, kc*center_rho)**2 * sp.jvp(1, kc*Rcycl)**2

    # To vectorize (for speed), we don't want an if statement, just remove Nans for cutoff frequencies
    return np.nan_to_num(power, copy=False, nan=0, posinf=0)
