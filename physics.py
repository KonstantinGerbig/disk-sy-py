import numpy as np
from constants import *


def f_Omega(R, M_star):
    return np.sqrt(Grav * M_star/ R**3)


def f_temperature(R, T_star, q, irr_angle):
    return T_star * (R_sun/R)**q * irr_angle**0.25

def f_sound_speed(temp):
    return np.sqrt(k_b * temp / (mu*m_p))


def f_viscosity(alpha, c_s, Omega):
    return alpha * c_s**2 / Omega

def f_similarity_solution(r, t, scaling_C, nu_1, gamma, t_s, l_gasevolution):
    if l_gasevolution == True:
        T = t/t_s + 1
    else:
        T = 1
    # Note that r is dimensionless: r = R/R_1
    prefactor = scaling_C / (3 * pi * nu_1)
    return prefactor * r**(-gamma) * T**(-((5./2)-gamma)/(2 - gamma)) * np.exp(-(r**(2-gamma))/T)

def f_pressure_grad(r, t, gamma, q, t_s, l_gasevolution):
    if l_gasevolution == True:
        T = t/t_s + 1
    else:
        T = 1
    return np.abs((0.5*q - gamma - 1.5)- (2-gamma) * r**(2-gamma) / T)


def f_growth_timescale(dtogr, St_pbb, St_dst, Omega):
    return np.log(St_pbb / St_dst) / (dtogr * Omega)

def f_Ep_Stokes(Sigma_g, dst_size, dst_dens):
    # assume dst_size = const. 
    # assume dust always in Epstein regime
    return pi/2 * dst_size * dst_dens / Sigma_g

def f_St_frag(frag_velocity, alpha, c_s):
    return  (1./3) * frag_velocity**2 / (alpha * c_s**2)

def f_St_drift(dtogr, pressure_grad, v_k, c_s):
    # stuff drifts away
    fudge_drift = 1
    return fudge_drift * dtogr * v_k**2 / (c_s**2 * pressure_grad)

def f_St_df(pressure_grad, v_k, c_s, frag_velocity):
    # differential drift
    N_df = 0.5 # equal sized collisions
    return frag_velocity * v_k / (pressure_grad * c_s**2 * (1 - N_df))

def f_St_pbb(dtogr, pressure_grad, frag_velocity, alpha, c_s, v_k, St_dst):
    St_frag = f_St_frag(frag_velocity, alpha, c_s)
    St_drift = f_St_drift(dtogr, pressure_grad, v_k, c_s)
    St_df = f_St_df(pressure_grad, v_k, c_s, frag_velocity)    
    St_pbb = np.amin([St_frag, St_drift, St_df])
    # prevent Stokes number that is smaller than dust stokes number
    St_pbb = np.amax([St_pbb, 2*St_dst])
    return St_pbb, St_drift, St_df, St_frag


def f_v_pbb(St, pressure_grad, dtogr, H, R, c_s):
    return (St / (St**2 + (1+dtogr)**2)) * (H/R) * pressure_grad * c_s


def f_escape_velocity(mass, radius):
    return np.sqrt(2*Grav*mass/radius)

def f_hill_velocity(mass, R, Omega, M_star):
    return R*Omega * (2*mass/(3*M_star))**(1./3)

def f_collision_cross_section(radius, v_esc, v_rel):
    return 4*pi *radius**2 * (1+ (v_esc/v_rel)**2)

def f_coll_timescale(Sigma_pls, m_pls, Omega, r_pls, R, M_star):
    v_hill = f_hill_velocity(m_pls, R, Omega, M_star)
    v_esc = f_escape_velocity(m_pls, r_pls)
    v_rel = v_hill
    sigma_coll_cross_sect = f_collision_cross_section(r_pls, v_esc, v_rel)
    collision_fudge = 1
    if not Sigma_pls == 0:
        coll_timescale = collision_fudge* np.sqrt(2*pi) * m_pls * v_hill / (Sigma_pls * Omega * sigma_coll_cross_sect * v_rel)
    else:
        coll_timescale = 10**10 * year # some large number
    return coll_timescale


def f_fragment_fraction(m_dst, m_pls, a_small, a_big, xi, pls_dens):
    m_1 = (4/3)*pi * pls_dens * (a_small/2)**3
    m_2 = (4/3)*pi * pls_dens * (a_big/2)**3
    
    frac = (m_2**(2-xi) - m_1**(2-xi))/(m_pls**(2-xi) - m_dst**(2-xi))
    return frac



class dust_growth():

    def __init__(self, R, VAR, CONST, PHYS, t):

        """
            Load current surface densities
        """
        
        sig_g = VAR.sig_g

        sig_dst = VAR.sig_dst
        sig_pbb = VAR.sig_pbb

        sig_dbr_total = VAR.sig_dbr_um + VAR.sig_dbr_mm + VAR.sig_dbr_dm +VAR.sig_dbr_m
        self.dtogr = (sig_dbr_total + sig_dst + sig_pbb)/sig_g
        
        # dust Stokes number
        self.st_dst = f_Ep_Stokes(sig_g, PHYS.dst_size, PHYS.dst_dens)

        # pressure gradient
        self.pressure_grad = f_pressure_grad(R/PHYS.R_1, t, PHYS.gamma, PHYS.q, CONST.t_s, PHYS.l_gasevolution)

        # pebble Stokes number
        self.st_pbb, self.st_drift, self.st_df, self.st_frag = f_St_pbb(self.dtogr, self.pressure_grad, PHYS.frag_velocity, PHYS.alpha, CONST.c_s, CONST.v_k, self.st_dst)

        # get dust growth timescale and flux
        if PHYS.l_dustgrowth == True:
            self.timescale_growth = f_growth_timescale(self.dtogr, self.st_pbb, self.st_dst, CONST.Omega)
            self.dust_growth_flux = sig_dst / self.timescale_growth
        else:
            self.timescale_growth = np.NaN
            self.dust_growth_flux = float(0)

class formation():

    def __init__(self, R, VAR, CONST, PHYS, DUST_GROWTH):

        """
            Load current surface densities
        """
        sig_pbb = VAR.sig_pbb
        


        self.v_pbb = f_v_pbb(DUST_GROWTH.st_pbb, DUST_GROWTH.pressure_grad, DUST_GROWTH.dtogr, CONST.H, R, CONST.c_s)


        if PHYS.l_formation == True:
            self.conversion_length =  (PHYS.trapping_distance * CONST.H) / PHYS.trapping_efficiency
            self.pls_form_flux_first = sig_pbb * self.v_pbb / self.conversion_length
            
            if PHYS.l_formation_second_generation == True:
                
                sig_dbr_dm = VAR.sig_dbr_dm
                sig_dbr_m = VAR.sig_dbr_m

                self.pls_form_flux_second = (sig_dbr_dm + sig_dbr_m) * self.v_pbb / self.conversion_length
                self.pls_form_flux_dm = sig_dbr_dm * self.v_pbb / self.conversion_length
                self.pls_form_flux_m = sig_dbr_m * self.v_pbb / self.conversion_length
            else:
                self.pls_form_flux_second = float(0)
                self.pls_form_flux_dm = float(0)
                self.pls_form_flux_m = float(0)
        else:
            self.conversion_length = np.NaN
            self.pls_form_flux_first = float(0)
            self.pls_form_flux_second = float(0)
            self.pls_form_flux_dm = float(0)
            self.pls_form_flux_m = float(0)

class collisions():
    
    def __init__(self, R, VAR, CONST, PHYS):

        """
            Load current surface densities
        """

        sig_pls_first = VAR.sig_pls_first
        sig_pls_second = VAR.sig_pls_second

        sig_pls_total = sig_pls_first + sig_pls_second

        

        if PHYS.l_collisions == True:
            """
            Load Fragment fractions
            """        

            frac_um = f_fragment_fraction(CONST.m_dst, CONST.m_pls, PHYS.a_0, PHYS.a_1, PHYS.xi, PHYS.pls_dens)
            frac_mm = f_fragment_fraction(CONST.m_dst, CONST.m_pls, PHYS.a_1, PHYS.a_2, PHYS.xi, PHYS.pls_dens)
            frac_dm = f_fragment_fraction(CONST.m_dst, CONST.m_pls, PHYS.a_2, PHYS.a_3, PHYS.xi, PHYS.pls_dens)
            frac_Dm = f_fragment_fraction(CONST.m_dst, CONST.m_pls, PHYS.a_3, PHYS.a_4, PHYS.xi, PHYS.pls_dens)
            frac_km = f_fragment_fraction(CONST.m_dst, CONST.m_pls, PHYS.a_4, PHYS.a_5, PHYS.xi, PHYS.pls_dens)

            """
            Get collision fractions
            """   

            self.timescale_collision = f_coll_timescale(sig_pls_total, CONST.m_pls, CONST.Omega, PHYS.r_pls, R, PHYS.M_star)
            self.collision_flux_pls_first = (frac_km - 1) * sig_pls_first/self.timescale_collision
            
            self.collision_flux_um = frac_um* sig_pls_total/self.timescale_collision
            self.collision_flux_mm = frac_mm* sig_pls_total/self.timescale_collision
            self.collision_flux_dm = frac_dm* sig_pls_total/self.timescale_collision
            self.collision_flux_m = frac_Dm* sig_pls_total/self.timescale_collision

            self.collision_flux_pls_second = (frac_km - 1)* sig_pls_second/self.timescale_collision
        else:
            self.timescale_collision = np.NaN
            self.collision_flux_pls_first = float(0)

            self.collision_flux_um = float(0)
            self.collision_flux_mm = float(0)
            self.collision_flux_dm = float(0)
            self.collision_flux_m = float(0)

            self.collision_flux_pls_second = float(0)

class drift():

    def __init__(self, R, VAR,  CONST, PHYS, DUST_GROWTH):
        """
            Load current surface densities
        """
        sig_pbb = VAR.sig_pbb
        sig_dbr_mm = VAR.sig_dbr_mm
        sig_dbr_dm = VAR.sig_dbr_dm
        sig_dbr_m = VAR.sig_dbr_m

        if PHYS.l_drift == True:

            # PEBBLES
            v_pbb = f_v_pbb(DUST_GROWTH.st_pbb, DUST_GROWTH.pressure_grad, DUST_GROWTH.dtogr, CONST.H, R, CONST.c_s)
            timescale_drift_pbb = R / v_pbb
            self.drift_flux_pbb = sig_pbb/timescale_drift_pbb

            # DEBRIS

            a_mm = 0.5*(PHYS.a_1 + PHYS.a_2)
            St_mm = f_Ep_Stokes(VAR.sig_g, a_mm, PHYS.dst_dens)
            v_dbr_mm = f_v_pbb(St_mm, DUST_GROWTH.pressure_grad, DUST_GROWTH.dtogr, CONST.H, R, CONST.c_s)
            timescale_drift_dbr_mm = R / v_dbr_mm
            self.drift_flux_dbr_mm = sig_dbr_mm / timescale_drift_dbr_mm

            a_dm = 0.5*(PHYS.a_2 + PHYS.a_3)
            St_dm = f_Ep_Stokes(VAR.sig_g, a_dm, PHYS.dst_dens)
            v_dbr_dm = f_v_pbb(St_dm, DUST_GROWTH.pressure_grad, DUST_GROWTH.dtogr, CONST.H, R, CONST.c_s)
            timescale_drift_dbr_dm = R / v_dbr_dm
            self.drift_flux_dbr_dm = sig_dbr_dm / timescale_drift_dbr_dm

            a_m = 0.5*(PHYS.a_3 + PHYS.a_4)
            St_m = f_Ep_Stokes(VAR.sig_g, a_m, PHYS.dst_dens)
            v_dbr_m = f_v_pbb(St_m, DUST_GROWTH.pressure_grad, DUST_GROWTH.dtogr, CONST.H, R, CONST.c_s)
            timescale_drift_dbr_m = R / v_dbr_m
            self.drift_flux_dbr_m = sig_dbr_m / timescale_drift_dbr_m

            
        else:
            self.drift_flux_pbb = float(0)
            self.drift_flux_dbr_mm = float(0)
            self.drift_flux_dbr_dm = float(0)
            self.drift_flux_dbr_m = float(0)


class gas_evolution():

    def __init__(self, R, t, dt, CONST, PHYS):

        self.sig_g = f_similarity_solution(R/PHYS.R_1, t+dt, CONST.scaling_C, CONST.nu_1, PHYS.gamma, CONST.t_s, PHYS.l_gasevolution)