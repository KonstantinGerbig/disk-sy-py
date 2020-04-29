import numpy as np
from constants import *
from physics import *
import pickle
import os
from mpl_style import *

class data():
    def __init__(self, VAR):
        import inspect
        attributes_all = inspect.getmembers(VAR, lambda x:not(inspect.isroutine(x)))
        attributes = [x for x in attributes_all if not(x[0].startswith('__') and x[0].endswith('__'))]
        for k,v in attributes:
             setattr(self, k, v)

    def append(self, VAR):
        import inspect
        import warnings
        attributes_all = inspect.getmembers(VAR, lambda x:not(inspect.isroutine(x)))
        attributes = [x for x in attributes_all if not(x[0].startswith('__') and x[0].endswith('__'))]
        for k,v in attributes:
            if hasattr(self, k):
                setattr(self, k, np.append(getattr(self, k), v))
            else:
                warnings.warn("Attribute not found")

    def write_to_file(self, R, diskpath):
        import inspect
        attributes_all = inspect.getmembers(self, lambda x:not(inspect.isroutine(x)))
        attributes = [x for x in attributes_all if not(x[0].startswith('__') and x[0].endswith('__'))]
        filename = str(np.round(R/AU, 2))
        filename = filename.replace('.', '_') + '.pkl'
        pickle.dump(attributes, open(diskpath + '/local/' + filename, 'wb'))


class var():
    """
        List of quantities that can be tracked
    """
    t = None
    dt = None
    
    sig_g = None
    sig_dst = None
    sig_pbb = None
    sig_pls_first = None
    sig_pls_second = None
    
    sig_dbr_um = None 
    sig_dbr_mm = None
    sig_dbr_dm = None
    sig_dbr_m = None

    sig_drift_loss = None
    
    #st_dst = None
    #st_frag = None
    #st_drift = None
    #st_df = None
    
    #v_pbb = None
    
    
    
    def __init__(self, **kwargs):
        """
        Write all in var object
        """
        import warnings
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                warnings.warn("Attribute not found")
    

    def renew(self, **kwargs):
        import warnings
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                warnings.warn("Attribute not found")



class var_nonevolving():
    def __init__(self, R, PHYS):
        self.Omega_1 = f_Omega(PHYS.R_1, PHYS.M_star)

    
        self.temp_1 = f_temperature(PHYS.R_1, PHYS.T_star, PHYS.q, PHYS.irr_angle)
        self.c_s_1 = f_sound_speed(self.temp_1)
        
        self.nu_1 = f_viscosity(PHYS.alpha, self.c_s_1, self.Omega_1)

        self.t_s = (1./ (3 * (2 - PHYS.gamma)**2)) * PHYS.R_1**2 / self.nu_1

        self.scaling_C = (3./2) * self.nu_1 * PHYS.M_disk_0 * (2-PHYS.gamma) / PHYS.R_1**2

        self.Omega = f_Omega(R, PHYS.M_star)
        self.temp = f_temperature(R, PHYS.T_star, PHYS.q, PHYS.irr_angle)
        self.c_s = f_sound_speed(self.temp)
        self.v_k = R * self.Omega
        self.H = self.c_s/self.Omega

        self.m_dst = (4/3)*pi * PHYS.dst_dens * (PHYS.dst_size/2)**3

        self.m_pls =(4/3)*pi * PHYS.pls_dens * PHYS.r_pls**3



class local():

    def __init__(self, R, PHYS, NUM):

        """
        Get nonevolving quantitities
        """

        self.CONST = var_nonevolving(R, PHYS)

        """
        Get initial conditions
        """

        self.t_init = float(0)

        # ICELINES TBA

        self.dst_tgr_init = PHYS.metallicity # initial dust-to-gas ratio
        self.Sigma_g_init =  f_similarity_solution(R/PHYS.R_1, 0, self.CONST.scaling_C, self.CONST.nu_1, PHYS.gamma, self.CONST.t_s, PHYS.l_gasevolution)
        self.Sigma_dst_init = self.dst_tgr_init * self.Sigma_g_init
        self.Sigma_pbb_init = float(0)
        self.Sigma_pls_init = float(0)
        self.Sigma_debris_init = float(0)
        self.Sigma_new_gen_pls_init = float(0)
        
    
    def initialize_arrays(self):
        dt_init = float("NaN")

        self.VAR = var(t = np.array([self.t_init]),
                        dt = np.array([dt_init]),
                        sig_g = np.array([self.Sigma_g_init]),
                        sig_dst = np.array([self.Sigma_dst_init]),
                        sig_pbb = np.array([self.Sigma_pbb_init]),
                        sig_pls_first = np.array([self.Sigma_pls_init]),
                        sig_pls_second = np.array([self.Sigma_new_gen_pls_init]),
                        sig_dbr_um = np.array([self.Sigma_debris_init]),
                        sig_dbr_mm = np.array([self.Sigma_debris_init]),
                        sig_dbr_dm = np.array([self.Sigma_debris_init]),
                        sig_dbr_m = np.array([self.Sigma_debris_init]),
                        sig_drift_loss = np.array([float(0)])
                        )
        
        # inititalize data
        self.DATA = data(self.VAR)
        
        

    def integrate(self, R, PHYS, NUM):
        
        t = self.t_init
        nt = int(0) # counts timesteps

        while t < NUM.t_max:
            """
                Break if too many max time step is reached
            """
            if nt >= NUM.nt_max:
                break


            """
                Dust growth
            """
            
            DUST_GROWTH = dust_growth(R, self.VAR, self.CONST, PHYS, t)

            dust_growth_flux = DUST_GROWTH.dust_growth_flux

            """
               Planetesimal formation
            """
            
            FORMATION = formation(R, self.VAR, self.CONST, PHYS, DUST_GROWTH)

            pls_form_flux_first = FORMATION.pls_form_flux_first
            pls_form_flux_second = FORMATION.pls_form_flux_second
            pls_form_flux_dm = FORMATION.pls_form_flux_dm
            pls_form_flux_m = FORMATION.pls_form_flux_m

            """
               Planetesimal Collisions
            """

            COLLISIONS = collisions(R, self.VAR, self.CONST, PHYS)

            collision_flux_pls_first = COLLISIONS.collision_flux_pls_first

            collision_flux_um = COLLISIONS.collision_flux_um
            collision_flux_mm = COLLISIONS.collision_flux_mm
            collision_flux_dm = COLLISIONS.collision_flux_dm
            collision_flux_m = COLLISIONS.collision_flux_m

            collision_flux_pls_second  = COLLISIONS.collision_flux_pls_second

            """
                Drift
            """
            DRIFT = drift(R, self.VAR, self.CONST, PHYS, DUST_GROWTH)

            drift_flux_pbb = DRIFT.drift_flux_pbb
            drift_flux_dbr_mm = DRIFT.drift_flux_dbr_mm
            drift_flux_dbr_dm = DRIFT.drift_flux_dbr_dm
            drift_flux_dbr_m = DRIFT.drift_flux_dbr_m

            """
                Get surface density fluxes
            """

            sig_dst_flux = - dust_growth_flux
            sig_pbb_flux = + dust_growth_flux - pls_form_flux_first - drift_flux_pbb
            sig_pls_first_flux = pls_form_flux_first + collision_flux_pls_first
    
    
            sig_dbr_um_flux = + collision_flux_um
            sig_dbr_mm_flux = + collision_flux_mm - drift_flux_dbr_mm
            sig_dbr_dm_flux = + collision_flux_dm - pls_form_flux_dm - drift_flux_dbr_dm
            sig_dbr_m_flux = + collision_flux_m - pls_form_flux_m - drift_flux_dbr_m
    
            sig_pls_second_flux = pls_form_flux_second + collision_flux_pls_second

            sig_drift_loss_flux = + drift_flux_pbb + drift_flux_dbr_mm + drift_flux_dbr_dm + drift_flux_dbr_m

            """
                Calculate timestep
            """

            # write in arrays
            sig_flux_pairs = [(self.VAR.sig_dst, sig_dst_flux),
                                (self.VAR.sig_pbb, sig_pbb_flux),
                                (self.VAR.sig_pls_first, sig_pls_first_flux),
                                (self.VAR.sig_pls_second, sig_pls_second_flux),
                                (self.VAR.sig_dbr_um, sig_dbr_um_flux),
                                (self.VAR.sig_dbr_mm, sig_dbr_mm_flux),
                                (self.VAR.sig_dbr_dm, sig_dbr_dm_flux),
                                (self.VAR.sig_dbr_m, sig_dbr_m_flux)]


            dt_old = NUM.t_max # initialize max time step for comparison
            for sig, flux in sig_flux_pairs:
                if np.abs(flux) != 0:
                    dt_current = NUM.C_max * sig / np.abs(flux)
                else:
                    dt_current = dt_old
                if dt_current < dt_old:
                    dt_old = dt_current

                
            dt = np.amax([NUM.dt_min, dt_old])
            

            """
                Make timestep
            """


            sig_dst = self.VAR.sig_dst + sig_dst_flux * dt

            sig_pbb = self.VAR.sig_pbb + sig_pbb_flux * dt

            sig_pls_first = self.VAR.sig_pls_first + sig_pls_first_flux * dt

            sig_pls_second = self.VAR.sig_pls_second + sig_pls_second_flux * dt

            sig_dbr_um = self.VAR.sig_dbr_um + sig_dbr_um_flux * dt
            sig_dbr_mm = self.VAR.sig_dbr_mm + sig_dbr_mm_flux * dt
            sig_dbr_dm = self.VAR.sig_dbr_dm + sig_dbr_dm_flux * dt
            sig_dbr_m = self.VAR.sig_dbr_m + sig_dbr_m_flux * dt

            sig_drift_loss = self.VAR.sig_drift_loss + sig_drift_loss_flux * dt

            """
                Check if a surface density falls below sig_min
                Set to zero if it does
                This prevents small time steps when surface densities get very small
            """
            
            if sig_dst < NUM.sig_min and sig_dst_flux < 0:
                sig_dst = float(0)
            if sig_pbb < NUM.sig_min and sig_pbb_flux < 0:
                sig_pbb = float(0)
            if sig_pls_first < NUM.sig_min and sig_pls_first_flux < 0:
                sig_pls_first = float(0)
            if sig_pls_second < NUM.sig_min and sig_pls_second_flux < 0:
                sig_pls_second = float(0)
            if sig_dbr_um < NUM.sig_min and sig_dbr_um_flux < 0:
                sig_dbr_um = float(0)
            if sig_dbr_mm < NUM.sig_min and sig_dbr_mm_flux < 0:
                sig_dbr_mm = float(0)
            if sig_dbr_dm < NUM.sig_min and sig_dbr_dm_flux < 0:
                sig_dbr_dm = float(0)
            if sig_dbr_m < NUM.sig_min and sig_dbr_m_flux < 0:
                sig_dbr_m = float(0)
    	    

            """
                Calculate gas surface density at new timestep
            """
            GAS = gas_evolution(R, t, dt, self.CONST, PHYS)
            sig_g = GAS.sig_g
            

            """
                Renew VAR
            """
            t = t + dt

            self.VAR.renew(t = t,
                                dt = dt,
                                sig_g = sig_g,
                                sig_dst = sig_dst,
                                sig_pbb = sig_pbb,
                                sig_pls_first = sig_pls_first,
                                sig_pls_second = sig_pls_second,
                                sig_dbr_um = sig_dbr_um,
                                sig_dbr_mm = sig_dbr_mm,
                                sig_dbr_dm = sig_dbr_dm,
                                sig_dbr_m = sig_dbr_m,
                                sig_drift_loss = sig_drift_loss
                        )

            # Append data
            if nt%NUM.save_nt == 0:
                self.DATA.append(self.VAR)

            
            if NUM.l_print_progress_local == True:
                if nt%1000 == 0:
                    print(str(nt) + ' nt | ' + str(int(self.VAR.t/year)) + ' yr | ' + str(np.round(100*self.VAR.t/NUM.t_max, 2)) + ' % |')
            nt += 1


    def plot_local(self, R, PHYS, NUM, plotpath):
        title = str(np.round(R/AU, 2))
        title = title.replace('.', '_')    

        plt.figure()
        plt.plot(self.DATA.t/year, self.DATA.dt/year)
        plt.xscale('log')
        plt.savefig(plotpath + '/local_' + title + '_a.png')
        plt.close()

        plt.figure()
        plt.plot(self.DATA.t/year, self.DATA.sig_dst)
        plt.plot(self.DATA.t/year, self.DATA.sig_pbb)
        plt.plot(self.DATA.t/year, self.DATA.sig_pls_first)
        plt.plot(self.DATA.t/year, self.DATA.sig_pls_second)
        plt.plot(self.DATA.t/year, self.DATA.sig_dbr_um)
        plt.plot(self.DATA.t/year, self.DATA.sig_dbr_mm)
        plt.plot(self.DATA.t/year, self.DATA.sig_dbr_dm)
        plt.plot(self.DATA.t/year, self.DATA.sig_dbr_m)

        plt.plot(self.DATA.t/year, self.DATA.sig_g)


        plt.ylim(2*NUM.sig_min, np.amax(self.DATA.sig_g))
        plt.xlim(10**3, NUM.t_max/year)

        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(plotpath + '/local_' + title + '_b.png')
        plt.close()
        return None