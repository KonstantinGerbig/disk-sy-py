from constants import *
#from physics import *

class params_phys():

    """
    Mass
    """
    M_star = 1* M_sun
    M_disk_0 = 0.02 * M_star

    """
    Temperature profile
    """
    T_star = 4000
    R_star = 1.25 * R_sun
    irr_angle = 0.1
    q = 0.5 # temperature profile exponent

    """
    gas evolution: viscous time and viscosity
    """
    l_gasevolution = True

    alpha = 1e-3
    gamma = 1 # radial dependance of viscosity, maybe couple via: gamma = 3/2 - q
    R_1 = 30 * AU # scaling radius

    """
    Icelines and dust composition
    """
    l_icelines = True

    temp_freeze_H2O = 140 # K Pinilla, Pohl 2017
    temp_freeze_NH3 = 80 # K
    temp_freeze_CO = 20 # K

    # ricci testi natta 2010. fig 3
    H20_fraction = 0.3
    CO_fraction = 0.2
    NH3_fraction = 0.01 # questionable


    """
    Dust properties
    """
    metallicity = 0.01
    dst_size = 10**(-4) # cm (1 micro meter)
    dst_dens = 1.2 # g/cm^3

    l_dustgrowth = True

    """
    Dust drift
    """
    l_drift = True

    trapping_efficiency = 0.999

    """
    Dynamical pebble stokes number
    """
    frag_velocity = 5*100 # cm/s

    """
    Planetesimal formation
    """
    l_formation = True
    l_formation_second_generation = True

    conversion_efficiency = 0.1
    # formation efficiency = trapping efficiency * conversion efficiency

    trapping_distance = 5 # in units of H

    """
    Planetesimal collisions
    """

    l_collisions = True

    # for timescale

    r_pls = 50 * 10**5 # cm 
    pls_dens = dst_dens

    # Fragment distribution
    xi = 1.83 # fragment_power_law

    a_0 = dst_size # 1 micro m 
    a_1 = 100 * a_0 # 0.1* mm
    a_2 = 100 * a_1 # 1 cm 
    a_3 = 100 * a_2 # 1 m
    a_4 = 100 * a_3 # 100 m 
    a_5 = 2*r_pls # 100 km

    def __init__(self, **kwargs):
        """
        Overwrite default params
        """
        import warnings
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                warnings.warn('Argument not found')
    

class params_num():
    """
    Numerical parameters for each disk
    """

    t_max = 10**8 * year # run_time

    C_max = 0.01 # "Courant" number

    dt_min = 10 * year # minimum time step
    #nt_max = 10**5 # max number of time steps

    nt_max = 10000

    save_nt = 5 # how often is VAR saved

    l_print_progress_local = False 

    integration_scheme = 'Euler'

    l_plotlocal = False # make local plots

    plotdir = 'plots'

    diskname = 'disk'


    # If a surface density falls below this value, it is set to zero
    # prevents time step issues
    sig_min = 10**(-10) # in g/cm^2.


    def __init__(self, **kwargs):
        """
        Overwrite default params with kwargs
        """
        import warnings
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                warnings.warn('Argument not found')

class params_population():

    l_run = True

    datadir = '../data'

    #Rs = np.array([0.1*AU, 1*AU, 3*AU, 5*AU, 7*AU, 10*AU, 20*AU, 30*AU, 40*AU, 60*AU, 100*AU])
    #Rs = np.linspace(1,100, 50)*AU
    Rs = np.logspace(np.log10(1*AU), np.log10(100*AU), 20)


    l_plotmass = False

    l_plot_profiles = False

    l_print_progress_disk = True

    def __init__(self, **kwargs):
        """
        Overwrite default params with kwargs
        """
        import warnings
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                warnings.warn('Argument not found')
