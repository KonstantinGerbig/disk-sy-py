import numpy as np
from constants import *
from params import params_phys, params_num
from local import *
from physics import *
import pickle
import os









class disk():
    def __init__(self, diskpath, plotpath):
        self.path = diskpath
        self.plotpath = plotpath
        return None
    
    def load_local_data(self, R):
        # pkl filenaming scheme must match the one in local.py:data.write_to_file
        filename = str(np.round(R/AU, 2))
        filename = filename.replace('.', '_') + '.pkl'
        try:
            data = pickle.load(open(self.path + '/local/' + filename, 'rb'))
            setattr(self, str(R/AU), data)
        except:
            import warnings
            warnings.warn("Couldn't load local data of R = " + str(R/AU) + " ... Check if pkl file exists!")
        return None

    def evaluation(self, Rs, NUM, PHYS, N_t_grid = 100, l_profiles = True, l_total = True, l_calc_mm_flux = True):
        # check if at least one of the evaluation flags is on
        if (l_profiles != True) and (l_total != True):
            import warnings
            warnings.warn("Disk evaluation called but don't know what to evaluate: No flag is set to true")
            return None


        # check if local data is loaded
        local_data_missing = False
        for R in Rs:
            if not hasattr(self, str(R/AU)):
                local_data_missing = True
        if local_data_missing == True:
            import warnings
            warnings.warn("Couldn't evaluate disk: Local data missing. Check if if load_local_data() was invoked")
            return None
        else:
            print('Evaluating disk ...')

        
        # create directory for radial profiles
        if l_profiles == True:
            if not os.path.exists(self.path + '/profiles'):
                os.mkdir(self.path + '/profiles')
        # create directory for total quantities
        if l_total == True:
            if not os.path.exists(self.path + '/total'):
                os.mkdir(self.path + '/total')

        # The below array is the temporal grid for evaluations.
        # This includes: 
        #  - the radial profile for different t for different species
        #  - Total mass evolution for different species
        # This is necessary since the dynamic time step leads to different t grids at each R
        
        t_grid = np.logspace(np.log10(1000*year), np.log10(NUM.t_max), N_t_grid)

        # save t grid to directories:
        if l_profiles == True:
            pickle.dump(t_grid, open(self.path + '/profiles/t_grid.pkl', 'wb'))
        if l_total == True:
            pickle.dump(t_grid, open(self.path + '/total/t_grid.pkl', 'wb'))


        # These are the species that are being evaluated.
        # The sig keys need to match the names of the attributes of local.var. order does not matter
        key_tuples = [('M_g', 'sig_g'),
            ('M_dst', 'sig_dst'),
            ('M_pbb', 'sig_pbb'),
            ('M_pls_first', 'sig_pls_first'),
            ('M_pls_second', 'sig_pls_second'),
            ('M_dbr_um', 'sig_dbr_um'),
            ('M_dbr_mm', 'sig_dbr_mm'),
            ('M_dbr_dm', 'sig_dbr_dm'),
            ('M_dbr_m', 'sig_dbr_m'),
            ('M_drift_loss', 'sig_drift_loss')]

        # initialize M dict
        if l_total == True:
            M_dict = {key:values for key, values in zip([k[0] for k in key_tuples], [np.array([]) for k in key_tuples])}

        # go through each timestep in temporal grid

        i_t_grid = 0 # iteration counter. Needed for filenaming scheme
        for t_approx in t_grid:

            # Initialize sig profile dictonary 
            # (gets overwritten each t_approx in t_grid so it is written to pkl file)
            sig_dict = {key:values for key, values in zip([k[1] for k in key_tuples], [np.array([]) for k in key_tuples])}

            # go through each R:
            for R in Rs:
                R_list_object = getattr(self, str(R/AU))

                # find t array
                for data_tuple in R_list_object[:]:
                    if data_tuple[0] == 't':
                        t_current_array = data_tuple[1]

                # find t in t_current_array closest to t_approx
                def find_nearest(array, value):
                    array = np.asarray(array)
                    idx = (np.abs(array - value)).argmin()
                    # array[idx] returns the nearest value
                    return idx
                
                i_approx = find_nearest(t_current_array, t_approx)

                # go through each species and get sig arrays
                for sig_key in sig_dict.keys():
                    # find correct sig array
                    for data_tuple in R_list_object[:]:
                        if data_tuple[0] == sig_key:
                            sig_dict[sig_key] = np.append(sig_dict[sig_key], data_tuple[1][i_approx])
            
            # save radial profiles (directory was created above)
            if l_profiles == True:
                pickle.dump(sig_dict, open(self.path + '/profiles/sig_dict_'+ str(i_t_grid) + '.pkl', 'wb'))

            """
            Calculate radial flux profile
            """

            if l_calc_mm_flux == True:

                opacity_kappa = 3.07107

                # this needs to be more sophisticated
                sig_mm_profile = sig_dict['sig_pbb'] + sig_dict['sig_dbr_mm']

                # optical_depth = opacity*sigma*cos(disk inclination)
                
                

                nu_mm = 0.1/c
                
                Flux_profile = np.array([])
                optical_depth_profile = np.array([])
                j = 0
                for R in Rs:
                    Omega = f_Omega(R, PHYS.M_star)

                    T = f_temperature(R, PHYS.T_star, PHYS.q, PHYS.irr_angle)

                    disk_inclination = np.arctan(f_sound_speed(T)/(Omega*R))

                    optical_depth = opacity_kappa * sig_mm_profile[j] * np.cos(disk_inclination)
                    optical_depth_profile = np.append(optical_depth_profile, optical_depth)

                    #Flux = f_Planck_function(T, nu_mm) * (1- np.exp(-optical_depth))* Omega

                    #Flux_profile = np.append(Flux_profile, Flux)

                    j += 1

                #if i_t_grid == len(t_grid)-1:
                    #print(sig_mm_profile)
                    #print(optical_depth_profile)
                    #print(Flux_profile)

                # save optical depth profiles
                pickle.dump(optical_depth_profile, open(self.path + '/profiles/tau_mm_'+ str(i_t_grid) + '.pkl', 'wb'))

            # Calc total quantities
            if l_total == True:

                # each species
                for sig_key in sig_dict.keys():
                    sig_current = sig_dict[sig_key]

                    # integrate
                    M_current = 0
                    for i_R in range(len(Rs)-1):
                        M_current = M_current + (Rs[i_R+1]-Rs[i_R])*(Rs[i_R]*sig_current[i_R] + Rs[i_R+1]*sig_current[i_R+1])
                    M_current = M_current * np.pi

                    # write M_current to M_dict
                    for key_tuple in key_tuples:
                        if key_tuple[1] == sig_key:
                            M_key = key_tuple[0]
                            M_dict[M_key] = np.append(M_dict[M_key], M_current)

            # increase iteration counter
            i_t_grid += 1

        # save M_dict (the directory is created above)
        if l_total == True:
            pickle.dump(M_dict, open(self.path + '/total/M_dict.pkl', 'wb'))

        return None






    def plot_M(self, PHYS):
        try:
            M_dict = pickle.load(open(self.path + '/total/M_dict.pkl', 'rb'))
            t_grid = pickle.load(open(self.path + '/total/t_grid.pkl', 'rb'))
        except:
            import warnings
            warnings.warn("Couldn't find M_dict file. Did you call disk.evaluation with l_total = True? Stopping plotting of mass evolution.")
            return None

        print('Plotting mass evolution')

        plt.figure()
        for key in M_dict.keys():
            key_label = key.replace('_', ' ')
            plt.plot(t_grid/year, M_dict[key]/PHYS.M_star, label = key_label)

        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc = 'best')
        plt.xlim(10**3, t_grid[-1]/year)
        plt.ylim(10**(-5), 10**(-3))

        plt.xlabel(r'$t$ [\rm{yr}]')
        plt.ylabel(r'$M \, / M_{\rm{star}}$')

        plt.savefig(self.plotpath + '/total_M.png')
        plt.close()
        
        return None









    def plot_profiles(self, PHYS, Rs, which_ts = 'last'):

        try:
            t_grid = pickle.load(open(self.path + '/profiles/t_grid.pkl', 'rb'))
        except:
            import warnings
            warnings.warn("Couldn't find t_grid file. Did you call disk.evaluaton with l_profiles = true")
            return None

        print('Plotting profiles ...')

        i_t = 0 # iteration counter. Needed for filenaming scheme
        for t in t_grid:
            
            # Initialize bool plot
            bool_plot_t = False

            if which_ts == 'last' or which_ts == 'Three':
                if t == t_grid[-1]:

                    bool_plot_t = True
            if which_ts == 'mid' or which_ts == 'Three':
                if i_t == int(len(t_grid)/2):
                    bool_plot_t = True

            if which_ts == 'first' or which_ts == 'Three':
                if i_t == 0:
                    bool_plot_t = True

            if which_ts == 'all':
                bool_plot_t = True


            if bool_plot_t == True:
                
                try:
                    sig_dict = pickle.load(open(self.path + '/profiles/sig_dict_'+ str(i_t) + '.pkl', 'rb'))
                except:
                    import warnings
                    warnings.warn("Couldn't find sig_dict file of  = " + str(i_t) + " AU. Did you call disk.evaluation with l_profiles = True?")
                    break

                plt.figure()
                for key in sig_dict.keys():
                    plt.plot(Rs/AU, sig_dict[key])

                plt.xlim(Rs[0]/AU, Rs[-1]/AU)
                plt.ylim(10**(-4), 10**3)

                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel(r'$R$ [AU]')
                plt.ylabel(r'$\Sigma \, [\rm{g}/\rm{cm}^2]$')
                plt.title("{:.2e}".format(np.round(t/year, 0)))
                plt.savefig(self.plotpath + '/profile_' + str(i_t) + '.png')
                plt.close()

            i_t += 1
        return None
    











    def print_attributes(self):
        import inspect
        attributes_all = inspect.getmembers(self, lambda x:not(inspect.isroutine(x)))
        attributes = [x for x in attributes_all if not(x[0].startswith('__') and x[0].endswith('__'))]
        print(attributes)
        return None













def disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS):
    """
    This is a wrapper for running the integrator at a range of semi-major axises in order to simulate an entire disk
    Arguments:
    ----------
    ARGS : instance of the input parameter object
    Keywords:
    ---------
    save : bool
          whether or not to write the data to disk
    Output:
    -------
    results : instance of the results object
    """

    Rs = PARAMS_POP.Rs

    """
    Make disk data directory
    """

    diskname = NUM.diskname
    print('')
    print('Simulating disk: ' + diskname)

    diskpath =  datadir + '/' + diskname
    if not os.path.exists(diskpath):
        os.mkdir(diskpath)
        if os.path.exists(diskpath):
            print('Created disk data directory')
    else:
        print('Disk data directory already exists. Nothing to create.')

    if not os.path.exists(diskpath + '/local'):
        os.mkdir(diskpath + '/local')
        if os.path.exists(diskpath + '/local'):
            print('Created local data directory')
    else:
        print('Local data directory already exists. Nothing to create.')


    if (NUM.l_plotlocal == True) or (PARAMS_POP.l_plotmass == True) or (PARAMS_POP.l_plot_profiles == True): 
        plotpath =  diskpath + '/' + NUM.plotdir
        # check if plot directory exist
        if not os.path.exists(plotpath):
            os.mkdir(plotpath)
            print('Created directory for disk plots')
        else:
            print('Plot directory already exists')
    else:
        plotpath = None


    if PARAMS_POP.l_run == True:
        print('Running Disk')

        i_R = 0
        for R in Rs:
            if PARAMS_POP.l_print_progress_disk == True:
                str(int((i_R+1)/len(Rs)*100))
                print('R = ' + str(np.round(R/AU, 3)) + ' AU | ' + str(int((i_R+1)/len(Rs)*100)) + ' %')
            LOCAL = local(R, PHYS, NUM)
            LOCAL.initialize_arrays()
            LOCAL.integrate(R, PHYS, NUM)
            LOCAL.DATA.write_to_file(R, diskpath)

            if NUM.l_plotlocal == True:
                print('Plotting local ...')
                LOCAL.plot_local(R, PHYS, NUM, plotpath)

            i_R += 1
                
    else:
        print('l_run = False: Not running disk model')
    

    # Inititialize disk object
    DISK = disk(diskpath, plotpath)

    for R in Rs:
        DISK.load_local_data(R)
    
    #DISK.print_attributes()

    # write diskvar
    DISK.evaluation(Rs, NUM, PHYS)

    if PARAMS_POP.l_plotmass == True:
        DISK.plot_M(PHYS)

    if PARAMS_POP.l_plot_profiles == True:
        DISK.plot_profiles(PHYS, Rs, which_ts = 'all')

    if PARAMS_POP.l_light_save == True:
        import shutil
        try:
            shutil.rmtree(diskpath + '/local')
            shutil.rmtree(diskpath + '/profiles')
        except:
            import warnings
            warnings.warn('Could not delete unnessary files')
    return None
