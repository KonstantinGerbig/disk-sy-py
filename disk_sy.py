import numpy as np
from params import params_population
import os
from disk_wrapper import *
from constants import *
import pickle

class disk_population():

    def __init__(self):
        return None

    def make_datadir(self, PARAMS_POP):
        # check if datadir exist, create if it doesn't exist
        datadir = os.getcwd() + '/' + self.PARAMS_POP.datadir
        if not os.path.exists(datadir):
            os.mkdir(self.PARAMS_POP.datadir)
            if os.path.exists(datadir):
                print('Created main data directory')
        else:
            print('Main data directory already exists. Nothing to create.')
        return datadir

    def print_params(self, path, NUM, PHYS):

        # write all attributes of NUM and PHYS in dictonaries
        import inspect

        phys_dict = {}
        for i in inspect.getmembers(PHYS):
            if not i[0].startswith('_'):
                phys_dict.update({i[0] : i[1]})

        num_dict = {}
        for i in inspect.getmembers(NUM):
            if not i[0].startswith('_'):
                num_dict.update({i[0] : i[1]})

        
        # Pickle dump the dictionaries

        pickle.dump(num_dict, open(path + '/num.pkl', 'wb'))
        pickle.dump(phys_dict, open(path + '/phys.pkl', 'wb'))
        return 0

    def high_res_default_disk(self):

        Rs = np.logspace(np.log10(1*AU), np.log10(100*AU), 30)
        self.PARAMS_POP = params_population(Rs = Rs, l_light_save = False, l_plotmass = True, l_plot_profiles = True)
        datadir = self.make_datadir(self.PARAMS_POP)

        diskname = 'default_high_res'
        NUM = params_num(diskname = diskname)
        PHYS = params_phys()
        disk_wrapper_function(datadir, self.PARAMS_POP, NUM, PHYS)

        return None

    
    def sample0(self):
        # generates a grid based sample population.
        # Numerical properties are kept constant

        self.PARAMS_POP = params_population()
        datadir = self.make_datadir(self.PARAMS_POP)

        trapping_efficiency_array = np.array([1-3e-3, 1-3e-2, 0.85])
        frag_velocity_array = np.array([1*100, 3*100, 10*100])
        alpha_array = np.array([1e-4, 1e-3, 1e-2])
        conversion_efficiency_array = np.array([0.1, 0.5, 0.8])
        M_disk_0_array = np.array([0.01, 0.02, 0.03])
        R_1_array = np.array([30*AU, 50*AU])
        M_star_array = np.array([1*M_sun, 3*M_sun])

        # add something that calculates how many disks are calculated

        i = 0
        for a in trapping_efficiency_array:
            for b in frag_velocity_array:
                for c in alpha_array:
                    for d in conversion_efficiency_array:
                        for e in M_disk_0_array:
                            for f in R_1_array:
                                for g in M_star_array:
                                    diskname = 'sample0_disk_' + str(i)
                                    NUM = params_num(diskname = diskname)
                                
                                    PHYS = params_phys(trapping_efficiency = a,
                                        frag_velocity = b,
                                        alpha = c,
                                        conversion_efficiency = d,
                                        M_disk_0 = e,
                                        R_1 = f,
                                        M_star = g)
                                    disk_wrapper_function(datadir, self.PARAMS_POP, NUM, PHYS)

                                    self.print_params(datadir + '/' + diskname, NUM, PHYS)
                                    i += 1
        
        return None


def population_wrapper():
    #PARAMS_POP = params_population()

    
    POPULATION = disk_population()

    #POPULATION.sample0()

    POPULATION.high_res_default_disk()



    """
    NUM = params_num(diskname = 'finetune')
    PHYS = params_phys(
        trapping_efficiency = 1-3e-3,
        frag_velocity = 10*100,
        alpha = 2e-3,
        R_1 = 30*AU,
        M_disk_0 = 0.025*M_sun,
        conversion_efficiency = 0.2
        )
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    """



    return None

population_wrapper()