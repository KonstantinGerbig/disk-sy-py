import numpy as np
from params import params_population
import os
from disk_wrapper import *
from constants import *

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

    def print_params(self, datadir, NUM, PHYS):
        
        return 0
    
    def sample0(self):
        # generates a grid based sample population.
        # Numerical properties are kept constant

        self.PARAMS_POP = params_population()
        datadir = self.make_datadir(self.PARAMS_POP)

        trapping_efficiency_array = np.array([1-3e-3, 1-3e-2])
        frag_velocity_array = np.array([3*100])
        alpha_array = np.array([2e-3])
        conversion_efficiency_array = np.array([0.1])

        i = 0
        for a in trapping_efficiency_array:
            for b in frag_velocity_array:
                for c in alpha_array:
                    for d in conversion_efficiency_array:
                        diskname = 'sample0_disk_' + str(i)
                        NUM = params_num(diskname = diskname)
                        PHYS = params_phys(trapping_efficiency = a, frag_velocity = b, alpha = c, conversion_efficiency = d)

                        self.print_params(self, datadir, NUM, PHYS)
                        disk_wrapper_function(datadir, self.PARAMS_POP, NUM, PHYS)
                        i += 1
        
        return None


def population_wrapper():
    #PARAMS_POP = params_population()

    
    POPULATION = disk_population()

    POPULATION.sample0()

    """
    
    NUM = params_num(diskname = 'disk1')
    PHYS = params_phys()
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'disk2')
    PHYS = params_phys(alpha = 1e-4)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'disk3')
    PHYS = params_phys(R_1 = 60*AU, M_disk_0 = 0.1*M_sun)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'disk4')
    PHYS = params_phys(conversion_efficiency = 0.01)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'disk5')
    PHYS = params_phys(trapping_efficiency = 0.99)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'disk6')
    PHYS = params_phys(frag_velocity = 10*100)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'disk7')
    PHYS = params_phys(frag_velocity = 1*100)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)
    

    NUM = params_num(diskname = 'visc1')
    PHYS = params_phys(alpha = 1e-2)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'visc2')
    PHYS = params_phys(alpha = 5e-3)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'visc3')
    PHYS = params_phys(alpha = 5e-4)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'visc4')
    PHYS = params_phys(alpha = 1e-4)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    

    NUM = params_num(diskname = 'visc5')
    PHYS = params_phys(alpha = 2e-3)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'visc6')
    PHYS = params_phys(alpha = 3e-3)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'visc7')
    PHYS = params_phys(alpha = 4e-3)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'drift1')
    PHYS = params_phys(trapping_efficiency = 1-0)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'drift2')
    PHYS = params_phys(trapping_efficiency = 1-1e-4)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'drift3')
    PHYS = params_phys(trapping_efficiency = 1-1e-3)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'drift4')
    PHYS = params_phys(trapping_efficiency = 1-5e-3)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    NUM = params_num(diskname = 'drift5')
    PHYS = params_phys(trapping_efficiency = 1-1e-2)
    disk_wrapper_function(datadir, PARAMS_POP, NUM, PHYS)

    

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