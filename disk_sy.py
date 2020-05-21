import numpy as np
from params import params_population
import os
from disk_wrapper import *
from constants import *


def population_wrapper():
    PARAMS_POP = params_population()

    # check if datadir exist, create if it doesn't exist
    datadir = os.getcwd() + '/' + PARAMS_POP.datadir
    if not os.path.exists(datadir):
        os.mkdir(PARAMS_POP.datadir)
        if os.path.exists(datadir):
            print('Created main data directory')
    else:
        print('Main data directory already exists. Nothing to create.')


    
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


    return None

population_wrapper()