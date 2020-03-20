import numpy as np
from params import params_population
import os
from disk_wrapper import *


def population_wrapper():
    PARAMS = params_population()

    # check if datadir exist, create if it doesn't exist
    datadir = os.getcwd() + '/' + PARAMS.datadir
    if not os.path.exists(datadir):
        os.mkdir(PARAMS.datadir)
        if os.path.exists(datadir):
            print('Created main data directory')
    else:
        print('Main data directory already exists. Nothing to create.')

    disk_wrapper_function(datadir, PARAMS)

    return None

population_wrapper()