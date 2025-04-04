"""
Script to create SOAP-like hdf5 files containing
halo luminosities integrated over filters.
Created by Anna Durrant on 04.4.2025.
"""

import numpy as np
from datetime import datetime
from swiftsimio.visualisation.smoothing_length.generate import generate_smoothing_lengths as gsl
import unyt
import yaml
import argparse
import os
import h5py
from swiftsimio.objects import cosmo_array, cosmo_factor, a

# Set simName
parser = argparse.ArgumentParser(
    description="Script to create particle .txt files, based on Kyle Oman's swiftgalaxy framework."
)

parser.add_argument(
    "BoxSize",
    type=int,
    help="Boxsize of the simulation in Mpc.",
)

parser.add_argument(
    "Resolution",
    type=int,
    help="Particle mass resolution of the simulation in log10(M/Msun).",
)

parser.add_argument(
    "--snaps",
    type=int,
    required=True,
    nargs='+',
    help="<Required> Snapshot number(s).",
)

parser.add_argument(
    "--mode",
    type=str,
    default="Thermal", # Thermal AGN feedback with non-equilibrium chemistry
    help="Simulation mode (default: Thermal).",
)

args = parser.parse_args()

sim = 'L{:03.0f}_m{:01.0f}'.format(args.BoxSize, args.Resolution)
simName = sim + '/' + args.mode

# Define filepaths from parameter file
dir_path = os.path.dirname(os.path.realpath(__file__))
with open(f'{dir_path}/SKIRT_parameters.yml','r') as stream:
    params = yaml.safe_load(stream)

simPath = params['InputFilepaths']['simPath'].format(simName=simName)
sampleFolder = params['OutputFilepaths']['sampleFolder'].format(simPath=simPath)

# Read in filters 
filterNames = [f for f in os.listdir(f'{dir_path}/Filters')]

luminosityFilters = {}

for filter in filterNames:
    luminosityFilters[filter.split('.')[0]] = np.genfromtxt(f'{dir_path}/Filters/{filter}')

def integrate_SED(
        SED,
        filters,
):
    integrated_vals = np.empty_like(filters,dtype=float)
    
    for idx, filter in enumerate(filters):
        # convolve SED with filter
        integrated_vals[idx] = SED * filter
        
    return integrated_vals

for snap in args.snaps:

    halo_IDs, Rstar = np.loadtxt(sampleFolder + '/sample_' + str(snap) + '.txt', unpack = True, usecols = [0, 2])
    halo_IDs = halo_IDs.astype(int)

    halo_luminosities_array = np.zeros( (np.max(halo_IDs), len(luminosityFilters)) )

    for ID in halo_IDs:

        halo_SED = params['OutputFilepaths']['haloSEDfiles'].format(simPath=simPath,snap_nr=snap)

        halo_luminosities_array[ID] = integrate_SED(halo_SED, luminosityFilters)


    # Open / create hdf5 file
    haloPropertiesFilepath = params['OutputFilepaths']['haloPropertiesFile'].format(simPath=simPath,snap_nr=snap)
    haloPropertiesFile = h5py.File(haloPropertiesFilepath,'a')

    try:
        SKIRTluminosities_dset = haloPropertiesFile['ExclusiveSphere/100kpc/Luminosities']
    except:
        # SKIRTluminosities dataset does not exist yet
        SKIRTluminosities_dset = haloPropertiesFile.create_dataset('ExclusiveSphere/100kpc/Luminosities')
        





