"""
Edit a ski file
"""

import numpy as np
import subprocess
import sys
import warnings
from datetime import datetime
import unyt
import os
import yaml

startTime = datetime.now()

# Global settings

SKIRTboxsize = unyt.unyt_quantity(100., 'kpc')
old_stars_tmin = unyt.unyt_quantity(10., 'Myr') # Minimum age in Myr for an evolved star particle. Also determines the TODDLERS averaging timescale
# Don't change this unless you know what you're doing :)

Npp = int(10**7.5) # Number of photon packets
binTreeMaxLevel = 36 # Max refinement level of the spatial grid

snapNum = sys.argv[1]
haloID = sys.argv[2]
Rstar = float(sys.argv[3])

txtFilePath = sys.argv[4]
SKIRTinputFilePath = sys.argv[5]

skifileversion = '4.0'

mediumXMLlineDict = {'4.0': (40, 197)} # Lines in the .ski file that describe the medium system

# Define filepaths from parameter file
dir_path = os.path.dirname(os.path.realpath(__file__))
with open(f'{dir_path}/SKIRT_parameters.yml','r') as stream:
    params = yaml.safe_load(stream)

# Edit ski file

def editSki(snapNum, haloID, Rstar):

    SKIRTinputFiles = SKIRTinputFilePath + 'snap' + snapNum + '_ID' + haloID

    skifilename = params['InutFilepaths']['skiFilepath'].format(skifileversion=skifileversion)

    skifilename_halo = 'snap' + snapNum + '_ID' + haloID + '.ski'


    subprocess.run(['cp', skifilename, skifilename_halo]) # copy the skirt file for each galaxy

    # Calculate max dust fraction based on particle data

    with warnings.catch_warnings():
        warnings.simplefilter('ignore') # Ignore warning if file is empty
        gas_file = np.atleast_2d(np.loadtxt(txtFilePath + 'snap' + snapNum + '_ID' + haloID + '_gas.txt')) # Calculate dust surface density from the 
        # original gas particle data, to avoid issues with negative dust masses due to TODDLERS dust subtraction
    
    if np.shape(gas_file) == (1, 10): # Only one gas particle

        maxDustFraction = 10**(-4.5)
    
    elif np.size(gas_file, axis = 1) > 0: # 2 or more gas particles

        dust_r = np.sqrt(gas_file[:, 0]**2 + gas_file[:, 1]**2 + gas_file[:, 2]**2) * 1e-3 # In kpc
        dust_m = np.sum(gas_file[:, 7:], axis = 1) # In Msun


        dustMasses_sorted = dust_m[np.argsort(dust_r)]

        idx_halfmass = np.min(np.argwhere((np.cumsum(dustMasses_sorted) / np.sum(dustMasses_sorted)) >= 0.5))

        dustHalfMassRadius = np.sort(dust_r)[idx_halfmass]

        dustHalfMass = (np.sum(dust_m) / 2.)

        SigmaDust = dustHalfMass / (np.pi * dustHalfMassRadius**2) # In solar masses / kpc^2

        maxDustFraction = np.clip(10**(-0.5 - np.log10(SigmaDust)), a_min = 10**(-6.5), a_max = 10**(-4.5))

        subprocess.run(['perl', '-pi', '-e', 's/maxLevel=\"0/maxLevel=\"' + str(binTreeMaxLevel) + '/g', skifilename_halo])

        subprocess.run(['perl', '-pi', '-e', 's#dust.txt#' + SKIRTinputFiles + '_dust.txt#g', skifilename_halo])

        subprocess.run(['perl', '-pi', '-e', 's/minX=\"-0/minX=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
        subprocess.run(['perl', '-pi', '-e', 's/maxX=\"0/maxX=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
        subprocess.run(['perl', '-pi', '-e', 's/minY=\"-0/minY=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
        subprocess.run(['perl', '-pi', '-e', 's/maxY=\"0/maxY=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
        subprocess.run(['perl', '-pi', '-e', 's/minZ=\"-0/minZ=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
        subprocess.run(['perl', '-pi', '-e', 's/maxZ=\"0/maxZ=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])

        subprocess.run(['perl', '-pi', '-e', 's/maxDustFraction=\"0/maxDustFraction=\"' + str(maxDustFraction) + '/g', skifilename_halo])

    else:
    
        # Change the skirt simulation to noMedium

        subprocess.run(['perl', '-pi', '-e', 's/simulationMode=\"DustEmission/simulationMode=\"NoMedium/g', skifilename_halo])

        with open(skifilename_halo, 'r+') as fp:
            # read and store all lines into list
            lines = fp.readlines()
            # move file pointer to the beginning of a file
            fp.seek(0)
            # truncate the file
            fp.truncate()

            # start writing lines
            # iterate line and line number
            for number, line in enumerate(lines):
                if number <= mediumXMLlineDict[skifileversion][0] or number >= mediumXMLlineDict[skifileversion][1]:
                    fp.write(line)


    subprocess.run(['perl', '-pi', '-e', 's/numPackets=\"0/numPackets=\"' + str(Npp) + '/g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's#old_stars#' + SKIRTinputFiles + '_old_stars#g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's#starforming_gas#' + SKIRTinputFiles + '_starforming_gas#g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/Period0/Period' + str(int(old_stars_tmin.to('Myr').value)) + '/g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's/radius=\"1 Rstar/radius=\"' + str(Rstar) + ' kpc' +  '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/radius=\"3 Rstar/radius=\"' + str(3. * Rstar) + ' kpc' +  '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/radius=\"5 Rstar/radius=\"' + str(5. * Rstar) + ' kpc' +  '/g', skifilename_halo])

    return None

editSki(snapNum, haloID, Rstar)

print('Elapsed time to edit ski file and calculate dust surface density:', datetime.now() - startTime)