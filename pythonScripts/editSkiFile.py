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

def delete_medium_system(ski_file):
    """
    Delete the entire block from <mediumSystem ...> to </mediumSystem> in a .ski file when 'NoMedium' is used
    """

    # Read all lines from the file
    with open(ski_file, 'r') as f:
        lines = f.readlines()

    new_lines = []
    skip = False

    for line in lines:

        # Detect the start tag
        if '<mediumSystem' in line:
            skip = True
            continue

        # Detect the end tag
        if '</mediumSystem>' in line:
            skip = False
            continue

        # Keep lines that are not within the mediumSystem block
        if not skip:
            new_lines.append(line)


    # Overwrite the original file with the new content
    with open(ski_file, 'w') as f:
        f.writelines(new_lines)

startTime = datetime.now()

# Global settings

old_stars_tmin = unyt.unyt_quantity(10., 'Myr') # Minimum age in Myr for an evolved star particle. Also determines the TODDLERS averaging timescale
# Don't change this unless you know what you're doing :)

Npp = int(10**4.5) # Number of photon packets
binTreeMaxLevel = 36 # Max refinement level of the spatial grid

snapNum = sys.argv[1]
haloID = sys.argv[2]
Rstar = float(sys.argv[3])

txtFilePath = sys.argv[4]
SKIRTinputFilePath = sys.argv[5]

f = open(txtFilePath + 'snap' + snapNum + '_' + 'ID' + haloID + '_stars.txt', 'r')
header = f.readline() # Read first header line
redshift = float(header.split(' ')[-1])
f.close()

scaleFactor = 1. / (1. + redshift) # Scale factor for the snapshot
SKIRTboxsize = unyt.unyt_quantity(min(100., 100. * 1.8 / 0.7 * scaleFactor), 'kpc') # Scale SKIRT box size akin to COLIBRE gravitational softening length

skifileversion = '5.0'


# Define filepaths from parameter file
dir_path = os.path.dirname(os.path.realpath(__file__))
with open(f'{dir_path}/../SKIRT_parameters.yml','r') as stream:
    params = yaml.safe_load(stream)

# Edit ski file

def editSki(snapNum, haloID, Rstar):

    SKIRTinputFiles = SKIRTinputFilePath + 'snap' + snapNum + '_ID' + haloID

    skifilename = params['SkirtFilepaths']['skiFilepath'].format(skifileversion=skifileversion)

    skifilename_halo = 'snap' + snapNum + '_ID' + haloID + '.ski'


    subprocess.run(['cp', skifilename, skifilename_halo]) # copy the skirt file for each galaxy

    # Calculate max dust fraction based on particle data

    subprocess.run(['perl', '-pi', '-e', 's/maxLevel=\"0/maxLevel=\"' + str(binTreeMaxLevel) + '/g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's#dust.txt#' + SKIRTinputFiles + '_dust.txt#g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's/minX=\"-0/minX=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/maxX=\"0/maxX=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/minY=\"-0/minY=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/maxY=\"0/maxY=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/minZ=\"-0/minZ=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/maxZ=\"0/maxZ=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])


    with warnings.catch_warnings():
        warnings.simplefilter('ignore') # Ignore warning if file is empty
        gas_file = np.atleast_2d(np.loadtxt(txtFilePath + 'snap' + snapNum + '_ID' + haloID + '_gas.txt')) # Calculate dust surface density from the 
        # original gas particle data, to avoid issues with negative dust masses due to TODDLERS dust subtraction

    if np.size(gas_file, axis = 1) > 0: # if there are gas particles

        dust_r = np.sqrt(gas_file[:, 0]**2 + gas_file[:, 1]**2 + gas_file[:, 2]**2) * 1e-3 # In kpc

        dust_m = np.sum(gas_file[:, 10:], axis = 1) # In Msun

        if (dust_m>0).sum() >= 2: # if there are more than 2 gas particles containing dust

            dustMasses_sorted = dust_m[np.argsort(dust_r)]

            idx_halfmass = np.min(np.argwhere((np.cumsum(dustMasses_sorted) / np.sum(dustMasses_sorted)) >= 0.5))

            dustHalfMassRadius = np.sort(dust_r)[idx_halfmass]

            dustHalfMass = (np.sum(dust_m) / 2.)

            SigmaDust = dustHalfMass / (np.pi * dustHalfMassRadius**2) # In solar masses / kpc^2

            maxDustFraction = np.clip(10**(-0.5 - np.log10(SigmaDust)), a_min = 10**(-6.5), a_max = 10**(-4.5))

            subprocess.run(['perl', '-pi', '-e', 's/maxDustFraction=\"0/maxDustFraction=\"' + str(maxDustFraction) + '/g', skifilename_halo])

        

        elif (dust_m>0).sum()==1: # if there are only one gas particle containing dust

            maxDustFraction = 10**(-4.5)

            subprocess.run(['perl', '-pi', '-e', 's/maxDustFraction=\"0/maxDustFraction=\"' + str(maxDustFraction) + '/g', skifilename_halo])



        else:# if there is gas particles but dust mass is all 0

            subprocess.run(['perl', '-pi', '-e', 's/simulationMode="[^"]*"/simulationMode="NoMedium"/g', skifilename_halo])

            delete_medium_system(skifilename_halo)
            
    else:

        # Change the skirt simulation to noMedium
        
        subprocess.run(['perl', '-pi', '-e', 's/simulationMode="[^"]*"/simulationMode="NoMedium"/g', skifilename_halo])

        delete_medium_system(skifilename_halo)




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