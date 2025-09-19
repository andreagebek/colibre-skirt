"""
Edit a ski file
"""

import numpy as np
import subprocess
import sys
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

Npp = int(10**7.5) # Number of photon packets
binTreeMaxLevel = 36 # Max refinement level of the spatial grid

snapNum = sys.argv[1]
haloID = sys.argv[2]
Rstar = unyt.unyt_quantity(float(sys.argv[3]), 'kpc')
Mdust = unyt.unyt_quantity(float(sys.argv[4]), 'Msun')

txtFilePath = sys.argv[5]
SKIRTinputFilePath = sys.argv[6]

Resolution = sys.argv[7]

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

    subprocess.run(['perl', '-pi', '-e', 's/maxLevel=\"0/maxLevel=\"' + str(binTreeMaxLevel) + '/g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's#dust.txt#' + SKIRTinputFiles + '_dust.txt#g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's/minX=\"-0/minX=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/maxX=\"0/maxX=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/minY=\"-0/minY=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/maxY=\"0/maxY=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/minZ=\"-0/minZ=\"' + str(-SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/maxZ=\"0/maxZ=\"' + str(SKIRTboxsize.to('pc').value / 2.) + '/g', skifilename_halo])

    if Mdust == 0: # No dust within 50-kpc exclusive sphere aperture -> run SKIRT without medium

        subprocess.run(['perl', '-pi', '-e', 's/simulationMode="[^"]*"/simulationMode="NoMedium"/g', skifilename_halo])

        delete_medium_system(skifilename_halo)

    else:

        SigmaDust = Mdust / (np.pi * Rstar**2) # Dust surface density

        if Resolution == '5':

            maxDustFraction = np.clip(10**(-0.5 - np.log10(SigmaDust)), a_min = 10**(-6.5), a_max = 10**(-4.5))
            # Slightly higher resolved spatial grid at m5 resolution

        else: # Same for m6 and m7 for now

            maxDustFraction = np.clip(10**(-1.5 - 0.75 * np.log10(SigmaDust)), a_min = 1e-6, a_max = 10**(-4.5))

        subprocess.run(['perl', '-pi', '-e', 's/maxDustFraction=\"0/maxDustFraction=\"' + str(maxDustFraction) + '/g', skifilename_halo])


    subprocess.run(['perl', '-pi', '-e', 's/numPackets=\"0/numPackets=\"' + str(Npp) + '/g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's#old_stars#' + SKIRTinputFiles + '_old_stars#g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's#starforming_gas#' + SKIRTinputFiles + '_starforming_gas#g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/Period0/Period' + str(int(old_stars_tmin.to('Myr').value)) + '/g', skifilename_halo])

    subprocess.run(['perl', '-pi', '-e', 's/radius=\"1 Rstar/radius=\"' + str(Rstar.to('kpc').value) + ' kpc' +  '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/radius=\"3 Rstar/radius=\"' + str(3. * Rstar.to('kpc').value) + ' kpc' +  '/g', skifilename_halo])
    subprocess.run(['perl', '-pi', '-e', 's/radius=\"5 Rstar/radius=\"' + str(5. * Rstar.to('kpc').value) + ' kpc' +  '/g', skifilename_halo])

    return None

editSki(snapNum, haloID, Rstar)

print('Elapsed time to edit ski file:', datetime.now() - startTime)