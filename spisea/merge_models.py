### Code used to incorporate new brown dwarf evolution/atmosphere models into pre-existing code.
### Messy, ask Caitlin Begbie if there are questions

import math
import logging
from numpy import genfromtxt
import numpy as np
import os
import glob
import pandas as pd
import pdb
import warnings
from astropy.table import Table, vstack, Column
from scipy import interpolate
import pylab as py
from spisea.utils import objects
from scipy.interpolate import RegularGridInterpolator
from spisea import exceptions
import re
import matplotlib
import matplotlib.pyplot as plt
from spisea import atmospheres

logger = logging.getLogger('evolution')

# Fetch root directory of evolution models.
try:
    models_dir = os.environ['SPISEA_MODELS']
    models_dir += '/evolution/'
except KeyError:
    warnings.warn("SPISEA_MODELS is undefined; functionality "
                  "will be SEVERELY crippled.")
    models_dir = ''

# Code to deconstruct current Phillips2020 evolution files and reconstruct iso files based on age

# Define input and output directories
input_dir = '/System/Volumes/Data/mnt/g/lu/models/evolution/Phillips2020/z00_mass'
output_dir = '/System/Volumes/Data/mnt/g/lu/models/evolution/Phillips2020/z00_age'
os.makedirs(output_dir, exist_ok=True)

# Combine data from all files
combined_data = []

for filename in os.listdir(input_dir):
    if filename.endswith('.txt'):
        file_path = os.path.join(input_dir, filename)
        table = Table.read(file_path, format='ascii')
        combined_data.append(table)

# Concatenate all data into a single table
all_data = combined_data[0]

# Add the remaining tables to the main table
for table in combined_data[1:]:
    all_data = vstack([all_data, table])  # Use vstack to stack tables

# Group data by age
ages = np.unique(all_data['Age'])

# Iterate over each unique age and filter data
for age in ages:
    age_group = all_data[all_data['Age'] == age]
    output_path = os.path.join(output_dir, f"ATMO_{age}.txt")

    # Save each age group as a .txt file
    age_group.write(output_path, format='ascii', overwrite=True)

print(f"New files created in: {output_dir}")

# Code to eradicate duplicate lines seen in above code
"""new_input_dir = output_dir

for filename2 in os.listdir(new_input_dir):
    if filename2.endswith('.txt'):
        file_path2 = os.path.join(new_input_dir, filename2)

    # Remove duplicates before loading as a table
    with open(file_path2, 'r') as f:
        lines = f.readlines()

    header = lines[0]
    data_lines = [line for line in lines[1:] if not line.startswith(header.split()[0])]

    with open(file_path2, 'w') as f:
            f.write(header)
            f.writelines(data_lines)

    # Now process with Astropy Table
    table2 = Table.read(file_path2, format='ascii')
    unique_table2 = unique(table2)
    unique_table2.write(file_path2, format='ascii', overwrite=True)

print('files overwritten successfully!')


# Reformatting files to match Parsec
def reformat():
    for file in os.listdir(output_dir):
        if file.endswith('.txt'):
            new_fp = os.path.join(output_dir, file)

        # Add metallicity column of all 0.0
        r_table = Table.read(new_fp, format='ascii')
        r_table.add_column(0.0, name='Z', index=0)

        # Duplicate mass to include current mass column
        r_table.add_column(r_table['Mass'], name='Mass_current')

        # Update age column to make it log scale
        r_table['Age'] = np.log10(r_table['Age'] * 1e9)  #originally in Gyr

        # Update Teff column to make it log scale
        r_table['Teff'] = np.log10(r_table['Teff'])

        # Reorder columns to match Parsec
        new_order = ['Z', 'Age', 'Mass', 'Mass_current', 'Luminosity', 'Teff', 'Gravity', 'Gaia_Gbp', 'Gaia_G', 'Gaia_Grp']
        t_new = r_table[new_order]
        t_new.write(new_fp, format='ascii', overwrite=True)

        print(f"{file} reformatted successfully!")

    return

# Transform reformatted .txt files to iso fits files
def ATMO_to_iso():
    i_dir = output_dir
    o_dir = '/System/Volumes/Data/mnt/g/lu/models/evolution/Phillips2020/iso' #output directory SPISEA will pull from

    for file in os.listdir(i_dir):
        if file.endswith('.txt'):
            fp = os.path.join(i_dir, file)

            try:
                age_str = file.split('_')[1].split('.txt')[0]
                age = float(age_str)
                log_age = np.log10(age * 1e9)
            except (IndexError, ValueError):
                print(f"{file} cannot be processed.")
                continue

            try:
                # load in .txt file as ascii table
                table = Table.read(fp, format='ascii')
            except Exception as e:
                print(f"{fp} could not be read: {e}")
                continue

            #create output file name
            o_file = os.path.join(o_dir, f'iso_{log_age}.fits')

            try:
                # write to output directory as fits file
                table.write(o_file, format='fits', overwrite=True)
                print(f"{o_file} created successfully!")
            except Exception as e:
                print(f"{o_file} could not be generated: {e}")
                continue

    return


"""
### CODE TO UNPACK MARLEY FILES ###
def Marley_deconstruct():
    """
    Code to deconstruct big Marley files into individual age files
    """
    in_dir = '/System/Volumes/Data/mnt/g/lu/models/evolution/Marley2021/Marley_age'
    out_dir = '/System/Volumes/Data/mnt/g/lu/models/evolution/Marley2021/iso'

    # Set path to each file
    for filename in os.listdir(in_dir):
        file_path = os.path.join(in_dir, filename)

        # Extract metallicity from filename
        match = re.search(r'nc([+-]\d+\.\d+)_co', filename)
        metallicity = match.group(1) if match else "unknown"

        # Deconstruct initial format of single column table
        with open(file_path, 'r') as f:
            lines = f.readlines()
            header = ['Age', 'Mass', 'log_L', 'Teff', 'logg', 'Radius']
            data = [line.strip().split() for line in lines[1:]]
            data = [row for row in data if len(row) > 1]
            table = Table(rows=data, names=header)

        # Find unique ages to create iso age_based files
        ages = np.unique(table['Age'])
    
        for age in ages:
            age_group = table[table['Age'] == age]
            out_path = os.path.join(out_dir, f"Marley_{metallicity}_{age}.txt")

            # Save each age group as a .txt file
            age_group.write(out_path, format='ascii', overwrite=True)

        print(f"New files created in: {out_dir}")
        
    return

# Reformatting files to match Parsec
def reformat_Marley():
    age_dir = '/System/Volumes/Data/mnt/g/lu/models/evolution/Marley2021/age_txt'
    for file in os.listdir(age_dir):
        if file.endswith('.txt'):
            new_fp = os.path.join(age_dir, file)

        # Add metallicity column
        table = Table.read(new_fp, format='ascii')
        
        # Extract metallicity from filename
        match = re.search(r'Marley_([+-]?\d+\.\d+)_\d+\.\d+.txt', file)
        metallicity = float(match.group(1)) if match else np.nan

        # Replace or add metallicity column
        if 'Z' in table.colnames:
            table.replace_column('Z', [metallicity] * len(table))
        else:
            table.add_column([metallicity] * len(table), name='Z', index=0)

        # Duplicate mass to include current mass column
        #table.add_column(table['Mass'], name='Mass_current')

        # Update age column to make it log scale
        #table['Age'] = np.log10(table['Age'] * 1e9)  #originally in Gyr

        # Update Teff column to make it log scale
        #table['Teff'] = np.log10(table['Teff'])

        # Reorder columns to match Parsec
        new_order = ['Z', 'Age', 'Mass', 'Mass_current', 'log_L', 'Teff', 'logg', 'Radius']
        t_new = table[new_order]
        t_new.write(new_fp, format='ascii', overwrite=True)

        print(f"{file} reformatted successfully!")

    return

def Marley_to_iso():
    in_dir = '/System/Volumes/Data/mnt/g/lu/models/evolution/Marley2021/age_txt'
    out_dir_1 = '/System/Volumes/Data/mnt/g/lu/models/evolution/Marley2021/iso/zm05'
    out_dir_2 = '/System/Volumes/Data/mnt/g/lu/models/evolution/Marley2021/iso/zp00'
    out_dir_3 = '/System/Volumes/Data/mnt/g/lu/models/evolution/Marley2021/iso/zp05'

    for file in os.listdir(in_dir):
        if file.endswith('.txt'):
            fp = os.path.join(in_dir, file)

            try:
                parts = file.split('_')
                metal_str = parts[1]  # Second part is metallicity
                age_str = parts[2].split('.txt')[0]
                
                age = float(age_str)
                metallicity = float(metal_str)
                log_age = np.log10(age * 1e9)
            except (IndexError, ValueError):
                print(f"{file} cannot be processed.")
                continue

            try:
                # load in .txt file as ascii table
                table = Table.read(fp, format='ascii')
            except Exception as e:
                print(f"{fp} could not be read: {e}")
                continue

            # Determine output directory based on metallicity
            if metallicity == -0.5:
                out_dir = out_dir_1
            elif metallicity == 0.0:
                out_dir = out_dir_2
            elif metallicity == 0.5:
                out_dir = out_dir_3
            else:
                print(f"Skipping {file}: Unexpected metallicity value ({metallicity}).")
                continue

            # Create output filename
            out_file = os.path.join(out_dir, f'iso_{log_age}.fits')

            try:
                # Write to output directory as FITS file
                table.write(out_file, format='fits', overwrite=True)
                print(f"{out_file} created successfully!")
            except Exception as e:
                print(f"Error writing {out_file}: {e}")
                continue
                
    return

### CREATING MERGED EVO MODEL WITH PHILLIPS ###
def get_phillips_isochrone(logAge, metallicity='solar'):
    """
    Load mass, effective temperature, log gravity, and log luminosity
    for the Phillips isochrones at given logAge. Code will quit if that
    logAge value doesn't exist (can make some sort of interpolation thing
    later).

    Note: mass is currently initial mass, not instantaneous mass
    
    Inputs:
    logAge - Logarithmic Age
    metallicity - in Z (def = solar of 0.014)
    """
    rootDir = models_dir + 'Phillips2020/iso/'
    metSuffix = 'z00/'
    if metallicity != 'solar':
        print( 'Non-solar Phillips 2020 metallicities not supported yet')
        return
    rootDir += metSuffix

    # List available isochrone files
    available_ages = []
    for filename in os.listdir(rootDir):
        if filename.startswith('iso_') and filename.endswith('.fits'):
            age_str = filename.split('_')[1].split('.fits')[0]  # Extract the age part from filename
            available_ages.append(float(age_str))
    
    # Find the closest available age from Phillips
    closest_age = min(available_ages, key=lambda x: abs(x - logAge))
    print(closest_age)
    
    # Load the corresponding isochrone
    isoFile = rootDir + f'iso_{closest_age}.fits'
    print(f"Loading isochrone for logAge = {closest_age}")

    # Check to see if isochrone exists
    if not os.path.exists(isoFile):
        print( f'Phillips isochrone for logAge = {closest_age} does\'t exist')
        print( 'Quitting')
        return
        
    data = Table.read(isoFile, format='fits')
    print("Available Columns:", data.keys())
    cols = data.keys()
    mass = data[cols[2]] #Note: this is initial mass, in M_sun
    logT = data[cols[5]] # K
    logL = data[cols[4]] # L_sun
    logg = data[cols[6]]
    mass_current = data[cols[3]] # Matches initial mass -- BD masses assumed to not change w current models
    phase = np.ones(len(mass), dtype=int)
    
    obj = objects.DataHolder()
    obj.mass = mass
    obj.logT = logT
    obj.logg = logg
    obj.logL = logL
    obj.mass_current = mass_current
    obj.phase = phase

    return obj

def get_Baraffe15_isochrone(logAge, metallicity='solar'):
    """
    Load mass, effective temperature, log gravity, and log luminosity
    for the Baraffe+15 isochrones at given logAge. Code will quit if that
    logAge value doesn't exist (can make some sort of interpolation thing
    later).

    ALSO interpolates isochrone to a finer mass grid.

    Inputs:
    logAge - Logarithmic Age
    metallicity - in Z (def = solar of 0.014)
    """
    rootDir = models_dir + 'Baraffe15/iso/'
    if metallicity != 'solar':
        print( 'Non-solar Baraffe+15 metallicities not supported yet')
        return

    # Check to see if isochrone exists
    isoFile = rootDir + 'iso_%.2f.fits' % logAge
    if not os.path.exists(isoFile):
        print( 'Baraffe+15 isochrone for logAge = {0:3.2f} does\'t exist'.format(logAge))
        print( 'Quitting')
        return
        
    data = Table.read(isoFile, format='fits')
    mass = data['Mass'] #Note: this is initial mass, in M_sun
    logT = np.log10(data['Teff']) # K
    logL = data['logL'] # L_sun
    logg = data['logG']

    # Interpolate isochrone to finer mass grid. Spacing
    # is one model every 0.02 M_sun down to 0.2 M_sun, then
    # one model every 0.005 M_sun down to 0.07 M_sun
    new_masses0 = np.arange(min(mass), 0.1, 0.005)
    new_masses2 = np.arange(0.1, max(mass), 0.02)
    
    #new_masses = np.concatenate((new_masses0, new_masses1, new_masses2))
    new_masses = np.concatenate((new_masses0, new_masses2))
    
    # Build interpolators in linear space
    f_logT = interpolate.interp1d(mass, 10**logT, kind='linear')
    f_logL = interpolate.interp1d(mass, 10**logL, kind='linear')
    f_logg = interpolate.interp1d(mass, 10**logg, kind='linear')

    # Do interpolation, convert back to logspace
    logT_interp = np.log10(f_logT(new_masses))
    logL_interp = np.log10(f_logL(new_masses))
    logg_interp = np.log10(f_logg(new_masses))

    # Hack, add new mass_current and phase columns.
    mass_current = np.array(new_masses)
    phase = np.ones(len(new_masses), dtype=int)

    # Test the interpolation, if desired
    test = False
    if test:
        py.figure(1, figsize=(10,10))
        py.clf()
        py.plot(mass, logT, 'k.', ms = 10, label='Orig')
        py.plot(new_masses, logT_interp, 'r.', ms=7, label='Interp')
        py.xlabel('Mass')
        py.ylabel('logT')
        py.legend()

        py.figure(2, figsize=(10,10))
        py.clf()
        py.plot(mass, logL, 'k.', ms = 10, label='Orig')
        py.plot(new_masses, logL_interp, 'r.', ms=7, label='Interp')
        py.xlabel('Mass')
        py.ylabel('logL')
        py.legend()

        py.figure(3, figsize=(10,10))
        py.clf()
        py.plot(mass, logg, 'k.', ms = 10, label='Orig')
        py.plot(new_masses, logg_interp, 'r.', ms=7, label='Interp')
        py.xlabel('Mass')
        py.ylabel('logg')
        py.legend()

        pdb.set_trace()
    
    # Make isochrone
    obj = objects.DataHolder()
    obj.mass = new_masses
    obj.logT = logT_interp
    obj.logg = logg_interp
    obj.logL = logL_interp
    obj.mass_current = mass_current
    obj.phase = phase

    return obj

def merge_isochrone_baraffe_phillips(logAge, metallicity='solar'):
    """
    Function to merge Baraffe+15 and Phillips 2020 models. Will take
    100% Phillips2020 between 0.01 - 0.07 M_sun, transition between
    0.07 - 0.075 M_sun, and take 100% Baraffe from 0.075 M_sun and up.

    Can only handle ages at which models already exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    """
    if metallicity != 'solar':
        print( 'Non-solar metallicity not supported yet')
        return

    # Get individual Baraffe and Phillips isochrones at desired age. Note
    # that this will also give the Baraffe models a finer mass sampling
    isoBaraffe = get_Baraffe15_isochrone(logAge, metallicity=metallicity)
    isoPhillips = get_phillips_isochrone(logAge, metallicity=metallicity)

    # Identify M >= 0.075 M_sun in Baraffe and M <= 0.07 M_sun in Phillips
    good_b = np.where(isoBaraffe.mass >= 0.075)
    good_p = np.where(isoPhillips.mass <= 0.07)

    # Sample between 0.4 M_sun and 0.5 M_sun in steps of 0.02 M_sun.
    # Will do linear combo of Baraffe and Phillips over this range
    mid_mass = np.arange(0.07, 0.075, 0.002)
    mid_logT = []
    mid_logL = []
    mid_logG = []
    mid_Mcurr = []
    mid_phase = []
    for mass in mid_mass:
        # Find the appropriate masses in Baraffe + Phillips to build from.
        # Baraffe has identical sampling over this range, and Phillips sampling
        # is very close. As a result, we will just take the closest mass model
        # to each mid_mass
        idx_b = np.where( abs(isoBaraffe.mass - mass)  == min(abs(isoBaraffe.mass - mass)) )
        idx_p = np.where( abs(isoPhillips.mass - mass)  == min(abs(isoPhillips.mass - mass)) )

        # Quality control check: we won't let the difference between model mass and
        # chosen mass to be >= 0.01 M_sun
        if ((isoPhillips.mass[idx_p] - mass) >= 0.002) | ((isoBaraffe.mass[idx_b] - mass) >= 0.002):
            print( 'WARNING: Baraffe or Phillips model interpolation between 0.01 - 0.08 M_sun may \
            be inaccurate. Check this!')
            pdb.set_trace()
    
        # Now, do the linear combo of models at this mass, weighted by distance from
        # 0.07 or 0.075 (whichever is appropriate)
        diff = 0.075 - 0.07
        weight_p = (0.075 - mass) / diff
        weight_b = 1.0 - weight_p
        print( 'Baraffe {0} and Phillips {1} at mass {2}'.format(weight_p, weight_b, mass))

        # Now, do the merge IN LINEAR SPACE!
        Teff = (10**isoPhillips.logT[idx_p] * weight_p) + \
            (10**isoBaraffe.logT[idx_b] * weight_b)
        L = (10**isoPhillips.logL[idx_p] * weight_p) + \
            (10**isoBaraffe.logL[idx_b] * weight_b)
        g = (10**isoPhillips.logg[idx_p] * weight_p) + \
            (10**isoBaraffe.logg[idx_b] * weight_b)
        mcurr = (isoPhillips.mass_current[idx_p] * weight_p) + \
                (isoBaraffe.mass_current[idx_b] * weight_b)
        phase = np.round((isoPhillips.mass_current[idx_p] * weight_p) + \
                         (isoBaraffe.mass_current[idx_b] * weight_b))
        
        mid_logT = np.concatenate((mid_logT, np.log10(Teff)))
        mid_logL = np.concatenate((mid_logL, np.log10(L)))
        mid_logG = np.concatenate((mid_logG, np.log10(g)))
        mid_Mcurr = np.concatenate((mid_Mcurr, mcurr))
        mid_phase = np.concatenate((mid_phase, phase))

    # Now, final isochrone will be combination of Baraffe at M>=0.075,
    # Phillips at M<=0.075, and the combination inbetween
    mass = np.concatenate((isoPhillips.mass[good_p], mid_mass, isoBaraffe.mass[good_b]))
    logT = np.concatenate((isoPhillips.logT[good_p], mid_logT, isoBaraffe.logT[good_b]))
    logL = np.concatenate((isoPhillips.logL[good_p], mid_logL, isoBaraffe.logL[good_b]))
    logG = np.concatenate((isoPhillips.logg[good_p], mid_logG, isoBaraffe.logg[good_b]))
    mcurr = np.concatenate((isoPhillips.mass_current[good_p], mid_Mcurr, isoBaraffe.mass_current[good_b]))
    phase = np.concatenate((isoPhillips.phase[good_p], mid_phase, isoBaraffe.phase[good_b]))

    # Also add a source flag
    source = np.concatenate( (['Phillips']*len(good_p[0]), 
                              ['Baraffe+Phillips']*len(mid_mass),
                              ['Baraffe']*len(good_b[0])) )

    iso = objects.DataHolder()
    iso.mass = mass
    iso.logL = logL
    iso.logg = logG
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.source = source

    return iso

def merge_isochrone_baraffe_pisa(logAge, metallicity='solar'):
    """
    Function to merge Baraffe+15 and Pisa 2011 models. Will take
    100% Baraffe+15 between 0.07 - 0.4 M_sun, transition between
    0.4 - 0.5 M_sun, and take 100% Pisa from 0.5 M_sun and up.

    Can only handle ages at which models already exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    """
    if metallicity != 'solar':
        print( 'Non-solar metallicity not supported yet')
        return

    # Get individual Baraffe and Pisa isochrones at desired age. Note
    # that this will also give the Baraffe models a finer mass sampling
    isoBaraffe = get_Baraffe15_isochrone(logAge, metallicity=metallicity)
    isoPisa = get_pisa_isochrone(logAge, metallicity=metallicity)

    # Identify M <= 0.4 M_sun in Baraffe and M >= 0.5 M_sun in Pisa
    good_b = np.where(isoBaraffe.mass <= 0.4)
    good_p = np.where(isoPisa.mass >= 0.5)

    # Sample between 0.4 M_sun and 0.5 M_sun in steps of 0.02 M_sun.
    # Will do linear combo of Baraffe and Pisa over this range
    mid_mass = np.arange(0.4, 0.5+0.01, 0.02)
    mid_logT = []
    mid_logL = []
    mid_logG = []
    mid_Mcurr = []
    mid_phase = []
    for mass in mid_mass:
        # Find the appropriate masses in Baraffe + Pisa to build from.
        # Baraffe has identical sampling over this range, and Pisa sampling
        # is very close. As a result, we will just take the closest mass model
        # to each mid_mass
        idx_b = np.where( abs(isoBaraffe.mass - mass)  == min(abs(isoBaraffe.mass - mass)) )
        idx_p = np.where( abs(isoPisa.mass - mass)  == min(abs(isoPisa.mass - mass)) )

        # Quality control check: we won't let the difference between model mass and
        # chosen mass to be >= 0.01 M_sun
        if ((isoPisa.mass[idx_p] - mass) >= 0.02) | ((isoBaraffe.mass[idx_b] - mass) >= 0.02):
            print( 'WARNING: Baraffe or Pisa model interpolation between 0.4 - 0.5 M_sun may \
            be inaccurate. Check this!')
            pdb.set_trace()
    
        # Now, do the linear combo of models at this mass, weighted by distance from
        # 0.4 or 0.5 (whichever is appropriate)
        diff = 0.5 - 0.4
        weight_b = (0.5 - mass) / diff
        weight_p = 1.0 - weight_b
        print( 'Baraffe {0} and Pisa {1} at mass {2}'.format(weight_b, weight_p, mass))

        # Now, do the merge IN LINEAR SPACE!
        Teff = (10**isoBaraffe.logT[idx_b] * weight_b) + \
            (10**isoPisa.logT[idx_p] * weight_p)
        L = (10**isoBaraffe.logL[idx_b] * weight_b) + \
            (10**isoPisa.logL[idx_p] * weight_p)
        g = (10**isoBaraffe.logg[idx_b] * weight_b) + \
            (10**isoPisa.logg[idx_p] * weight_p)
        mcurr = (isoBaraffe.mass_current[idx_b] * weight_b) + \
                (isoPisa.mass_current[idx_p] * weight_p)
        phase = np.round((isoBaraffe.mass_current[idx_b] * weight_b) + \
                         (isoPisa.mass_current[idx_p] * weight_p))
        
        mid_logT = np.concatenate((mid_logT, np.log10(Teff)))
        mid_logL = np.concatenate((mid_logL, np.log10(L)))
        mid_logG = np.concatenate((mid_logG, np.log10(g)))
        mid_Mcurr = np.concatenate((mid_Mcurr, mcurr))
        mid_phase = np.concatenate((mid_phase, phase))

    # Now, final isochrone will be combination of Baraffe at M<=0.4,
    # Pisa at M>=0.5, and the combination inbetween
    mass = np.concatenate((isoBaraffe.mass[good_b], mid_mass, isoPisa.mass[good_p]))
    logT = np.concatenate((isoBaraffe.logT[good_b], mid_logT, isoPisa.logT[good_p]))
    logL = np.concatenate((isoBaraffe.logL[good_b], mid_logL, isoPisa.logL[good_p]))
    logG = np.concatenate((isoBaraffe.logg[good_b], mid_logG, isoPisa.logg[good_p]))
    mcurr = np.concatenate((isoBaraffe.mass_current[good_b], mid_Mcurr, isoPisa.mass_current[good_p]))
    phase = np.concatenate((isoBaraffe.phase[good_b], mid_phase, isoPisa.phase[good_p]))

    # Also add a source flag
    source = np.concatenate( (['Baraffe']*len(good_b[0]), ['Baraffe+Pisa']*len(mid_mass),
                              ['Pisa']*len(good_p[0])) )

    iso = objects.DataHolder()
    iso.mass = mass
    iso.logL = logL
    iso.logg = logG
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.source = source

    return iso

"""def plot_merged_isochrone(iso):
    ""
    Function to plot the merged isochrone model to visualize the transition
    between Phillips 2020 and Baraffe+15 models.

    Parameters:
    -----------
    iso : DataHolder object
        Merged isochrone returned by merge_isochrone_baraffe_phillips().
    ""

    # Define colors for different sources
    color_map = {'Phillips': 'blue', 'Baraffe': 'red', 'Baraffe+Phillips': 'pink'}

    # Create figure and subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Mass vs. logT (Effective Temperature)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[0].scatter(iso.mass[mask], iso.logT[mask], label=source, color=color_map[source], alpha=0.7)
        axes[0].plot(iso.mass[mask], iso.logT[mask], label='linear representation')
    axes[0].set_xlabel("Mass ($M_{\odot}$)")
    axes[0].set_ylabel("$\log T_{\mathrm{eff}}$ (K)")
    axes[0].set_title("Mass vs. Effective Temperature")
    axes[0].axvline(0.07, linestyle="--", color="gray", alpha=0.5)  # Transition marker
    axes[0].axvline(0.08, linestyle="--", color="gray", alpha=0.5)  # Transition marker
    axes[0].set_xlim(0.0, 0.1)
    axes[0].legend()

    # Mass vs. logL (Luminosity)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[1].scatter(iso.mass[mask], iso.logL[mask], label=source, color=color_map[source], alpha=0.7)
        axes[1].plot(iso.mass[mask], iso.logL[mask], label='linear representation')
    axes[1].set_xlabel("Mass ($M_{\odot}$)")
    axes[1].set_ylabel("$\log L$ ($L_{\odot}$)")
    axes[1].set_title("Mass vs. Luminosity")
    axes[1].axvline(0.07, linestyle="--", color="gray", alpha=0.5)
    axes[1].axvline(0.08, linestyle="--", color="gray", alpha=0.5)
    axes[1].set_xlim(0.0, 0.1)
    axes[1].legend()

    # Mass vs. logg (Surface Gravity)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[2].scatter(iso.mass[mask], iso.logg[mask], label=source, color=color_map[source], alpha=0.7)
        axes[2].plot(iso.mass[mask], iso.logg[mask], label='linear representation')
    axes[2].set_xlabel("Mass ($M_{\odot}$)")
    axes[2].set_ylabel("$\log g$ (cm/s²)")
    axes[2].set_title("Mass vs. Surface Gravity")
    axes[2].axvline(0.07, linestyle="--", color="gray", alpha=0.5)
    axes[2].axvline(0.08, linestyle="--", color="gray", alpha=0.5)
    axes[2].set_xlim(0.0, 0.1)
    axes[2].legend()

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()
"""

def plot_merged_isochrone(iso):
    """
    Function to plot the merged isochrone model to visualize the transition
    between Phillips 2020 and Baraffe+15 models.

    Parameters:
    -----------
    iso : DataHolder object
        Merged isochrone returned by merge_isochrone_baraffe_phillips().
    """

    # Define colors for different sources
    color_map = {'Phillips': 'blue', 'Baraffe': 'red', 'Baraffe+Phillips': 'pink'}

    # Create figure and subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Mass vs. logT (Effective Temperature)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[0].scatter(iso.mass[mask], iso.logT[mask], label=source, color=color_map[source], alpha=0.7)
        axes[0].plot(iso.mass[mask], iso.logT[mask], color=color_map[source])  # Fixed plot function
    axes[0].set_xlabel("Mass ($M_{\odot}$)")
    axes[0].set_ylabel("$\log T_{\mathrm{eff}}$ (K)")
    axes[0].set_title("Mass vs. Effective Temperature")
    axes[0].axvline(0.07, linestyle="--", color="gray", alpha=0.5)  # Transition marker
    axes[0].axvline(0.08, linestyle="--", color="gray", alpha=0.5)  # Transition marker
    axes[0].set_xlim(0.0, 0.1)
    axes[0].legend()

    # Mass vs. logL (Luminosity)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[1].scatter(iso.mass[mask], iso.logL[mask], label=source, color=color_map[source], alpha=0.7)
        axes[1].plot(iso.mass[mask], iso.logL[mask], color=color_map[source])  # Fixed plot function
    axes[1].set_xlabel("Mass ($M_{\odot}$)")
    axes[1].set_ylabel("$\log L$ ($L_{\odot}$)")
    axes[1].set_title("Mass vs. Luminosity")
    axes[1].axvline(0.07, linestyle="--", color="gray", alpha=0.5)
    axes[1].axvline(0.08, linestyle="--", color="gray", alpha=0.5)
    axes[1].set_xlim(0.0, 0.1)
    axes[1].legend()

    # Mass vs. logg (Surface Gravity)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[2].scatter(iso.mass[mask], iso.logg[mask], label=source, color=color_map[source], alpha=0.7)
        axes[2].plot(iso.mass[mask], iso.logg[mask], color=color_map[source])  # Fixed plot function
    axes[2].set_xlabel("Mass ($M_{\odot}$)")
    axes[2].set_ylabel("$\log g$ (cm/s²)")
    axes[2].set_title("Mass vs. Surface Gravity")
    axes[2].axvline(0.07, linestyle="--", color="gray", alpha=0.5)
    axes[2].axvline(0.08, linestyle="--", color="gray", alpha=0.5)
    axes[2].set_xlim(0.0, 0.1)
    axes[2].legend()

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()

def plot_merged_bpp_isochrone(iso):
    """
    Function to plot the merged isochrone model to visualize the transition
    between Phillips 2020, Baraffe+15, and Pisa models.

    Parameters:
    -----------
    iso : DataHolder object
        Merged isochrone returned by merge_isochrone_baraffe_phillips().
    """

    # Define colors for different sources
    color_map = {'Phillips': 'blue', 
                 'Baraffe': 'red', 
                 'Pisa':'green', 
                 'Baraffe+Phillips': 'pink', 
                 'Baraffe+Pisa': 'cyan'
                }

    # Create figure and subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Mass vs. logT (Effective Temperature)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[0].scatter(iso.mass[mask], iso.logT[mask], label=source, color=color_map[source], alpha=0.7)
        axes[0].plot(iso.mass[mask], iso.logT[mask], color=color_map[source])  # Fixed plot function
    axes[0].set_xlabel("Mass ($M_{\odot}$)")
    axes[0].set_ylabel("$\log T_{\mathrm{eff}}$ (K)")
    axes[0].set_title("Mass vs. Effective Temperature")
    axes[0].axvline(0.07, linestyle="--", color="gray", alpha=0.5)  # Transition marker
    axes[0].axvline(0.08, linestyle="--", color="gray", alpha=0.5)  # Transition marker
    axes[0].legend()

    # Mass vs. logL (Luminosity)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[1].scatter(iso.mass[mask], iso.logL[mask], label=source, color=color_map[source], alpha=0.7)
        axes[1].plot(iso.mass[mask], iso.logL[mask], color=color_map[source])  # Fixed plot function
    axes[1].set_xlabel("Mass ($M_{\odot}$)")
    axes[1].set_ylabel("$\log L$ ($L_{\odot}$)")
    axes[1].set_title("Mass vs. Luminosity")
    axes[1].axvline(0.07, linestyle="--", color="gray", alpha=0.5)
    axes[1].axvline(0.08, linestyle="--", color="gray", alpha=0.5)
    axes[1].legend()

    # Mass vs. logg (Surface Gravity)
    for source in np.unique(iso.source):
        mask = iso.source == source
        axes[2].scatter(iso.mass[mask], iso.logg[mask], label=source, color=color_map[source], alpha=0.7)
        axes[2].plot(iso.mass[mask], iso.logg[mask], color=color_map[source])  # Fixed plot function
    axes[2].set_xlabel("Mass ($M_{\odot}$)")
    axes[2].set_ylabel("$\log g$ (cm/s²)")
    axes[2].set_title("Mass vs. Surface Gravity")
    axes[2].axvline(0.07, linestyle="--", color="gray", alpha=0.5)
    axes[2].axvline(0.08, linestyle="--", color="gray", alpha=0.5)
    axes[2].legend()

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()

"""def merge_isochrone_pisa_baraffe_phillips(logAge, metallicity='solar', rotation=True, iso_in=None): #?????
    ""
    Function to merge Pisa 2011, Baraffe, and Phillips 2020 models. Solar metallicity is
    Z = 0.015 for Pisa 2011 and Z = 0.015 for Phillips 2020.

    If iso_in = None, will take Pisa models to smallest available mass,
    then switch to Baraffe/Phillips. If iso_in is defined, then will take this
    isochrone to highest mass and switch to Ekstrom

    Can only handle ages at which models already exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    ""
    # Get individual Ekstrom and Pisa isochrones at desired age
    isoBaraffePhillips = merge_Ekstrom_isochrone(logAge, metallicity=metallicity,
                                       rotation=rotation)

    if iso_in == None:
        isoPisa = get_pisa_isochrone(logAge, metallicity=metallicity)
        # Create array specifying source for Pisa isochrone
        isoPisa.source = np.array(['Pisa']*len(isoPisa.mass))
    else:
        # iso_in isn't really a Pisa isochrone (it could be Baraffe-Pisa merge),
        # but name it isoPisa for simplicity anyway.
        isoPisa = iso_in
        
    # Take Pisa isochrone as high up in mass as it goes, then switch to Ekstrom.
    # Will trim Ekstrom isochrone here
    max_Pisa = max(isoPisa.mass)
    good = np.where(isoEkstrom.mass > max_Pisa)
    isoEkstrom.mass = isoEkstrom.mass[good]
    isoEkstrom.logT = isoEkstrom.logT[good]
    isoEkstrom.logg = isoEkstrom.logg[good]
    isoEkstrom.logL = isoEkstrom.logL[good]
    isoEkstrom.mass_current = isoEkstrom.mass_current[good]
    isoEkstrom.phase = isoEkstrom.phase[good]
    isoEkstrom.logT_WR = isoEkstrom.logT_WR[good]
    
    # Make array containing Ekstrom source ID
    isoEkstrom.source = np.array(['Ekstrom']*len(isoEkstrom.mass))

    # Combine the arrays
    M = np.append(isoPisa.mass, isoEkstrom.mass)
    logT = np.append(isoPisa.logT, isoEkstrom.logT)
    logg = np.append(isoPisa.logg, isoEkstrom.logg)
    logL = np.append(isoPisa.logL, isoEkstrom.logL)
    mcurr = np.append(isoPisa.mass_current, isoEkstrom.mass_current)
    phase = np.append(isoPisa.phase, isoEkstrom.phase)
    logT_WR = np.append(isoPisa.logT, isoEkstrom.logT_WR)
    source = np.append(isoPisa.source, isoEkstrom.source)

    iso = objects.DataHolder()
    iso.mass = M
    iso.logL = logL
    iso.logg = logg
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.logT_WR = logT_WR
    iso.source = source

    return iso
"""

def merge_isochrone_baraffe_pisa(logAge, metallicity='solar'):
    """
    Function to merge Baraffe+15 and Pisa 2011 models. Will take
    100% Baraffe+15 between 0.07 - 0.4 M_sun, transition between
    0.4 - 0.5 M_sun, and take 100% Pisa from 0.5 M_sun and up.

    Can only handle ages at which models already exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    """
    if metallicity != 'solar':
        print( 'Non-solar metallicity not supported yet')
        return

    # Get individual Baraffe and Pisa isochrones at desired age. Note
    # that this will also give the Baraffe models a finer mass sampling
    isoBaraffe = get_Baraffe15_isochrone(logAge, metallicity=metallicity)
    isoPisa = get_pisa_isochrone(logAge, metallicity=metallicity)

    # Identify M <= 0.4 M_sun in Baraffe and M >= 0.5 M_sun in Pisa
    good_b = np.where(isoBaraffe.mass <= 0.4)
    good_p = np.where(isoPisa.mass >= 0.5)

    # Sample between 0.4 M_sun and 0.5 M_sun in steps of 0.02 M_sun.
    # Will do linear combo of Baraffe and Pisa over this range
    mid_mass = np.arange(0.4, 0.5+0.01, 0.02)
    mid_logT = []
    mid_logL = []
    mid_logG = []
    mid_Mcurr = []
    mid_phase = []
    for mass in mid_mass:
        # Find the appropriate masses in Baraffe + Pisa to build from.
        # Baraffe has identical sampling over this range, and Pisa sampling
        # is very close. As a result, we will just take the closest mass model
        # to each mid_mass
        idx_b = np.where( abs(isoBaraffe.mass - mass)  == min(abs(isoBaraffe.mass - mass)) )
        idx_p = np.where( abs(isoPisa.mass - mass)  == min(abs(isoPisa.mass - mass)) )

        # Quality control check: we won't let the difference between model mass and
        # chosen mass to be >= 0.01 M_sun
        if ((isoPisa.mass[idx_p] - mass) >= 0.02) | ((isoBaraffe.mass[idx_b] - mass) >= 0.02):
            print( 'WARNING: Baraffe or Pisa model interpolation between 0.4 - 0.5 M_sun may \
            be inaccurate. Check this!')
            pdb.set_trace()
    
        # Now, do the linear combo of models at this mass, weighted by distance from
        # 0.4 or 0.5 (whichever is appropriate)
        diff = 0.5 - 0.4
        weight_b = (0.5 - mass) / diff
        weight_p = 1.0 - weight_b
        print( 'Baraffe {0} and Pisa {1} at mass {2}'.format(weight_b, weight_p, mass))

        # Now, do the merge IN LINEAR SPACE!
        Teff = (10**isoBaraffe.logT[idx_b] * weight_b) + \
            (10**isoPisa.logT[idx_p] * weight_p)
        L = (10**isoBaraffe.logL[idx_b] * weight_b) + \
            (10**isoPisa.logL[idx_p] * weight_p)
        g = (10**isoBaraffe.logg[idx_b] * weight_b) + \
            (10**isoPisa.logg[idx_p] * weight_p)
        mcurr = (isoBaraffe.mass_current[idx_b] * weight_b) + \
                (isoPisa.mass_current[idx_p] * weight_p)
        phase = np.round((isoBaraffe.mass_current[idx_b] * weight_b) + \
                         (isoPisa.mass_current[idx_p] * weight_p))
        
        mid_logT = np.concatenate((mid_logT, np.log10(Teff)))
        mid_logL = np.concatenate((mid_logL, np.log10(L)))
        mid_logG = np.concatenate((mid_logG, np.log10(g)))
        mid_Mcurr = np.concatenate((mid_Mcurr, mcurr))
        mid_phase = np.concatenate((mid_phase, phase))

    # Now, final isochrone will be combination of Baraffe at M<=0.4,
    # Pisa at M>=0.5, and the combination inbetween
    mass = np.concatenate((isoBaraffe.mass[good_b], mid_mass, isoPisa.mass[good_p]))
    logT = np.concatenate((isoBaraffe.logT[good_b], mid_logT, isoPisa.logT[good_p]))
    logL = np.concatenate((isoBaraffe.logL[good_b], mid_logL, isoPisa.logL[good_p]))
    logG = np.concatenate((isoBaraffe.logg[good_b], mid_logG, isoPisa.logg[good_p]))
    mcurr = np.concatenate((isoBaraffe.mass_current[good_b], mid_Mcurr, isoPisa.mass_current[good_p]))
    phase = np.concatenate((isoBaraffe.phase[good_b], mid_phase, isoPisa.phase[good_p]))

    # Also add a source flag
    source = np.concatenate( (['Baraffe']*len(good_b[0]), ['Baraffe+Pisa']*len(mid_mass),
                              ['Pisa']*len(good_p[0])) )

    iso = objects.DataHolder()
    iso.mass = mass
    iso.logL = logL
    iso.logg = logG
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.source = source

    return iso

def merge_isochrone_baraffe_pisa_phillips(logAge, metallicity='solar', rotation=True, iso_in=None):
    """
    Merges Phillips, Baraffe, and Pisa isochrones to create a complete evolutionary track.
    
    Inputs:
    logAge - Logarithmic Age
    metallicity - Metallicity (default is 'solar')
    rotation - Boolean indicating whether to include rotation in the models (default is True)
    
    Outputs:
    iso - An object containing the merged isochrone data.
    """
    # Merge Phillips with Baraffe
    print(f"Merging Phillips with Baraffe for logAge = {logAge}")
    isoBaraffePhillips = merge_isochrone_baraffe_phillips(logAge, metallicity=metallicity)
    
    # Merge Baraffe+Phillips with Pisa (depending on the logAge)
    if logAge <= 7.4:
        print(f"Merging Baraffe+Phillips with Pisa for logAge = {logAge}")
        iso = merge_isochrone_baraffe_pisa(logAge, metallicity=metallicity, iso_in=isoBaraffePhillips)
    else:
        # If logAge > 7.4, you may want to merge with different models (Parsec/Ekstrom)
        print(f"Merging Baraffe+Phillips with higher mass models for logAge = {logAge}")
        iso = merge_isochrone_baraffe_pisa(logAge, metallicity=metallicity, iso_in=isoBaraffePhillips)

    return iso


def merge_all_isochrones_phillips_baraffe_pisa_ekstrom_parsec(metallicity='solar', rotation=True, plot=False):
    """
    Make evolutionary isochrones containing a continuous distribution of 
    masses from the PMS to the MS. The models used are the following:

    PMS
    Phillips from 0.01 - 0.07 M_sun
    Baraffe+15 from 0.07 - 0.4 M_sun
    Pisa 2011 from 0.5 - top of grid (~7 M_sun)

    MS (M > Pisa 2011)
    Ekstrom+12 for logAge < 7.4
    Parsec V1.2s for logAge > 7.4

    BD stars 
    Phillips for 6.0 <= logAge <= 10.0

    metallicity = 'solar' --> Ekstrom+12 z014, Pisa2011 z015, Parsec z015, Phillips z00

    if plot = True, will make plots of merged isochrones in 'plots' directory,
    which must already exist
    
    Code is expected to be run in merged model working directory.
    """
    # Root data directory for Ekstrom+12 isochrones
    rootDirE = models_dir + 'Ekstrom2012/iso/'
    metalPart = 'z014/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return
    rotPart = 'rot/'
    if not rotation:
        rotPart = 'norot/'
    rootDirE += metalPart+rotPart

    # Root data directory for the Baraffe isochrones
    rootDirBaraffe = models_dir + 'Baraffe15/iso/'

    # Root data directory for Pisa isochrones
    rootDirPisa = models_dir + 'Pisa2011/iso/'
    metSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return
    rootDirPisa += metSuffix

    # Root data directory for Parsec isochrones
    rootDirParsec = models_dir + 'ParsecV1.2s/iso/'
    metalSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return        
    rootDirParsec += metalSuffix

    # Root data directory for Phillips isochrones
    rootDirPhillips = models_dir + 'Phillips2020/iso/'
    mSuffix = 'z00/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return
    rootDirPhillips +=mSuffix

    # Search both directories for iso_*.dat files
    isoFilesE = glob.glob(rootDirE + 'iso_*.dat')
    isoFilesB = glob.glob(rootDirBaraffe + 'iso_*.fits')
    isoFilesPi = glob.glob(rootDirPisa + 'iso_*.dat')
    isoFilesPa = glob.glob(rootDirParsec + 'iso_*')
    isoFilesPh = glob.glob(rootDirPhillips + 'iso_*.fits')

    # Output of merged isochrones
    if rotation == True:
        outDir = models_dir + 'merged/phillips_baraffe_pisa_ekstrom_parsec/{0}_rot_test/'.format(metSuffix[:-1])
        # Raise error if wanting Phillips data?
    else:
        outDir = models_dir + 'merged/phillips_baraffe_pisa_ekstrom_parsec/{0}_norot/'.format(metSuffix[:-1])
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    # Isolate the iso*.dat file names
    for ii in range(len(isoFilesE)):
        isoFilesE[ii] = isoFilesE[ii].split('/')[-1]

    for ii in range(len(isoFilesB)):
        isoFilesB[ii] = isoFilesB[ii].split('/')[-1]    

    for ii in range(len(isoFilesPi)):
        isoFilesPi[ii] = isoFilesPi[ii].split('/')[-1]

    for ii in range(len(isoFilesPa)):
        isoFilesPa[ii] = isoFilesPa[ii].split('/')[-1]

    for ii in range(len(isoFilesPh)):
        isoFilesPh[ii] = isoFilesPh[ii].split('/')[-1]

    # Loop through the Pisa isochrones, adding the MS and Baraffe models
    # as appropriate
    for ii in range(len(isoFilesPi)):
        isoFilePi = isoFilesPi[ii]

        logAgeStr = isoFilePi.replace('iso_', '').replace('.dat', '')
        logAge = float(logAgeStr)
        
        #-----PRE-MAIN SEQUENCE----#
        # Merge with the Baraffe+15 models from 0.07 - 0.4 M_sun. Includes
        # transition region between 0.4 - 0.5 M_sun in which we shift
        # from 100% Baraffe to 100% Pisa

        print( 'Merging isochrones Pisa + Baraffe + Phillips from ', isoFilePi)
        isoPMS = merge_isochrone_phillips_baraffe_pisa(logAge, metallicity=metallicity)

        #--------MAIN SEQUENCE-------#
        # Case where logAge <= 7.4, we merge with Ekstrom. Otherwise, merge
        # which parsec
        if logAge <= 7.4:
            if isoFilePi not in isoFilesE:
                print( 'Skipping isochrones from ', isoFilePi)
                print( 'PROBLEM WITH PISA OR EKSTROM')
                pdb.set_trace()

            print( 'Merging isochrones Phillips+Pisa+Ekstrom from ', isoFilePi)
            iso = merge_isochrone_phillips_pisa_Ekstrom(logAge, metallicity=metallicity,
                                               rotation=rotation, iso_in=isoPMS)
        else:
            if isoFilePi not in isoFilesPa:
                print( 'Skipping isochrones from ', isoFilePi)
                print( 'PROBLEM WITH PISA OR PARSEC')
                pdb.set_trace()
                
            print( 'Merging isochrones Phillips+Pisa+Parsec from ', isoFilePi)
            iso = merge_isochrone_phillips_pisa_parsec(logAge, metallicity=metallicity,
                                              iso_in=isoPMS)

        # Make test plot, if desired. These are put in plots directory
        if plot:
            # Make different models different colors
            phillips_ind = np.where(iso.source == 'Phillips')
            merge_pb = np.where(iso.source == 'Phillips+Baraffe')
            baraffe_ind = np.where(iso.source == 'Baraffe')
            merge_ind = np.where(iso.source == 'Baraffe+Pisa')
            pisa_ind = np.where(iso.source == 'Pisa')
            MS_ind = np.where( (iso.source == 'Ekstrom') | (iso.source == 'Parsec'))
            #Extract age
            logAge = isoFilePi.split('_')[1][:-4]

            # Temporarily turn off interactive plot. Turn back on at end
            py.ioff()
            
            py.figure(1)
            py.clf()
            py.plot(iso.logT[phillips_ind], iso.logL[phillips_ind], 'c-', label = 'Phillips',
                    linewidth=2)
            py.plot(iso.logT[merge_pb], iso.logL[merge_pb], 'y-', label = 'Phillips+Baraffe',
                    linewidth=2)
            py.plot(iso.logT[baraffe_ind], iso.logL[baraffe_ind], 'k-', label = 'Baraffe',
                    linewidth=2)
            py.plot(iso.logT[merge_ind], iso.logL[merge_ind], 'r-', label = 'Baraffe+Pisa',
                    linewidth=2)
            py.plot(iso.logT[pisa_ind], iso.logL[pisa_ind], 'g-', label = 'Pisa',
                    linewidth=2)
            py.plot(iso.logT[MS_ind], iso.logL[MS_ind], 'b-',
                    label = 'Ekstrom/Parsec', linewidth=2)
            py.xlabel('log Teff')
            py.ylabel('log L')
            py.title('Log Age = {0:s}'.format(logAge))
            py.legend(loc=3)
            py.axis([4.9, 3.2, -3.0, 8])
            py.savefig(outDir+'plots/iso_'+logAge+'.png')

            py.ion()
            
            
        _out = open(outDir + isoFilePi, 'w')

        hdr_fmt = '%12s  %10s  %10s  %10s %10s %12s %5s %-10s\n'
        _out.write(hdr_fmt % 
                   ('# M_init', 'log T', 'log L', 'log g', 'logT_WR', 'M_curr', 'phase', 'Source'))
        _out.write(hdr_fmt % 
                   ('# (Msun)', '(Kelvin)', '(Lsun)', '(cgs)', '(Kelvin)', '(Msun)', '()', '()'))

        for kk in range(len(iso.mass)):
            _out.write('%12.6f  %10.4f  %10.4f  %10.4f  %10.4f %12.6f %5d %-10s\n' %
                       (iso.mass[kk], iso.logT[kk], iso.logL[kk],
                        iso.logg[kk], iso.logT_WR[kk],
                        iso.mass_current[kk], iso.phase[kk], iso.source[kk]))

        _out.close()
    return

def get_pisa_isochrone(logAge, metallicity='solar'):
    """
    Load mass, effective temperature, log gravity, and log luminosity
    for the Pisa isochrones at given logAge. Code will quit if that
    logAge value doesn't exist (can make some sort of interpolation thing
    later).

    Note: mass is currently initial mass, not instantaneous mass
    
    Inputs:
    logAge - Logarithmic Age
    metallicity - in Z (def = solar of 0.014)
    """
    rootDir = models_dir + 'Pisa2011/iso/'
    metSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar Pisa 2011 metallicities not supported yet')
        return
    rootDir += metSuffix

    # Check to see if isochrone exists
    isoFile = rootDir + 'iso_%.2f.dat' % logAge
    if not os.path.exists(isoFile):
        print( 'Pisa isochrone for logAge = {0:3.2f} does\'t exist'.format(logAge))
        print( 'Quitting')
        return
        
    data = Table.read(isoFile, format='ascii')
    cols = data.keys()
    mass = data[cols[2]] #Note: this is initial mass, in M_sun
    logT = data[cols[1]] # K
    logL = data[cols[0]] # L_sun
    logg = data[cols[3]]

    # Hack: Make mass_current and phase
    mass_current = np.array(mass)
    phase = np.ones(len(mass), dtype=int)
    
    obj = objects.DataHolder()
    obj.mass = mass
    obj.logT = logT
    obj.logg = logg
    obj.logL = logL
    obj.mass_current = mass_current
    obj.phase = phase

    return obj

def get_phillips_pisa_isochrone(logAge, metallicity='solar',
                                warn_interpolated=True):
    """
    Load mass, Teff, log g, and log L for Phillips/Pisa isochrones at logAge.

    Quits if that logAge does not exist. Mass is initial mass, not
    instantaneous mass.

    Parameters
    ----------
    logAge : float
        Logarithmic age.
    metallicity : str, optional
        Metallicity label. Only ``'solar'`` (Z = 0.014) is supported.
    warn_interpolated : bool, optional
        If True, print a warning when any isochrone points come from an
        interpolated region. Default is True.

    Returns
    -------
    obj : `objects.DataHolder` or None
        Isochrone with mass, logT, logg, logL, mass_current, phase, and
        interpolated arrays. None if the file is missing or metallicity
        is unsupported.
    """
    rootDir = models_dir + 'merged/phillips_pisa'
    metSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar Phillips 2020 & Pisa 2011 metallicities not supported yet')
        return None
    rootDir += metSuffix

    # Check to see if isochrone exists
    isoFile = rootDir + 'iso_%.2f.fits' % logAge
    if not os.path.exists(isoFile):
        print( 'Phillips/Pisa isochrone for logAge = {0:3.2f} does\'t exist'.format(logAge))
        print( 'Quitting')
        return None

    data = Table.read(isoFile, format='fits')
    cols = data.keys()
    mass = data[cols[0]] #Note: this is initial mass, in M_sun
    logT = data[cols[2]] # K
    logL = data[cols[1]] # L_sun
    logg = data[cols[3]]
    interpolated = data['interpolated']

    # warning if data point is coming from interpolated region
    if warn_interpolated and np.any(interpolated):
        print( f"Warning: {np.sum(interpolated)} of {len(interpolated)} data points are interpolated.")

    # Hack: Make mass_current and phase
    mass_current = np.array(mass)
    phase = np.ones(len(mass), dtype=int)
    
    obj = objects.DataHolder()
    obj.mass = mass
    obj.logT = logT
    obj.logg = logg
    obj.logL = logL
    obj.mass_current = mass_current
    obj.phase = phase
    obj.interpolated = interpolated

    return obj

def merge_isochrone_baraffe_pisa(logAge, metallicity='solar'):
    """
    Function to merge Baraffe+15 and Pisa 2011 models. Will take
    100% Baraffe+15 between 0.07 - 0.4 M_sun, transition between
    0.4 - 0.5 M_sun, and take 100% Pisa from 0.5 M_sun and up.

    Can only handle ages at which models already exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    """
    if metallicity != 'solar':
        print( 'Non-solar metallicity not supported yet')
        return

    # Get individual Baraffe and Pisa isochrones at desired age. Note
    # that this will also give the Baraffe models a finer mass sampling
    isoBaraffe = get_Baraffe15_isochrone(logAge, metallicity=metallicity)
    isoPisa = get_pisa_isochrone(logAge, metallicity=metallicity)

    # Identify M <= 0.4 M_sun in Baraffe and M >= 0.5 M_sun in Pisa
    good_b = np.where(isoBaraffe.mass <= 0.4)
    good_p = np.where(isoPisa.mass >= 0.5)

    # Sample between 0.4 M_sun and 0.5 M_sun in steps of 0.02 M_sun.
    # Will do linear combo of Baraffe and Pisa over this range
    mid_mass = np.arange(0.4, 0.5+0.01, 0.02)
    mid_logT = []
    mid_logL = []
    mid_logG = []
    mid_Mcurr = []
    mid_phase = []
    for mass in mid_mass:
        # Find the appropriate masses in Baraffe + Pisa to build from.
        # Baraffe has identical sampling over this range, and Pisa sampling
        # is very close. As a result, we will just take the closest mass model
        # to each mid_mass
        idx_b = np.where( abs(isoBaraffe.mass - mass)  == min(abs(isoBaraffe.mass - mass)) )
        idx_p = np.where( abs(isoPisa.mass - mass)  == min(abs(isoPisa.mass - mass)) )

        # Quality control check: we won't let the difference between model mass and
        # chosen mass to be >= 0.01 M_sun
        if ((isoPisa.mass[idx_p] - mass) >= 0.02) | ((isoBaraffe.mass[idx_b] - mass) >= 0.02):
            print( 'WARNING: Baraffe or Pisa model interpolation between 0.4 - 0.5 M_sun may \
            be inaccurate. Check this!')
            pdb.set_trace()
    
        # Now, do the linear combo of models at this mass, weighted by distance from
        # 0.4 or 0.5 (whichever is appropriate)
        diff = 0.5 - 0.4
        weight_b = (0.5 - mass) / diff
        weight_p = 1.0 - weight_b
        print( 'Baraffe {0} and Pisa {1} at mass {2}'.format(weight_b, weight_p, mass))

        # Now, do the merge IN LINEAR SPACE!
        Teff = (10**isoBaraffe.logT[idx_b] * weight_b) + \
            (10**isoPisa.logT[idx_p] * weight_p)
        L = (10**isoBaraffe.logL[idx_b] * weight_b) + \
            (10**isoPisa.logL[idx_p] * weight_p)
        g = (10**isoBaraffe.logg[idx_b] * weight_b) + \
            (10**isoPisa.logg[idx_p] * weight_p)
        mcurr = (isoBaraffe.mass_current[idx_b] * weight_b) + \
                (isoPisa.mass_current[idx_p] * weight_p)
        phase = np.round((isoBaraffe.mass_current[idx_b] * weight_b) + \
                         (isoPisa.mass_current[idx_p] * weight_p))
        
        mid_logT = np.concatenate((mid_logT, np.log10(Teff)))
        mid_logL = np.concatenate((mid_logL, np.log10(L)))
        mid_logG = np.concatenate((mid_logG, np.log10(g)))
        mid_Mcurr = np.concatenate((mid_Mcurr, mcurr))
        mid_phase = np.concatenate((mid_phase, phase))

    # Now, final isochrone will be combination of Baraffe at M<=0.4,
    # Pisa at M>=0.5, and the combination inbetween
    mass = np.concatenate((isoBaraffe.mass[good_b], mid_mass, isoPisa.mass[good_p]))
    logT = np.concatenate((isoBaraffe.logT[good_b], mid_logT, isoPisa.logT[good_p]))
    logL = np.concatenate((isoBaraffe.logL[good_b], mid_logL, isoPisa.logL[good_p]))
    logG = np.concatenate((isoBaraffe.logg[good_b], mid_logG, isoPisa.logg[good_p]))
    mcurr = np.concatenate((isoBaraffe.mass_current[good_b], mid_Mcurr, isoPisa.mass_current[good_p]))
    phase = np.concatenate((isoBaraffe.phase[good_b], mid_phase, isoPisa.phase[good_p]))

    # Also add a source flag
    source = np.concatenate( (['Baraffe']*len(good_b[0]), ['Baraffe+Pisa']*len(mid_mass),
                              ['Pisa']*len(good_p[0])) )

    iso = objects.DataHolder()
    iso.mass = mass
    iso.logL = logL
    iso.logg = logG
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.source = source

    return iso
    
def merge_isochrone_phillips_baraffe_pisa(logAge, metallicity='solar'):
    """
    Function to merge Phillips 2020, Baraffe+15, and Pisa 2011 models.
    
    - Uses Phillips for M <= 0.07 M_sun
    - Transitions from Phillips to Baraffe between 0.07 - 0.075 M_sun
    - Uses Baraffe for 0.075 - 0.4 M_sun
    - Transitions from Baraffe to Pisa between 0.4 - 0.5 M_sun
    - Uses Pisa for M >= 0.5 M_sun
    
    Can only handle ages where models exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    """
    if metallicity != 'solar':
        print('Non-solar metallicity not supported yet')
        return

    # Load individual isochrones
    isoPhillips = get_phillips_isochrone(logAge, metallicity=metallicity)
    isoBaraffe = get_Baraffe15_isochrone(logAge, metallicity=metallicity)
    isoPisa = get_pisa_isochrone(logAge, metallicity=metallicity)

    # Identify mass ranges
    good_p = np.where(isoPhillips.mass <= 0.07)
    good_b = np.where((isoBaraffe.mass >= 0.075) & (isoBaraffe.mass <= 0.4))
    good_pisa = np.where(isoPisa.mass >= 0.5)

    # Transition Phillips → Baraffe (0.07 - 0.075 M_sun)
    mid_mass_phillips_baraffe = np.arange(0.07, 0.075, 0.002)
    mid_logT_phillips_baraffe, mid_logL_phillips_baraffe, mid_logG_phillips_baraffe = [], [], []
    mid_Mcurr_phillips_baraffe, mid_phase_phillips_baraffe = [], []
    
    for mass in mid_mass_phillips_baraffe:
        idx_p = np.argmin(np.abs(isoPhillips.mass - mass))
        idx_b = np.argmin(np.abs(isoBaraffe.mass - mass))
        
        weight_p = (0.075 - mass) / 0.005
        weight_b = 1.0 - weight_p

        Teff = (10**isoPhillips.logT[idx_p] * weight_p) + (10**isoBaraffe.logT[idx_b] * weight_b)
        L = (10**isoPhillips.logL[idx_p] * weight_p) + (10**isoBaraffe.logL[idx_b] * weight_b)
        g = (10**isoPhillips.logg[idx_p] * weight_p) + (10**isoBaraffe.logg[idx_b] * weight_b)
        mcurr = (isoPhillips.mass_current[idx_p] * weight_p) + (isoBaraffe.mass_current[idx_b] * weight_b)
        phase = np.round((isoPhillips.phase[idx_p] * weight_p) + (isoBaraffe.phase[idx_b] * weight_b))

        mid_logT_phillips_baraffe.append(np.log10(Teff))
        mid_logL_phillips_baraffe.append(np.log10(L))
        mid_logG_phillips_baraffe.append(np.log10(g))
        mid_Mcurr_phillips_baraffe.append(mcurr)
        mid_phase_phillips_baraffe.append(phase)

    # Transition Baraffe → Pisa (0.4 - 0.5 M_sun)
    mid_mass_baraffe_pisa = np.arange(0.4, 0.5 + 0.01, 0.02)
    mid_logT_baraffe_pisa, mid_logL_baraffe_pisa, mid_logG_baraffe_pisa = [], [], []
    mid_Mcurr_baraffe_pisa, mid_phase_baraffe_pisa = [], []

    for mass in mid_mass_baraffe_pisa:
        idx_b = np.argmin(np.abs(isoBaraffe.mass - mass))
        idx_p = np.argmin(np.abs(isoPisa.mass - mass))
        
        weight_b = (0.5 - mass) / 0.1
        weight_p = 1.0 - weight_b

        Teff = (10**isoBaraffe.logT[idx_b] * weight_b) + (10**isoPisa.logT[idx_p] * weight_p)
        L = (10**isoBaraffe.logL[idx_b] * weight_b) + (10**isoPisa.logL[idx_p] * weight_p)
        g = (10**isoBaraffe.logg[idx_b] * weight_b) + (10**isoPisa.logg[idx_p] * weight_p)
        mcurr = (isoBaraffe.mass_current[idx_b] * weight_b) + (isoPisa.mass_current[idx_p] * weight_p)
        phase = np.round((isoBaraffe.phase[idx_b] * weight_b) + (isoPisa.phase[idx_p] * weight_p))

        mid_logT_baraffe_pisa.append(np.log10(Teff))
        mid_logL_baraffe_pisa.append(np.log10(L))
        mid_logG_baraffe_pisa.append(np.log10(g))
        mid_Mcurr_baraffe_pisa.append(mcurr)
        mid_phase_baraffe_pisa.append(phase)

    # Merge all mass regimes
    mass = np.concatenate((isoPhillips.mass[good_p], mid_mass_phillips_baraffe, isoBaraffe.mass[good_b], 
                           mid_mass_baraffe_pisa, isoPisa.mass[good_pisa]))
    
    logT = np.concatenate((isoPhillips.logT[good_p], mid_logT_phillips_baraffe, isoBaraffe.logT[good_b], 
                           mid_logT_baraffe_pisa, isoPisa.logT[good_pisa]))
    
    logL = np.concatenate((isoPhillips.logL[good_p], mid_logL_phillips_baraffe, isoBaraffe.logL[good_b], 
                           mid_logL_baraffe_pisa, isoPisa.logL[good_pisa]))
    
    logG = np.concatenate((isoPhillips.logg[good_p], mid_logG_phillips_baraffe, isoBaraffe.logg[good_b], 
                           mid_logG_baraffe_pisa, isoPisa.logg[good_pisa]))
    
    mcurr = np.concatenate((isoPhillips.mass_current[good_p], 
                            mid_Mcurr_phillips_baraffe,
                            isoBaraffe.mass_current[good_b], 
                            mid_Mcurr_baraffe_pisa, 
                            isoPisa.mass_current[good_pisa]))
    
    phase = np.concatenate((isoPhillips.phase[good_p], mid_phase_phillips_baraffe, isoBaraffe.phase[good_b], 
                            mid_phase_baraffe_pisa, isoPisa.phase[good_pisa]))

    # Source flag
    source = np.concatenate((['Phillips']*len(good_p[0]), ['Phillips+Baraffe']*len(mid_mass_phillips_baraffe),
                             ['Baraffe']*len(good_b[0]), ['Baraffe+Pisa']*len(mid_mass_baraffe_pisa),
                             ['Pisa']*len(good_pisa[0])))

    # Create final merged isochrone object
    iso = objects.DataHolder()
    iso.mass = mass
    iso.logL = logL
    iso.logg = logG
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.source = source

    return iso

def merge_isochrone_pisa_Ekstrom(logAge, metallicity='solar', rotation=True, iso_in=None):
    """
    Function to merge Pisa 2011 and Ekstrom 2012 models. Solar metallicity is
    Z = 0.015 for Pisa 2011 and Z = 0.014 for Ekstrom+12.

    If iso_in = None, will take Pisa models to highest available mass,
    then switch to Ekstrom. If iso_in is defined, then will take this
    isochrone to highest mass and switch to Ekstrom

    Can only handle ages at which models already exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    """
    # Get individual Ekstrom and Pisa isochrones at desired age
    isoEkstrom = get_Ekstrom_isochrone(logAge, metallicity=metallicity,
                                       rotation=rotation)

    if iso_in == None:
        isoPisa = get_pisa_isochrone(logAge, metallicity=metallicity)
        # Create array specifying source for Pisa isochrone
        isoPisa.source = np.array(['Pisa']*len(isoPisa.mass))
    else:
        # iso_in isn't really a Pisa isochrone (it could be Baraffe-Pisa merge),
        # but name it isoPisa for simplicity anyway.
        isoPisa = iso_in
        
    # Take Pisa isochrone as high up in mass as it goes, then switch to Ekstrom.
    # Will trim Ekstrom isochrone here
    max_Pisa = max(isoPisa.mass)
    good = np.where(isoEkstrom.mass > max_Pisa)
    isoEkstrom.mass = isoEkstrom.mass[good]
    isoEkstrom.logT = isoEkstrom.logT[good]
    isoEkstrom.logg = isoEkstrom.logg[good]
    isoEkstrom.logL = isoEkstrom.logL[good]
    isoEkstrom.mass_current = isoEkstrom.mass_current[good]
    isoEkstrom.phase = isoEkstrom.phase[good]
    isoEkstrom.logT_WR = isoEkstrom.logT_WR[good]
    
    # Make array containing Ekstrom source ID
    isoEkstrom.source = np.array(['Ekstrom']*len(isoEkstrom.mass))

    # Combine the arrays
    M = np.append(isoPisa.mass, isoEkstrom.mass)
    logT = np.append(isoPisa.logT, isoEkstrom.logT)
    logg = np.append(isoPisa.logg, isoEkstrom.logg)
    logL = np.append(isoPisa.logL, isoEkstrom.logL)
    mcurr = np.append(isoPisa.mass_current, isoEkstrom.mass_current)
    phase = np.append(isoPisa.phase, isoEkstrom.phase)
    logT_WR = np.append(isoPisa.logT, isoEkstrom.logT_WR)
    source = np.append(isoPisa.source, isoEkstrom.source)

    iso = objects.DataHolder()
    iso.mass = M
    iso.logL = logL
    iso.logg = logg
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.logT_WR = logT_WR
    iso.source = source

    return iso

def merge_isochrone_phillips_pisa_Ekstrom(logAge, metallicity='solar', rotation=True, iso_in=None):
    """
    Function to merge Phillips 2020, Pisa 2011, and Ekstrom 2012 models. Solar metallicity is
    Z = 0.015 for Pisa 2011 and Phillips 2020 and Z = 0.014 for Ekstrom+12.

    If iso_in = None, will take Pisa/Phillips models to highest available mass,
    then switch to Ekstrom. If iso_in is defined, then will take this
    isochrone to highest mass and switch to Ekstrom

    Can only handle ages at which models already exist:
    logAge = 6.0 - 8.0, delta logAge = 0.01
    """
    # Get individual Ekstrom and Pisa isochrones at desired age
    isoEkstrom = get_Ekstrom_isochrone(logAge, metallicity=metallicity,
                                       rotation=rotation)

    if iso_in == None:
        isoPhPisa = get_phillips_pisa_isochrone(logAge, metallicity=metallicity)
        # Create array specifying source for Pisa isochrone
        isoPhPisa.source = np.array(['Phillips/Pisa']*len(isoPhPisa.mass))
    else:
        # iso_in isn't really a Pisa isochrone (it could be Baraffe-Pisa merge),
        # but name it isoPisa for simplicity anyway.
        isoPhPisa = iso_in
        
    # Take Pisa isochrone as high up in mass as it goes, then switch to Ekstrom.
    # Will trim Ekstrom isochrone here
    max_PhPisa = max(isoPhPisa.mass)
    good = np.where(isoEkstrom.mass > max_PhPisa)
    isoEkstrom.mass = isoEkstrom.mass[good]
    isoEkstrom.logT = isoEkstrom.logT[good]
    isoEkstrom.logg = isoEkstrom.logg[good]
    isoEkstrom.logL = isoEkstrom.logL[good]
    isoEkstrom.mass_current = isoEkstrom.mass_current[good]
    isoEkstrom.phase = isoEkstrom.phase[good]
    isoEkstrom.logT_WR = isoEkstrom.logT_WR[good]
    
    # Make array containing Ekstrom source ID
    isoEkstrom.source = np.array(['Ekstrom']*len(isoEkstrom.mass))

    # Combine the arrays
    M = np.append(isoPhPisa.mass, isoEkstrom.mass)
    logT = np.append(isoPhPisa.logT, isoEkstrom.logT)
    logg = np.append(isoPhPisa.logg, isoEkstrom.logg)
    logL = np.append(isoPhPisa.logL, isoEkstrom.logL)
    mcurr = np.append(isoPhPisa.mass_current, isoEkstrom.mass_current)
    phase = np.append(isoPhPisa.phase, isoEkstrom.phase)
    logT_WR = np.append(isoPhPisa.logT, isoEkstrom.logT_WR)
    source = np.append(isoPhPisa.source, isoEkstrom.source)

    iso = objects.DataHolder()
    iso.mass = M
    iso.logL = logL
    iso.logg = logg
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.logT_WR = logT_WR
    iso.source = source

    # add interpolation flag
    if hasattr(isoPhPisa, 'interpolated'):
        interpolated = np.append(isoPhPisa.interpolated, np.full(len(isoEkstrom.mass), False))
        iso.interpolated = interpolated
        
        # fix phase values in interpolated region (assign same as min-mass Parsec object)
        min_ekstrom_mass_idx = np.argmin(isoEkstrom.mass)
        min_ekstrom_phase = isoEkstrom.phase[min_ekstrom_mass_idx]
        interp_mask = iso.interpolated == True
        iso.phase[interp_mask] = min_ekstrom_phase

    return iso

def get_Ekstrom_isochrone(logAge, metallicity='solar', rotation=True):
    """
    Load mass, effective temperature, log gravity, and log luminosity
    for the Ekstrom isochrones at given logAge. Code will quit if that
    logAge value doesn't exist (can make some sort of interpolation thing
    later). Also interpolate model to finer mass grid

    Note: mass is currently initial mass, not instantaneous mass
    
    Inputs:
    logAge - Logarithmic Age
    metallicity - in Z (def = solar of 0.014)
    """
    rootDir = models_dir + 'Ekstrom2012/iso/'
    metSuffix = 'z014/'
    if metallicity != 'solar':
        print( 'Non-solar Ekstrom+12 metallicities not supported yet')
        return
    rotSuffix = 'rot/'
    if not rotation:
        rotSuffix = 'norot/'
        
    rootDir += metSuffix + rotSuffix

    # Check to see if isochrone exists
    isoFile = rootDir + 'iso_%.2f.dat' % logAge
    if not os.path.exists(isoFile):
        print( 'Ekstrom isochrone for logAge = {0:3.2f} does\'t exist'.format(logAge))
        print( 'Quitting')
        return
        
    data = Table.read(isoFile, format='ascii')
    cols = data.keys()
    mass = data[cols[2]] #Note: this is initial mass, in M_sun
    logT = data[cols[7]] # K
    logL = data[cols[6]] # L_sun
    logT_WR = data[cols[8]] # K; if this doesn't equal logT, we have a WR star
    mass_curr = data[cols[3]]
    phase = np.ones(len(data), dtype=int)

    # Need to calculate log g from mass and R
    R_sun = 7.*10**10 #cm
    M_sun = 2.*10**33 #g
    G_const = 6.67*10**-8 #cgs
    
    radius = data[cols[19]] #R_sun
    logg = np.log10( (G_const * np.array(mass).astype(float) * M_sun) /
                     (np.array(radius).astype(float) * R_sun)**2 )
    
    # Interpolate isochrone to finer mass grid on main-ish sequence
    # (1-60 M_sun, or the highest mass in the model); don't want to
    # completely redo all sampling, just this region
    if max(mass) > 60:
        new_masses = np.arange(1, 60+0.1, 0.5)
    else:
        new_masses = np.arange(1, max(mass), 0.5)
    mass_grid = np.append(new_masses, mass)
    mass_grid.sort() # Make sure grid is in proper order

    # Build interpolators in linear space
    f_logT = interpolate.interp1d(mass, 10**logT, kind='linear')
    f_logL = interpolate.interp1d(mass, 10**logL, kind='linear')
    f_logT_WR = interpolate.interp1d(mass, 10**logT_WR, kind='linear')
    f_logg = interpolate.interp1d(mass, 10**logg, kind='linear')
    f_Mcurr = interpolate.interp1d(mass, mass_curr, kind='linear')
    f_phase = interpolate.interp1d(mass, phase, kind='linear')

    # Do interpolation, convert back to logspace
    logT_interp = np.log10(f_logT(mass_grid))
    logL_interp = np.log10(f_logL(mass_grid))
    logT_WR_interp = np.log10(f_logT_WR(mass_grid))
    logg_interp = np.log10(f_logg(mass_grid))
    Mcurr_interp = f_Mcurr(mass_grid)
    phase_interp = f_phase(mass_grid)
    
    # Make isochrone
    obj = objects.DataHolder()
    obj.mass = mass_grid
    obj.logT = logT_interp
    obj.logg = logg_interp
    obj.logL = logL_interp
    obj.logT_WR = logT_WR_interp
    obj.mass_current = Mcurr_interp
    obj.phase = phase_interp

    return obj

def merge_isochrone_pisa_parsec(logAge, metallicity='solar', iso_in=None):
    """
    Function to merge Pisa 2011 and ParsecV1.2s models. Solar metallicity is
    Z = 0.015.

    If iso_in = None, will take Pisa models to highest available mass,
    then switch to Parsec. If iso_in is defined, then will take this
    isochrone to highest mass and switch to Parsec

    Can only handle ages at which both sets of models already exist:
    logAge = 6.6 - 8.0, delta logAge = 0.01
    """
    isoParsec = get_parsec_isochrone(logAge, metallicity=metallicity)
    # Make Parsec source array
    isoParsec.source = np.array(['Parsec']*len(isoParsec.mass))

    # Define isoPisa based on iso_in input
    if iso_in != None:
        isoPisa = iso_in
    else:
        isoPisa = get_pisa_isochrone(logAge, metallicity=metallicity)
        isoPisa.source = np.array(['Pisa']*len(isoPisa.mass))
        
    # Use Pisa model as high up as it goes, then switch to Parsec
    max_Pisa = max(isoPisa.mass)
    good = np.where(isoParsec.mass > max_Pisa)
    isoParsec.mass = isoParsec.mass[good]
    isoParsec.logT = isoParsec.logT[good]
    isoParsec.logg = isoParsec.logg[good]
    isoParsec.logL = isoParsec.logL[good]
    isoParsec.mass_current = isoParsec.mass_current[good]
    isoParsec.phase = isoParsec.phase[good]
    isoParsec.source = isoParsec.source[good]
    
    # Combine the arrays
    M = np.append(isoPisa.mass, isoParsec.mass)
    logT = np.append(isoPisa.logT, isoParsec.logT)
    logg = np.append(isoPisa.logg, isoParsec.logg)
    logL = np.append(isoPisa.logL, isoParsec.logL)
    mcurr = np.append(isoPisa.mass_current, isoParsec.mass_current)
    phase = np.append(isoPisa.phase, isoParsec.phase)
    logT_WR = np.append(isoPisa.logT, isoParsec.logT)
    source = np.append(isoPisa.source, isoParsec.source)

    iso = objects.DataHolder()
    iso.mass = M
    iso.logL = logL
    iso.logg = logg
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.logT_WR = logT_WR
    iso.source = source

    return iso

def merge_isochrone_phillips_pisa_parsec(logAge, metallicity='solar', iso_in=None):
    """
    Function to merge Phillips 2020, Pisa 2011,  and ParsecV1.2s models. 
    Solar metallicity is Z = 0.015.

    If iso_in = None, will take Phillips/Pisa models to highest available mass,
    then switch to Parsec. If iso_in is defined, then will take this
    isochrone to highest mass and switch to Parsec

    Can only handle ages at which both sets of models already exist:
    logAge = 6.6 - 8.0, delta logAge = 0.01
    """
    isoParsec = get_parsec_isochrone(logAge, metallicity=metallicity)
    # Make Parsec source array
    isoParsec.source = np.array(['Parsec']*len(isoParsec.mass))

    # Define isoPhPisa based on iso_in input
    if iso_in != None:
        isoPhPisa = iso_in
    else:
        isoPhPisa = get_phillips_pisa_isochrone(logAge, metallicity=metallicity)
        isoPhPisa.source = np.array(['Phillips/Pisa']*len(isoPhPisa.mass))
        
    # Use Phillips/Pisa model as high up as it goes, then switch to Parsec
    max_PhPisa = max(isoPhPisa.mass)
    good = np.where(isoParsec.mass > max_PhPisa)
    isoParsec.mass = isoParsec.mass[good]
    isoParsec.logT = isoParsec.logT[good]
    isoParsec.logg = isoParsec.logg[good]
    isoParsec.logL = isoParsec.logL[good]
    isoParsec.mass_current = isoParsec.mass_current[good]
    isoParsec.phase = isoParsec.phase[good]
    isoParsec.source = isoParsec.source[good]
    
    
    # Combine the arrays
    M = np.append(isoPhPisa.mass, isoParsec.mass)
    logT = np.append(isoPhPisa.logT, isoParsec.logT)
    logg = np.append(isoPhPisa.logg, isoParsec.logg)
    logL = np.append(isoPhPisa.logL, isoParsec.logL)
    mcurr = np.append(isoPhPisa.mass_current, isoParsec.mass_current)
    phase = np.append(isoPhPisa.phase, isoParsec.phase)
    logT_WR = np.append(isoPhPisa.logT, isoParsec.logT)
    source = np.append(isoPhPisa.source, isoParsec.source)

    iso = objects.DataHolder()
    iso.mass = M
    iso.logL = logL
    iso.logg = logg
    iso.logT = logT
    iso.mass_current = mcurr
    iso.phase = phase
    iso.logT_WR = logT_WR
    iso.source = source

    # add interpolation flag
    if hasattr(isoPhPisa, 'interpolated'):
        interpolated = np.append(isoPhPisa.interpolated, np.full(len(isoEkstrom.mass), False))
        iso.interpolated = interpolated

        # fix phase values in interpolated region (assign same as min-mass Parsec object)
        min_parsec_mass_idx = np.argmin(isoParsec.mass)
        min_parsec_phase = isoParsec.phase[min_parsec_mass_idx]
        interp_mask = iso.interpolated == True
        iso.phase[interp_mask] = min_parsec_phase

    return iso

def make_parsec_iso(logAge, metallicity='solar'):
    """
    Make parsec isochrone in Popstar iso object
    """
    isoParsec = get_parsec_isochrone(logAge, metallicity=metallicity)
    isoParsec.source = np.array(['Parsec']*len(isoParsec.mass))

    iso = objects.DataHolder()
    iso.mass = isoParsec.mass
    iso.logL = isoParsec.logL
    iso.logg = isoParsec.logg
    iso.logT = isoParsec.logT
    iso.mass_current = isoParsec.mass_current
    iso.phase = isoParsec.phase
    iso.logT_WR = isoParsec.logT
    iso.source = isoParsec.source

    return iso
def get_phillips_parsec_isochrone(logAge, metallicity='solar'):
    """
    Load mass, effective temperature, log gravity, and log luminosity
    for the Phillips/Parsec isochrones at given logAge. Code will quit if that
    logAge value doesn't exist (can make some sort of interpolation thing
    later).

    Note: mass is currently initial mass, not instantaneous mass
    
    Inputs:
    logAge - Logarithmic Age
    metallicity - in Z (def = solar of 0.014)
    """
    rootDir = models_dir + 'merged/phillips_parsec/'
    print(rootDir)
    metSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar Phillips 2020/Parsec 2011 metallicities not supported yet')
        return
    rootDir += metSuffix
    print(rootDir)

    # Check to see if isochrone exists
    #isoFile = rootDir + 'iso_%.2f.fits' % logAge
    isoFile = os.path.join(rootDir, f'iso_{logAge:.2f}.fits') # not a valid path!
    print(isoFile)
    print(os.path.exists(isoFile))
    if not os.path.exists(isoFile):
        print( 'Phillips/Parsec isochrone for logAge = {0:3.2f} does\'t exist'.format(logAge))   # formatting in this is mismatched!
        print( 'Quitting')  #raise IO Error instead and merge with top print statement
        return
        
    data = Table.read(isoFile, format='fits')
    cols = data.keys()
    mass = data[cols[0]] #Note: this is initial mass, in M_sun
    logT = data[cols[2]] # K
    logL = data[cols[1]] # L_sun
    logg = data[cols[3]]
    # current mass
    Mcurr = data[cols[3]] #issue!
    phase = np.ones(len(data), dtype=int)
    
    obj = objects.DataHolder()
    obj.mass = mass
    obj.logT = logT
    obj.logg = logg
    obj.logL = logL
    obj.mass_current = Mcurr
    obj.phase = phase
    
    return obj

def make_phillips_parsec_iso(logAge, metallicity='solar'):
    """
    Make Phillips/Parsec isochrone in Popstar iso object
    """
    isoPhParsec = get_phillips_parsec_isochrone(logAge, metallicity=metallicity)
    isoPhParsec.source = np.array(['Phillips+Parsec']*len(isoPhParsec.mass))

    iso = objects.DataHolder()
    iso.mass = isoPhParsec.mass
    iso.logL = isoPhParsec.logL
    iso.logg = isoPhParsec.logg
    iso.logT = isoPhParsec.logT
    iso.mass_current = isoPhParsec.mass_current
    iso.phase = isoPhParsec.phase
    iso.logT_WR = isoPhParsec.logT
    iso.source = isoPhParsec.source

    return iso

def get_parsec_isochrone(logAge, metallicity='solar'):
    """
    Load mass, effective temperature, log gravity, and log luminosity
    for the Parsec isochrones at given logAge. Code will quit if that
    logAge value doesn't exist (can make some sort of interpolation thing
    later).

    Note: mass is currently initial mass, not instantaneous mass
    
    Inputs:
    logAge - Logarithmic Age
    metallicity - in Z (def = solar of 0.014)
    """
    rootDir = models_dir + 'ParsecV1.2s/iso/'
    metSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar Parsec 2011 metallicities not supported yet')
        return
    rootDir += metSuffix

    # Check to see if isochrone exists
    isoFile = rootDir + 'iso_%.2f.dat' % logAge
    if not os.path.exists(isoFile):
        print( 'Parsec isochrone for logAge = {0:3.2f} does\'t exist'.format(logAge))
        print( 'Quitting')
        return
        
    data = Table.read(isoFile, format='ascii')
    cols = data.keys()
    mass = data[cols[2]] #Note: this is initial mass, in M_sun
    logT = data[cols[5]] # K
    logL = data[cols[4]] # L_sun
    logg = data[cols[6]]
    Mcurr = data[cols[3]]
    phase = np.ones(len(data), dtype=int)
    
    obj = objects.DataHolder()
    obj.mass = mass
    obj.logT = logT
    obj.logg = logg
    obj.logL = logL
    obj.mass_current = Mcurr
    obj.phase = phase
    
    return obj

#### CREATING FINAL ISOCHRONES ####

def merge_all_isochrones_phillips_baraffe_pisa_ekstrom_parsec2(logAge_arr=np.arange(6.0, 10.1, 0.01),
                                                               metallicity='solar',
                                                               rotation=True, plot=False):
    """
    Make evolutionary isochrones containing a continuous distribution of 
    masses from the PMS to the MS. The models used are the following:

    PMS (logAge < 8)
    Phillips 2020 from 0.01 to 0.07 M_sun
    Baraffe+15 from 0.075 - 0.4 M_sun
    Pisa 2011 from 0.5 - top of grid (~7 M_sun)

    MS (M > Pisa 2011)
    Ekstrom+12 for logAge < 7.4
    Parsec V1.2s for logAge > 7.4 

    metallicity = 'solar' --> Ekstrom+12 z014, Pisa2011 z015, Parsec z015

    if plot = True, will make plots of merged isochrones in 'plots' directory,
    which must already exist
    
    Code is expected to be run in merged model working directory.
    """
    # Root data directory for Ekstrom+12 isochrones
    rootDirE = models_dir + 'Ekstrom2012/iso/'
    metalPart = 'z014/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return
    rotPart = 'rot/'
    if not rotation:
        rotPart = 'norot/'
    rootDirE += metalPart+rotPart

    # Root data directory for the Baraffe isochrones
    rootDirBaraffe = models_dir + 'Baraffe15/iso/'

    # Root data directory for Pisa isochrones
    rootDirPisa = models_dir + 'Pisa2011/iso/'
    metSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return
    rootDirPisa += metSuffix

    # Root data directory for Parsec isochrones
    rootDirParsec = models_dir + 'ParsecV1.2s/iso/'
    metalSuffix = 'z015/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return        
    rootDirParsec += metalSuffix

    # Root data directory for Phillips isochrones
    rootDirPhillips = models_dir + 'Phillips2020/iso/'
    metalSuffix = 'z00/'
    if metallicity != 'solar':
        print( 'Non-solar metallicities not supported yet')
        return        
    rootDirPhillips += metalSuffix

    # Search both directories for iso_*.dat files
    isoFilesE = glob.glob(rootDirE + 'iso_*.dat')
    isoFilesB = glob.glob(rootDirBaraffe + 'iso_*.fits')
    isoFilesPi = glob.glob(rootDirPisa + 'iso_*.dat')
    isoFilesPa = glob.glob(rootDirParsec + 'iso_*')
    isoFilesPh = glob.glob(rootDirPhillips + 'iso_*.fits')

    # Output of merged isochrones
    if rotation == True:
        outDir = models_dir + 'merged/phillips_baraffe_pisa_ekstrom_parsec/{0}_rot/'.format(metSuffix[:-1])
        #outDir = models_dir + 'merged/test/{0}_rot/'.format(metSuffix[:-1])
    else:
        outDir = models_dir + 'merged/phillips_baraffe_pisa_ekstrom_parsec/{0}_norot/'.format(metSuffix[:-1])
        #outDir = models_dir + 'merged/test/{0}_norot/'.format(metSuffix[:-1])
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    # Isolate the iso*.dat file names
    for ii in range(len(isoFilesE)):
        isoFilesE[ii] = isoFilesE[ii].split('/')[-1]

    for ii in range(len(isoFilesB)):
        isoFilesB[ii] = isoFilesB[ii].split('/')[-1]    

    for ii in range(len(isoFilesPi)):
        isoFilesPi[ii] = isoFilesPi[ii].split('/')[-1]

    for ii in range(len(isoFilesPa)):
        isoFilesPa[ii] = isoFilesPa[ii].split('/')[-1]

    for ii in range(len(isoFilesPh)):
        isoFilesPh[ii] = isoFilesPh[ii].split('/')[-1]

    # Loop through the Pisa isochrones, adding the MS and Baraffe/Phillips models
    # as appropriate
    for ii in range(len(logAge_arr)):
        logAge = logAge_arr[ii]
        iso_name = 'iso_{0:.2f}.dat'.format(logAge)

        #-----PRE-MAIN SEQUENCE (only for logAge < 8)----#
        if logAge <= 8.0:
            # Merge with the Phillips 2020 models from 0.01 - 0.07 M_sun and
            # Baraffe+15 models from 0.075 - 0.4 M_sun. Includes
            # transition region between 0.4 - 0.5 M_sun in which we shift
            # from 100% Baraffe to 100% Pisa and 0.07 - 0.075 M_sun in which we 
            #shift from  100% Phillips to 100% Baraffe
        
            print( 'Merging isochrones Pisa + Phillips + Baraffe from ', iso_name)
            isoPMS = merge_isochrone_phillips_baraffe_pisa(logAge, metallicity=metallicity)

            #--------MAIN SEQUENCE-------#
            # Case where logAge <= 7.4, we merge with Ekstrom. Otherwise, merge
            # which parsec
            if logAge <= 7.4:
                if iso_name not in isoFilesE:
                    print( 'Skipping isochrones from ', iso_name)
                    print( 'PROBLEM WITH PISA OR EKSTROM')
                    pdb.set_trace()

                print( 'Merging isochrones Phillips+Pisa+Ekstrom from ', iso_name)
                iso = merge_isochrone_phillips_pisa_Ekstrom(logAge, metallicity=metallicity,
                                               rotation=rotation, iso_in=isoPMS)
            else:
                if iso_name not in isoFilesPa:
                    print( 'Skipping isochrones from ', iso_name)
                    print( 'PROBLEM WITH PISA OR PARSEC')
                    pdb.set_trace()
                
                print( 'Merging isochrones Phillips+Pisa+Parsec from ', iso_name)
                iso = merge_isochrone_phillips_pisa_parsec(logAge, metallicity=metallicity,
                                              iso_in=isoPMS)
        else:
            # If logAge > 8.0, just take the Phillips/Parsec model for the entire isochrone
            print( 'Making Phillips/Parsec from ', iso_name)
            # change to fits format
            iso = make_phillips_parsec_iso(logAge, metallicity=metallicity)
            
        # Make test plot, if desired. These are put in plots directory
        if plot:
            # Make different models different colors
            phillips_ind = np.where(iso.source == 'Phillips')
            merge_pb = np.where(iso.source == 'Phillips+Baraffe')
            baraffe_ind = np.where(iso.source == 'Baraffe')
            merge_ind = np.where(iso.source == 'Baraffe+Pisa')
            pisa_ind = np.where(iso.source == 'Pisa')
            MS_ind = np.where( (iso.source == 'Ekstrom') | (iso.source == 'Parsec'))

            # Temporarily turn off interactive plot. Turn back on at end
            py.ioff()
            
            py.figure(1)
            py.clf()
            py.plot(iso.logT[phillips_ind], iso.logL[phillips_ind], 'c-', label = 'Phillips',
                    linewidth=2)
            py.plot(iso.logT[merge_pb], iso.logL[merge_pb], 'y-', label = 'Phillips+Baraffe',
                    linewidth=2)
            py.plot(iso.logT[baraffe_ind], iso.logL[baraffe_ind], 'k-', label = 'Baraffe',
                    linewidth=2)
            py.plot(iso.logT[merge_ind], iso.logL[merge_ind], 'r-', label = 'Baraffe+Pisa',
                    linewidth=2)
            py.plot(iso.logT[pisa_ind], iso.logL[pisa_ind], 'g-', label = 'Pisa',
                    linewidth=2)
            py.plot(iso.logT[MS_ind], iso.logL[MS_ind], 'b-',
                    label = 'Ekstrom/Parsec', linewidth=2)
            py.xlabel('log Teff')
            py.ylabel('log L')
            py.title('Log Age = {0}'.format(logAge))
            py.legend(loc=3)
            py.axis([4.9, 3.2, -3.0, 8])
            py.savefig(outDir+'plots/iso_{0:.2f}.png'.format(logAge))

            py.ion()
            
            
        _out = open(outDir + iso_name, 'w')

        hdr_fmt = '%12s  %10s  %10s  %10s %10s %12s %5s %-10s\n'
        _out.write(hdr_fmt % 
                   ('# M_init', 'log T', 'log L', 'log g', 'logT_WR', 'M_curr', 'phase', 'Source'))
        _out.write(hdr_fmt % 
                   ('# (Msun)', '(Kelvin)', '(Lsun)', '(cgs)', '(Kelvin)', '(Msun)', '()', '()'))

        for kk in range(len(iso.mass)):
            _out.write('%12.6f  %10.4f  %10.4f  %10.4f  %10.4f %12.6f %5d %-10s\n' %
                       (iso.mass[kk], iso.logT[kk], iso.logL[kk],
                        iso.logg[kk], iso.logT_WR[kk],
                        iso.mass_current[kk], iso.phase[kk], iso.source[kk]))

        _out.close()
    return


##### CREATING MERGED ATMOSPHERE MODEL #####
def create_merged_models(cdbs_path, plot=False):
    """
    for 1200 - 1000 K, merge Meisner and BTSettl!
    
    From 1000 K - 1200 K, merge the Meisner and BTSettl atmospheres.
    More like BTSettl near 1200 K, more like Meisner near 1000 K

    cdbs_path is path to cdbs directory (including cdbs). Will make new directory
    in cdbs/grid named "merged_meisner_BTsettl" with with merged models. If plot = True, will plot
    the normalized merged model plus the original Meisner and BTSettl models

    Temp 1000 - 1200, steps of 20; logg 2.5 - 5.5, steps of 0.5, metallicity covering
    Meisner range (2.5 -- 5.5, in steps of 0.5). (So, this is one model at 5250)?

    Note: metallicity directories created by hand
    
    Creates new directory "merged_meisner_BTSettl" in cdbs/grid with new spectrum + catalog file.
    Also includes the atlas 5500K model and phoenix 5000K model, for interpolation purposes
    """
    #Setting logg sampling
    logg_arr = np.arange(2.5, 5.5+0.1, 0.5)

    # Setting metallicity sub-directories
    atlas_dir = ['ckm25', 'ckm20', 'ckm15', 'ckm10', 'ckm05', 'ckp00', 'ckp02', 'ckp05']
    phoenix_dir = ['phoenixm30', 'phoenixm20', 'phoenixm15', 'phoenixm10', 'phoenixm05', 'phoenixm00', 'phoenixp05', 'phoenixp05']
    output_dir = ['mergedm25', 'mergedm20', 'mergedm15', 'mergedm10', 'mergedm05', 'mergedp00', 'mergedp02', 'mergedp05']
    metal_arr = [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.2, 0.5]
    
    assert len(atlas_dir) == len(phoenix_dir) == len(output_dir)

    # Make new cdbs merged directory, if it doesn't already exist
    newPath = '{0}/grid/merged_atlas_phoenix'.format(cdbs_path)
    if not os.path.exists(newPath):
        os.mkdir(newPath)

    # For each metallicity, create merged model
    for ii in range(len(output_dir)):
        # Make metallicity dir for merged model
        final_dir = '{0}/{1}'.format(newPath, output_dir[ii])
        if not os.path.exists(final_dir):
            os.mkdir(final_dir)
        
        # Extract the relevant ATLAS and PHEONIX models near 5250 K
        atlas_path = '{0}/grid/ck04models/{1}/'.format(cdbs_path, atlas_dir[ii])
        atlas_hdu = fits.open('{0}/{1}_5250.fits'.format(atlas_path, atlas_dir[ii]))
        atlas_5250 = atlas_hdu[1].data
        phoenix_path = '{0}/grid/phoenix_v16_rebin/{1}'.format(cdbs_path, phoenix_dir[ii])
        phoenix_hdu1 = fits.open('{0}/{1}_05200.fits'.format(phoenix_path, phoenix_dir[ii]))
        phoenix_hdu2 = fits.open('{0}/{1}_05300.fits'.format(phoenix_path, phoenix_dir[ii]))
        phoenix_5200 = phoenix_hdu1[1].data
        phoenix_5300 = phoenix_hdu2[1].data
        print('Done reading input models')

        # Trim both models to the wavelength region we want: VRIJHKL (0.25 - 4.2 mircons)
        good = np.where( (atlas_5250['Wavelength'] > 2500) &
                        (atlas_5250['Wavelength'] < 52000) )
        atlas_5250 = atlas_5250[good]
        good = np.where( (phoenix_5200['Wavelength'] > 2500) &
                        (phoenix_5200['Wavelength'] < 52000) )
        phoenix_5200 = phoenix_5200[good]
        phoenix_5300 = phoenix_5300[good]

        #----------------------------#
        # For each logg val, create phoenix spectrum for 5250 K,
        # degrade resolution to match atlas, then average with
        # 5250 K atlas model to get merged model at 5250 K
        #----------------------------#
        new_model = []
        phoenix_5250_f = []
        atlas_5250_f = []
        for logg in logg_arr:
            # Create phoenix models at 5250 K
            low = phoenix_5200['g{0:2.1f}'.format(logg)]
            high = phoenix_5300['g{0:2.1f}'.format(logg)]

            arr = np.transpose([low, high])
            phoenix_5250 = np.mean(arr, axis=1)

            phoenix_5250_rebin = rebin_spec(phoenix_5200['Wavelength'], phoenix_5250,
                                            atlas_5250['Wavelength'])

            # Store phoenix rebinned average spectrum
            phoenix_5250_f.append(phoenix_5250_rebin)

            # Store atlas spectrum for later
            grav = str(logg).split('.')        
            atlas_5250_f.append(atlas_5250['g'+grav[0]+grav[1]])
        
            # Now, final 5250 K model will be average of atlas and phoenix
            # models (since exactly inbetween 5000 K and 5500 K)
            arr = np.transpose([phoenix_5250_rebin, atlas_5250['g'+grav[0]+grav[1]]])
            final = np.mean(arr, axis=1)

            new_model.append(final)
            print('Done with logg {0:2.1f}'.format(logg))

        # Create fits table with new model. First column is wavelength, followed by
        # fluxes for different logg
        wave = atlas_5250['Wavelength'] # Still same wavelength array as before
        c0 = fits.Column(name='Wavelength', format='D', array=wave)
        c1 = fits.Column(name='g0.0', format='E', array=new_model[0])
        c2 = fits.Column(name='g0.5', format='E', array=new_model[1])
        c3 = fits.Column(name='g1.0', format='E', array=new_model[2])
        c4 = fits.Column(name='g1.5', format='E', array=new_model[3])
        c5 = fits.Column(name='g2.0', format='E', array=new_model[4])
        c6 = fits.Column(name='g2.5', format='E', array=new_model[5])
        c7 = fits.Column(name='g3.0', format='E', array=new_model[6])
        c8 = fits.Column(name='g3.5', format='E', array=new_model[7])
        c9 = fits.Column(name='g4.0', format='E', array=new_model[8])
        c10 = fits.Column(name='g4.5', format='E', array=new_model[9])
        c11 = fits.Column(name='g5.0', format='E', array=new_model[10])

        cols = fits.ColDefs([c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11])
        tbhdr = phoenix_hdu1[1].header
        prihdr = phoenix_hdu1[0].header 
        prihdu = fits.PrimaryHDU(header=prihdr) 
        tbhdu = fits.BinTableHDU.from_columns(cols, header=tbhdr)
        # Add TUNIT keysto tbhdu
        tbhdu = make_merged_header(tbhdu)

        # Make final hdu table, save it
        finalhdu = fits.HDUList([prihdu, tbhdu])
        finalhdu.writeto('{0}/{1}_05250.fits'.format(final_dir, output_dir[ii]), clobber=True)
        
        # Test plot, if desired
        if plot == True:
            py.figure(1, figsize=(10,10))
            py.clf()
            # Plot merged model
            py.semilogy(wave, wave * new_model[8], 'r-', label='Merged')
            # Plot atlas 5250 K model
            py.semilogy(wave, wave * atlas_5250_f[8], 'b-', label = 'Atlas')
            # Plot phoenix 5250 K model
            py.semilogy(wave, wave * phoenix_5250_f[8], 'g-', label = 'Phoenix')
            py.legend()
            py.xlabel('Wavelength (Angstrom)')
            py.ylabel(r'log ($\lambda$ F$_{\lambda}$)')
            py.axis([5000, 40000, 2*10**9, 4*10**10])
            py.title('Merged 5250 K Spectrum')
            py.savefig('Merged_atlas_phoenix_{0}.png'.format(output_dir[ii]))
            
            pdb.set_trace()
            py.close('all')

        # Copy the atlas 5500 K and phoenix 5000 K models into the
        # merged directory to accompany the merged model. This is for
        # interpolation purposes for synphot
        cmd = 'cp {0}/{1}_5500.fits {2}'.format(atlas_path, atlas_dir[ii], final_dir)
        cmd2 = 'cp {0}/{1}_05000.fits {2}'.format(phoenix_path, phoenix_dir[ii], final_dir)
        os.system(cmd)
        os.system(cmd2)

        cmd = 'mv {0}/{1}_5500.fits {0}/{2}_05500.fits'.format(final_dir, atlas_dir[ii], output_dir[ii])
        cmd2 = 'mv {0}/{1}_05000.fits {0}/{2}_05000.fits'.format(final_dir, phoenix_dir[ii], output_dir[ii])
        os.system(cmd)
        os.system(cmd2)
        
        # Close all open hdu
        atlas_hdu.close()
        phoenix_hdu1.close()
        phoenix_hdu2.close()

        print('Done {0}'.format(output_dir[ii]))
    
    # Make catalog.fits table for new merged model, plus the atlas and
    # phoenix models at the high and low extremes of the temp range. Note temp
    # of merged model goes first, for filename purposes
    catalog, filenames = make_merged_catalog(output_dir, [5250, 5000, 5500], logg_arr, metal_arr)
    catalog.write(cdbs_path+'/grid/merged_atlas_phoenix/catalog.fits', overwrite=True)

    return



def create_merged_meisner_btsettl(cdbs_path, plot=False):
    """
    From 1000 K - 1200 K, merge the BTSettl_CIFITS2011_2015 and Meisner 2023 atmospheres.
    BTSettl at 1200 K, Meisner near 1000 K

    cdbs_path is path to cdbs directory (including cdbs). Will make new directory
    in cdbs/grid named "merged_meisner_btsettl" with with merged models. If plot = True, will plot
    the normalized merged model plus the original models

    
    (Only 1 truely merged model at 3500 K, with log g from 2.5 - 5.5)?
    
    Creates new directory "merged_meisner_btsettl" in cdbs/grid with new spectrum + catalog file.
    (Also includes the atlas 5500K model and phoenix 5000K model, for interpolation purposes)?
    """
    # Set log g sampling
    logg_solar_arr = np.arange(3.5, 5.5+0.1, 0.5)
    logg_nonsolar_arr = np.array([3.5, 5.0])

    # Setting metallicity sub-directories
    meisner_dir = ['mm05', 'mp00']
    output_dir = ['mergedm05', 'mergedp00']
    metal_arr = [-0.5, 0]
    
    assert len(meisner_dir) == len(output_dir)

    # Make new cdbs merged directory, if it doesn't already exist
    newPath = '{0}/grid/merged_meisner_BTSettl'.format(cdbs_path)
    if not os.path.exists(newPath):
        os.mkdir(newPath)

   # For each metallicity, create merged model
    for ii in range(len(output_dir)):
        # Make metallicity dir for merged model
        final_dir = '{0}/{1}'.format(newPath, output_dir[ii])
        if not os.path.exists(final_dir):
            os.mkdir(final_dir)

        # For each gravity, extract BTSettl and Meisner models at 1100 K and
        # average them to get the merged model. Also write BTSettl models at
        # 1200 K and Meisner models at 1000 K
        if metal_arr[ii] == 0:
            logg_tmp = logg_solar_arr
            BT_func = atmospheres.get_BTSettl_atmosphere
        else:
            BT_func = atmospheres.get_BTSettl_atmosphere
            logg_tmp = logg_nonsolar_arr

        for jj in logg_tmp:
            BT_atmo = BT_func(temperature=1100, gravity=jj, rebin=True)
            meisner_atmo = atmospheres.get_Meisner2023_atmosphere(temperature=1100, gravity=jj, rebin=True)

            # Check to make sure wavelengths are the same between atmospheres
            print("BT_func assigned to:", BT_func)
            print(f"BT_func({1100}, {jj}, rebin=True) -> {BT_atmo}")
            diff = BT_atmo.wave - meisner_atmo.wave
            if np.sum(diff > 0):
                print('Wavelength mismatch problem!')
                pdb.set_trace()

            merged_flux = np.mean(np.array([meisner_atmo.flux, BT_atmo.flux]), axis=0)
            
            # Create fits file with new atmosphere
            c0 = fits.Column(name='Wavelength', format='D', array=BT_atmo.wave)
            c1 = fits.Column(name='Flux', format='E', array=merged_flux)

            cols = fits.ColDefs([c0, c1])
            tbhdu = fits.BinTableHDU.from_columns(cols)

            prihdu = fits.PrimaryHDU()
            tbhdu.header['TUNIT1'] = 'ANGSTROM'
            tbhdu.header['TUNIT2'] = 'FLAM'
            hdu_new = fits.HDUList([prihdu, tbhdu])

            # Save fits file to merged Meisner-BTSettl direcotry
            hdu_new.writeto('{0}/{1}_03500_{2}.fits'.format(final_dir, output_dir[ii], jj), overwrite=True)

            # Make test plot, if desired
            if plot:
                py.figure(1, figsize=(10,10))
                py.clf()
                py.semilogy(BT_atmo.wave, BT_atmo.wave * merged_flux, 'r-', label='Merged')
                py.semilogy(BT_atmo.wave, BT_atmo.wave * BT_atmo.flux, 'b-', label = 'BTSettl')
                py.semilogy(meisner_atmo.wave, meisner_atmo.wave * meisner_atmo.flux, 'g-', label = 'Meisner')
                py.legend()
                py.xlabel('Wavelength (Angstrom)')
                py.ylabel(r'log ($\lambda$ F$_{\lambda}$)')
                py.xlim(5000, 40000)
                py.ylim(10**8, 10**10)
                py.title('Merged 1100 K Spectrum, Z = {1}, logg = {0}'.format(jj, metal_arr[ii]))
                py.savefig('Merged_Meisner_BTSettl_{1}_{0}.png'.format(jj, metal_arr[ii]))

            # Now to bring over the 1000 K Meisner and 1200 K BTSettl models
            BT_atmo = BT_func(temperature=1200, gravity=jj, rebin=True)
            meisner_atmo = atmospheres.get_Meisner2023_atmosphere(temperature=1000, gravity=jj, rebin=True)

            # Create new fits files for these as well
            c0_b = fits.Column(name='Wavelength', format='D', array=BT_atmo.wave)
            c1_b = fits.Column(name='Flux', format='E', array=BT_atmo.flux)
            c0_m = fits.Column(name='Wavelength', format='D', array=meisner_atmo.wave)
            c1_m = fits.Column(name='Flux', format='E', array=meisner_atmo.flux)
        
            cols_b = fits.ColDefs([c0_b, c1_b])
            tbhdu_b = fits.BinTableHDU.from_columns(cols_b)
            cols_m = fits.ColDefs([c0_m, c1_m])
            tbhdu_m = fits.BinTableHDU.from_columns(cols_m)
            
            prihdu = fits.PrimaryHDU()
            tbhdu_b.header['TUNIT1'] = 'ANGSTROM'
            tbhdu_b.header['TUNIT2'] = 'FLAM'
            tbhdu_m.header['TUNIT1'] = 'ANGSTROM'
            tbhdu_m.header['TUNIT2'] = 'FLAM'
            
            hdu_newb = fits.HDUList([prihdu, tbhdu_b])
            hdu_newm = fits.HDUList([prihdu, tbhdu_m])
            
            # Save fits file to merged Meisner-BTSettl direcotry
            hdu_newb.writeto('{0}/{1}_03200_{2}.fits'.format(final_dir, output_dir[ii], jj), clobber=True)
            hdu_newp.writeto('{0}/{1}_03800_{2}.fits'.format(final_dir, output_dir[ii], jj), clobber=True)
            hdu_new.close()
            hdu_newb.close()
            hdu_newp.close()
            
        print('Done {0}'.format(output_dir[ii]))
            
    return

def make_catalog_merged_models2(path='/g/lu/models/cdbs/grid/merged_BTSettl_phoenix/'):
    """
    Make cdbs catalog.fits file for merged BTSettl/phoenix direcotry. path should
    point to this directory.

    Writes catalog.fits file in the cdbs directory
    """
    output_dir = ['mergedm25', 'mergedm20', 'mergedm15', 'mergedm10', 'mergedm05', 'mergedp00', 'mergedp02', 'mergedp05']
    metal_arr = [-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.2, 0.5]

    index_str = []
    name_str = []
    for ii in range(len(output_dir)):
        files = glob.glob('{0}/{1}/*.fits'.format(path, output_dir[ii]))
    
        # Extract parameters for each atmosphere from the filename,
        # construct columns for catalog file
        for name in files:
            final = name.split('/')[-1]
            tmp = final.split('_')
            temp = float(tmp[1]) # In kelvin
            logg = float(tmp[2][:-5])

            index_str.append('{0},{1},{2:3.2f}'.format(int(temp), metal_arr[ii], logg))
            name_str.append('{0}/{1}[Flux]'.format(output_dir[ii], final))

    # Make catalog
    catalog = Table([index_str, name_str], names = ('INDEX', 'FILENAME'))

    # Create catalog.fits file in directory with the models
    catalog.write(path+'catalog.fits', format = 'fits', overwrite=True)
    
    return