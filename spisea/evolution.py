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
import astropy.units as u
import astropy.constants as c
from astropy import coordinates as coords

logger = logging.getLogger('evolution')

# Fetch root directory of evolution models.
try:
    models_dir = os.environ['SPISEA_MODELS']
    models_dir += '/evolution/'
except KeyError:
    warnings.warn("SPISEA_MODELS is undefined; functionality "
                  "will be SEVERELY crippled.")
    models_dir = ''

# Function to get installed evo grid number
def get_installed_grid_num(input_models_dir):
    """
    Get installed grid number
    """
    # Define the installed model grid number
    file_name = input_models_dir + '/grid_version.txt'

    # Read in the file. In the case where it doesn't
    # exist, then grid version is assumed to be 1.0
    # (since this didn't always exist)
    try:
        file1 = open(file_name, 'r')
        read = file1.readlines()
        evo_grid_num = float(read[1])
        file1.close()
    except FileNotFoundError:
        evo_grid_num = 1.0

    return evo_grid_num

# Function to check evo grid version number
def check_evo_grid_number(required_num, input_models_dir):
    """
    Check if installed grid meets the required
    grid version number. Installed grid number must
    be greater than or equal to this number
    """

    # Get installed gridnumber
    grid_num = get_installed_grid_num(input_models_dir)

    # Check: is installed grid number < required_num?
    # If not, raise mismatch error
    if grid_num < required_num:
        raise exceptions.ModelMismatch(required_num, grid_num, 'evolution')

    return grid_num

class StellarEvolution(object):
    """
    Base Stellar evolution class.

    Parameters
    ----------
    model_dir: path
        Directory path to evolution model files

    age_list: list
        List of ages

    mass_list: list
        List of masses

    z_list: list
        List of metallicities

    model_version_name: string
        Name of the model class plus additional details like version
        numbers and rotation if relevant.
    """
    def __init__(self, model_dir, age_list, mass_list, z_list):
        self.model_dir = model_dir
        self.z_list = z_list
        self.mass_list = mass_list
        self.age_list = age_list
        self.external_evol = False
        self.model_version_name = "None"

        return

class Geneva(StellarEvolution):
    def __init__(self):
        r"""
        Define intrinsic properties for Geneva stellar models.
        """
        self.model_version_name = "Geneva"
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]

        # define metallicity parameters for Geneva models
        z_list = [0.01, 0.02, 0.03]

        # populate list of isochrone ages (log scale)
        age_list = [round(5.5 + 0.01*i, 2) for i in range(190)]
        age_list += [round(7.4 + 0.05*i, 2) for i in range(12)]
        age_list += [round(math.log10(1.e8*x), 2) for x in range(1, 10)]
        age_list += [round(math.log10(1.e9*x), 2) for x in range(1, 10)]
        age_list = age_list

        # specify location of model files
        self.model_dir = models_dir + 'geneva/'
        
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)

        self.z_solar = 0.02
        self.z_file_map = {0.01: 'z01/', 0.02: 'z02/', 0.03: 'z03/'}

        # Define required evo_grid number
        self.evo_grid_min = 1.0

    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Geneva collection.
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)

        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # convert age (in yrs) to log scale and find nearest value in grid
        log_age = np.log10(age)

        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_' + str(self.age_list[age_idx]) + '.fits'

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file

        # return isochrone data
        return genfromtxt(full_iso_file, comments='#')


#---------------------------------------#
# Now for the Ekstrom+12 Geneva models
#---------------------------------------#

class Ekstrom12(StellarEvolution):
    """
    Evolution models from
    `Ekstrom et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.146E/abstract>`_.

    Downloaded from `website <http://obswww.unige.ch/Recherche/evoldb/index/Isochrone/>`_.

    Parameters
    ----------
    rot: boolean, optional
        If true, then use rotating Ekstrom models. Default is true.
    """
    def __init__(self, rot=True):
        if rot:
            self.model_version_name = "Ekstrom12-rot"
        else:
            self.model_version_name = "Ekstrom12-norot"
        # define metallicity parameters for Ekstrom+12 models
        self.z_list = [0.014]

        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.0, 8.0+0.005, 0.01)

        # Specify location of model files
        self.model_dir = models_dir+'Ekstrom2012/'

        # Specifying metallicity
        self.z_solar = 0.014
        self.z_file_map = {0.014: 'z014/'}

        # Specify rotation or not
        self.rot = rot

        # Define required evo_grid number
        self.evo_grid_min = 1.0

    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Ekstrom+12 Geneva collection.
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)

        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            raise ValueError(f'Requested metallicity z_solar * 10^{metallicity} = {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        if self.rot:
            full_iso_file = self.model_dir + 'iso/' + z_dir + 'rot/' + iso_file
        else:
            full_iso_file = self.model_dir + 'iso/' + z_dir + 'norot/' + iso_file

        # Return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col4', 'Z')
        iso.rename_column('col1', 'logAge')
        iso.rename_column('col3', 'mass')
        iso.rename_column('col6', 'mass_current')
        iso.rename_column('col7', 'logL')
        iso.rename_column('col8', 'logT')
        iso.rename_column('col22', 'logg')
        iso.rename_column('col9', 'logT_WR')

        # Add isWR column
        isWR = Column([False] * len(iso), name='isWR')
        idx_WR = np.where(iso['logT'] != iso['logT_WR'])
        isWR[idx_WR] = True
        iso.add_column(isWR)

        # Add a phase column... everything is just a star.
        iso.add_column( Column(np.ones(len(iso)), name = 'phase'))

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

    def format_isochrones(input_iso_dir):
        r"""
        Parse iso.fits (filename hardcoded) file downloaded from Ekstrom+12
        models, create individual isochrone files for the different ages.

        input_iso_directory should lead to
            Ekstrom2012/iso/<metallicity>
        directory, where iso.fits file should be located.

        Creates two new directories, rot and norot, which contain their
        respective isochrones.
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into metallicity direcotry, read iso.fits file
        os.chdir(input_iso_dir)

        print( 'Read Input: this is slow')
        iso = Table.read('iso.fits')
        print( 'Done'    )

        ages_all = iso['col1']

        # Extract the unique ages
        age_arr = np.unique(ages_all)

        # For each unique age, extract the proper rows and make corresponding
        # table. Be sure to separate rotating from non-rotating, and put in
        # separate subdirectories

        # First make the rot and norot directories, if they don't exit
        if os.path.exists('rot'):
            pass
        else:
            os.mkdir('rot')
            os.mkdir('norot')

        print( 'Making individual isochrone files')
        for age in age_arr:
            good = np.where(ages_all == age)

            # Identify rot vs. non-rot
            idx_r = np.where(iso[good]['col2'] == 'r')
            idx_n = np.where(iso[good]['col2'] == 'n')

            tmp_r = iso[good][idx_r]
            tmp_n = iso[good][idx_n]

            # Write tables
            tmp_r.write('rot/iso_{0:4.2f}.fits'.format(age))
            tmp_n.write('norot/iso_{0:4.2f}.fits'.format(age))

        # Return to starting directory
        os.chdir(start_dir)

        return

    def create_iso(fileList, ageList, rot=True):
        """
        Given a set of isochrone files downloaded from
        the server, put in correct
        iso.fits format for parse_iso code.

        fileList: list of downloaded isochrone files (could be one)

        ageList: list of lists of ages associated with each file in filelist.
        MUST BE IN SAME ORDER AS ISOCHRONES IN FILE! Also needs to be in logAge

        rot = TRUE: assumes that models are rotating, will add appropriate column

        This code writes the individual files, which is then easiest to combine by hand
        in aquamacs
        """
        # Read each file in fileList individually, add necessary columns
        for i in range(len(fileList)):
            t = Table.read(fileList[i],format='ascii')
            ages = ageList[i]

            # Find places where new models start; mass here is assumed to be 0.8
            start = np.where(t['M_ini'] == 0.8)

            # Now, each identified start is assumed to be associated with the
            # corresponding age in ages
            if len(start[0]) != len(ages):
                print( 'Ages mismatched in file! Quitting...')
                return

            age_arr = np.zeros(len(t))


            for j in range(len(start[0])):
                low_ind = start[0][j]
                # Deal with case at end of file
                if (j == len(start[0])-1):
                    high_ind = len(t)
                else:
                    high_ind = start[0][j+1]

                ind = np.arange(low_ind, high_ind, 1)
                age_arr[ind] = ages[j]

            # Add ages_arr column to column 1 in ischrone, as well as column
            # signifying rotation
            col_age = Column(age_arr, name = 'logAge')
            rot_val = np.chararray(len(t))
            rot_val[:] = 'r'
            if not rot:
                rot_val[:] = 'n'

            col_rot = Column(rot_val, name='Rot')

            t.add_column(col_rot, index=0)
            t.add_column(col_age, index=0)

            t.write('tmp'+str(i)+'.fits')

        return

#---------------------------------------#
# Now for the Parsec version 1.2s models
#---------------------------------------#

class Parsec(StellarEvolution):
    """
    Evolution models from
    `Bressan et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract>`_,
    version 1.2s.

    Downloaded from `here <http://stev.oapd.inaf.it/cgi-bin/cmd>_`

    Notes
    -----
    Evolution model parameters used in download:

    * n_Reimers parameter (mass loss on RGB) = 0.2
    * photometric system: HST/WFC3 IR channel
    * bolometric corrections OBC from Girardi+08, based on ATLAS9 ODFNEW models
    * Carbon star bolometric corrections from Aringer+09
    * no dust
    * no extinction
    * Chabrier+01 mass function
    """
    def __init__(self):
        r"""
        Define intrinsic properties for the Parsec version 1.2s stellar
        models.
        """
        # populate list of model masses (in solar masses)
        self.model_version_name = "Parsec1.2s"
        #mass_list = [(0.1 + i*0.005) for i in range(181)]

        # define metallicity parameters for Parsec models
        self.z_list = [0.005, 0.015, 0.04]

        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.6, 10.12+0.005, 0.01)
        self.age_list = np.append(6.40, self.age_list)

        # Specify location of model files
        self.model_dir = models_dir+'ParsecV1.2s/'

        # Specifying metallicity
        self.z_solar = 0.015
        self.z_file_map = {0.005: 'z005/', 0.015: 'z015/', 0.04: 'z04/'}

        # Define required evo_grid number
        self.evo_grid_min = 1.0

    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Parsec version 1.2s
        collection.
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)

        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file

        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col1', 'Z')
        iso.rename_column('col2', 'logAge')
        iso.rename_column('col3', 'mass')
        iso.rename_column('col4', 'mass_current')
        iso.rename_column('col5', 'logL')
        iso.rename_column('col6', 'logT')
        iso.rename_column('col7', 'logg')
        iso.rename_column('col15', 'phase')
        iso['logT_WR'] = iso['logT']

        # Parsec doesn't identify WR stars, so identify all as "False"
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso


    def format_isochrones(input_iso_dir, metallicity_list):
        r"""
        Parse isochrone file downloaded from Parsec version 1.2 for different
        metallicities, create individual isochrone files for the different ages.

        input_iso_dir: points to ParsecV1.2s/iso directory. Assumes metallicity
        subdirectories already exist with isochrone files downloaded in them
        (isochrones files expected to start with "output*")

        metallicity_list format: absolute (vs. relative to solar),
        z + <digits after decimal>: e.g. Z = 0.014 --> z014
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)

        # Work on each metallicity isochrones individually
        for metal in metallicity_list:
            # More into metallicity directory, read isochrone file
            os.chdir(metal)

            isoFile = glob.glob('output*')
            print( 'Read Input: this is slow')
            iso = Table.read(isoFile[0], format='fits')
            print( 'Done')

            ages_all = iso['col2']

            # Extract the unique ages
            age_arr = np.unique(ages_all)

            # For each unique age, extract the proper rows and make corresponding
            # table
            print( 'Making individual isochrone files')
            for age in age_arr:
                good = np.where(ages_all == age)
                tmp = iso[good]

                #Write table
                tmp.write('iso_{0:4.2f}.fits'.format(age))

            # Move back into iso directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return

class Phillips2020(StellarEvolution):
    """
    Evolution models from 
    `Phillips et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020A%26A...637A..38P/abstract>`_.

    Downloaded from `here <https://noctis.erc-atmo.eu/fsdownload/zyU96xA6o/phillips2020>_`

    Notes
    -----
    Evolution model parameters used in download:

    * Assume chemical equilibrium
    * Solar metallicity
    * For young BDs
    CHANGE!!!
    """

    def __init__(self):
        r"""
        Define intrinsic properties for the Phillips brown dwarf stellar models
        """
        # specify location of model files
        self.model_dir = models_dir + 'Phillips2020/'

        # specifying metallicity
        self.z_list = [0.015]
        self.z_file_map = {0.015: 'z00'}
        self.z_solar = 0.015

        # populate list of isochrone ages (log scale)
        self.age_list = np.array([6.0, 6.102564101374986, 6.205128204112327, 6.307692307498837, 
                         6.41025641041525, 6.512820513399242, 6.615384615079123, 6.717948717593952,
                         6.820512820432055, 6.923076923041517,7.025641027630915, 7.128205127364955,
                         7.230769230806862, 7.333333333326906, 7.43589743645667, 7.538461537884313,
                         7.64102564151456, 7.743589743588999, 7.846153846357851, 7.948717948583871,
                         8.051282050278353, 8.153846153753816, 8.256410257173776, 8.35897435818861,
                         8.461538460830697, 8.564102564448598, 8.666666666328634, 8.769230769062688,
                         8.87179487160516, 8.974358974304664, 9.076923076299368, 9.179487179310316,
                         9.282051282920115, 9.384615385138028, 9.487179487213739, 9.589743589709178,
                         9.692307692157154, 9.794871794783068, 9.89743589751841, 10.0])

        # define required evo_grid number
        self.evo_grid_min = 3.0

    def isochrone(self, age= 1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Phillips2020 collection.
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)
        
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # Specify subdirectory for metallicity
        iso_path = os.path.join(self.model_dir, 'iso', z_dir)

        # Find nearest age in grid to input grid by parsing through available files
        close_age = self.age_list[np.argmin(abs(self.age_list - log_age))]
        close_file = iso_path+f'/iso_{close_age:.15f}.fits'
        print(f"Found nearest age {close_age} for requested age of {log_age}")
        
        # Make sure the closest file exists
        if not os.path.exists(close_file):
            raise FileNotFoundError(f"Isochrone file not found: {close_file}.")
        
        # return isochrone data
        iso = Table.read(close_file, format='fits')
        iso.rename_column('Z', 'Z')
        iso.rename_column('Age', 'logAge')
        iso.rename_column('Mass', 'mass')
        iso.rename_column('Mass_current', 'mass_current')
        iso.rename_column('Luminosity', 'logL')
        iso.rename_column('Teff', 'logT')
        iso.rename_column('Gravity', 'logg')
        iso['logT_WR'] = iso['logT']

        # Phillips doesn't identify WR stars, so identify all as "False"
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)
        
        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = metallicity

        return iso

    def format_isochrones(input_iso_dir, metallicity_list):
        r"""
        Change
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)
        
        # Work on each metallicity isochrones individually
        for metal in metallicity_list:
            # More into metallicity directory, read isochrone file
            os.chdir(metal)

            isoFile = glob.glob('output*')
            print( 'Read Input: this is slow')
            iso = Table.read(isoFile[0], format='fits')
            print( 'Done')
    
            ages_all = iso['col2']

            # Extract the unique ages
            age_arr = np.unique(ages_all)

            # For each unique age, extract the proper rows and make corresponding
            # table
            print( 'Making individual isochrone files')
            for age in age_arr:
                good = np.where(ages_all == age)
                tmp = iso[good]

                #Write table
                tmp.write('iso_{0:4.2f}.fits'.format(age))

            # Move back into iso directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return
        

class Marley2021(StellarEvolution):
    """
    Evolution models from 
    `Marley et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021ApJ...920...85M/abstract>`_.

    Downloaded from `here <https://zenodo.org/records/5063476>_`

    Notes
    -----
    Evolution model parameters used in download:

    * 
    CHANGE!!!
    """
    def __init__(self):
        r"""
        Define intrinsic properties for the Marley brown dwarf stellar models.
        """
        # populate list of model masses (in solar masses)
        #mass_list = [(0.1 + i*0.005) for i in range(181)]
        
        # define metallicity parameters for Parsec models
        self.z_solar = 0.0142
        self.z_list = [self.z_solar * (10.**m) for m in [-0.5, 0.0, 0.5]]
        
        # populate list of isochrone ages (log scale)
        self.age_list = [10.0, 7.0, 8.0, 9.0, 6.0, 7.176091259055681, 8.176091259055681, 9.176091259055681, 6.301029995663981, 
                         7.301029995663981, 8.301029995663981, 9.301029995663981, 6.477121254719663, 7.477121254719663, 
                         8.477121254719663, 9.477121254719663, 6.6020599913279625, 7.6020599913279625, 8.602059991327963, 
                         9.602059991327963, 6.778151250383644, 7.778151250383644, 8.778151250383644, 9.778151250383644, 
                         6.903089986991944, 7.903089986991944, 8.903089986991944, 9.903089986991944]
        
        # Specify location of model files
        self.model_dir = models_dir+'Marley2021/'

        # Specifying metallicity
        self.z_file_map = {
            self.z_list[0]: 'zm05/', 
            self.z_list[1]: 'zp00/', 
            self.z_list[2]: 'zp05/'
        }

        # Define required evo_grid number
        self.evo_grid_min = 3.0

    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Marley2021 collection.
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)
        
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))
        
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # Find closest age in grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = f'iso_{self.age_list[age_idx]}.fits'

        # Create path to iso file
        full_iso_file = os.path.join(self.model_dir, 'iso', z_dir, iso_file)

        print(f"Found nearest age file as {full_iso_file} for requested age of {log_age}")
        
        # Make sure the closest file exists
        #if not os.path.exists(close_file):
            #raise FileNotFoundError(f"Isochrone file not found: {close_file}.")
        
        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('Z', 'Z')
        iso.rename_column('Age', 'logAge')
        iso.rename_column('Mass', 'mass')
        iso.rename_column('Mass_current', 'mass_current')
        iso.rename_column('log_L', 'logL')
        iso.rename_column('Teff', 'logT')
        iso.rename_column('logg', 'logg')
        iso.rename_column('Radius', 'radius')
        iso['logT_WR'] = iso['logT']

        # Marley doesn't identify WR stars, so identify all as "False"
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)
        
        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(z_defined / self.z_solar)

        return iso

    def format_isochrones(input_iso_dir, metallicity_list):
        r"""
        Parse isochrone files downloaded from Marley 2021 for different
        metallicities, create individual isochrone files for the different ages.
    
        input_iso_dir: points to Marley2021/iso directory. Assumes metallicity
        subdirectories already exist with isochrone files downloaded in them
        (isochrones files expected to start with "output*")

        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)
        
        # Work on each metallicity isochrones individually
        for metal in metallicity_list:
            # More into metallicity directory, read isochrone file
            os.chdir(metal)

            isoFile = glob.glob('output*')
            print( 'Read Input: this is slow')
            iso = Table.read(isoFile[0], format='fits')
            print( 'Done')
    
            ages_all = iso['Age']

            # Extract the unique ages
            age_arr = np.unique(ages_all)

            # For each unique age, extract the proper rows and make corresponding
            # table
            print( 'Making individual isochrone files')
            for age in age_arr:
                good = np.where(ages_all == age)
                tmp = iso[good]

                #Write table
                tmp.write('iso_{0:4.2f}.fits'.format(age))

            # Move back into iso directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return

#---------------------------------------#
# Now for the Pisa (Tognelli+11) models
#---------------------------------------#

class Pisa(StellarEvolution):
    """
    Evolution models from
    `Tognelli et al. 2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...533A.109T/abstract>`_.

    Downloaded `online <http://astro.df.unipi.it/stellar-models/index.php?m=1>`_

    Notes
    ------
    Parameters used in download:

    * Y = middle value of 3 provided (changes for different metallicities)
    * mixing length = 1.68
    * Deuterium fraction: 2*10^-5 for Z = 0.015, 0.03; 4*10^-4 for 0.005
    """
    def __init__(self):
        r"""
        Define intrinsic properties for the Pisa (Tognelli+11) stellar
        models.
        """
        self.model_version_name = "Pisa"
        # define metallicity parameters for Pisa models
        self.z_list = [0.015]

        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.0, 8.01+0.005, 0.01)

        # Specify location of model files
        self.model_dir = models_dir+'Pisa2011/'

        # Specifying metallicity
        self.z_solar = 0.015
        self.z_file_map = {0.015: 'z015/'}

        # Define required evo_grid number
        self.evo_grid_min = 1.0

    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Pisa (Tognelli+11)
        collection.
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)

        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')
            return
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds for evolution model. Available z-vals: {self.z_list}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file

        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col1', 'logL')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'mass')
        iso.rename_column('col4', 'logg')
        iso['logT_WR'] = iso['logT']

        # Pisa models are too low for WR phase, add WR column with all False
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)

        # Add columns for current mass and phase.
        iso.add_column( Column(np.zeros(len(iso)), name = 'phase'))
        iso.add_column( Column(iso['mass'], name = 'mass_current'))

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

    def format_isochrones(input_iso_dir, metallicity_list):
        r"""
        Rename the isochrone files extracted from Pisa (Tognelli+11) to fit
        naming/directory scheme

        input_iso_dir: points to Pisa2011/iso directory. Individual
        metallicity directories with the downloaded isochrones are
        expected to already exist there

        metallicity_list is the list of metallicities on which function
        is to be run.

        format for metallicity_list : absolute (vs. relative to sun)
        'z' + <digits after decimal>, e.g Z = 0.015 --> z015.
        """
        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)
        # Work on each metallicity directory individually
        for metal in metallicity_list:
            # Move into directory, check to see if files are already formatted
            os.chdir(metal)

            if os.path.exists('iso_6.00.fits'):
                print( 'Files in {0:s} already formatted'.format(metal))
            else:
                # Create a ReadMe with the original file names to preserve the
                # model details

                cmd = "ls *.FITS > ReadMe"
                os.system(cmd)

                # Collect all filenames in a list, rename files one
                # by one
                isoFile_list = glob.glob('*.FITS')
                for File in isoFile_list:
                    name = File.split('_')
                    # Extract iso age from filename
                    age = float(name[1][1:])
                    logAge = np.log10(age * 10**6)

                    cmd = "mv {0:s} iso_{1:4.2f}.fits".format(File, logAge)
                    os.system(cmd)

            # Return to overhead directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return

    def make_isochrone_grid(metallicity=0.015):
        """
        Create isochrone grid of given metallicity with time sampling = 0.01
        in logAge (hardcoded). This interpolates the downloaded isochrones
        when necessary. Builds upon the online iscohrone grid.

        Note: format of metallicity is important. After decimal point, must match
        the format of the metallcity directory (i.e., 0.015 matches directory z015,
        while 0.0150 would not)
        """
        logAge_arr = np.arange(6.0, 8.0+0.005, 0.01)

        count = 0
        for logAge in logAge_arr:
            # Could interpolate using evolutionary tracks, but less accurate.
            make_isochrone_pisa_interp(logAge, metallicity=metallicity)

            count += 1

            print( 'Done {0} of {1} models'.format(count, (len(logAge_arr))))

        return

#==============================#
# Baraffe+15 models
#==============================#
class Baraffe15(StellarEvolution):
    """
    Evolution models published in
    `Baraffe et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract>`_.

    Downloaded from `BHAC15 site <http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks>`_.
    """
    def __init__(self):
        self.model_version_name = "Baraffe15"
        # define metallicity parameters for Baraffe models
        self.z_list = [0.015]

        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(6.0, 8.0+0.005, 0.01)

        # Specify location of model files
        self.model_dir = models_dir+'Baraffe15/'

        # Specifying metallicity
        self.z_solar = 0.015
        self.z_file_map = {0.015: 'z015/'}

        # Define required evo_grid number
        self.evo_grid_min = 1.0

    def isochrone(self, age=5.e7, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Baraffe+15
        collection.
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)

        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file

        # Read isochrone, get in proper format
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('Mass', 'mass')
        iso.rename_column('logG', 'logg')
        iso['logT'] = np.log10(iso['Teff'])

        # Pisa models are too low for WR phase, add WR column with all False
        iso['logT_WR'] = iso['logT']
        isWR = Column([False] * len(iso), name='isWR')
        iso.add_column(isWR)

        # Add columns for current mass and phase.
        iso.add_column( Column(np.zeros(len(iso)), name = 'phase'))
        iso.add_column( Column(iso['mass'], name = 'mass_current'))

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

    def tracks_to_isochrones(self, tracksFile):
        r"""
        Create isochrones at desired age sampling (6.0 < logAge < 8.0,
        steps of 0.01; hardcoded) from the Baraffe+15 tracks downloaded
        online.

        tracksFile: tracks.dat file downloaded from Baraffe+15, with format
        modified to be read in python

        Writes isochrones in iso/ subdirectory off of work directory. Will
        create this subdirectory if it doesn't already exist
        """
        tracks = Table.read(tracksFile, format='ascii')

        age_arr = np.arange(6.0, 8.0+0.005, 0.01)
        #age_arr = [6.28]

        # Loop through the masses, interpolating track over time at each.
        # Resample track properties at hardcoded ages
        masses = np.unique(tracks['col1'])

        mass_interp = []
        age_interp = []
        Teff_interp = []
        logL_interp = []
        logG_interp = []
        print( 'Begin looping over masses')
        cnt=0
        for mass in masses:
            idx = np.where(tracks['col1'] == mass)
            tmp = tracks[idx]

            # First, extract Teff, logL, and logG, eliminating
            # duplicated inputs (these crash the interpolator)
            good_Teff = np.where( np.diff(tmp['col3']) != 0 )
            good_logG = np.where( np.diff(tmp['col5']) != 0 )
            good_logL = np.where( np.diff(tmp['col4']) != 0 )

            # Interpolate Teff, logL, and logG using linear interpolator
            tck_Teff = interpolate.interp1d(tmp['col2'], tmp['col3'])
            tck_logL = interpolate.interp1d(tmp['col2'], tmp['col4'])
            tck_logG = interpolate.interp1d(tmp['col2'], tmp['col5'])


            Teff = tck_Teff(age_arr)
            logL = tck_logL(age_arr)
            logG = tck_logG(age_arr)

            # Test interpolation if desired
            test=False
            if test:
                py.figure(1, figsize=(10,10))
                py.clf()
                py.plot(tmp['col2'], tmp['col3'], 'k.', ms=8)
                py.plot(age_arr, Teff, 'r-', linewidth=2)
                py.xlabel('logAge')
                py.ylabel('Teff')
                py.savefig('test_Teff.png')

                py.figure(2, figsize=(10,10))
                py.clf()
                py.plot(tmp['col2'], tmp['col4'], 'k.', ms=8)
                py.plot(age_arr, logL, 'r-', linewidth=2)
                py.xlabel('logAge')
                py.ylabel('logL')
                py.savefig('test_logL.png')

                py.figure(3, figsize=(10,10))
                py.clf()
                py.plot(tmp['col2'], tmp['col5'], 'k.', ms=8)
                py.plot(age_arr, logG, 'r-', linewidth=2)
                py.xlabel('logAge')
                py.ylabel('logG')
                py.savefig('test_logG.png')

                pdb.set_trace()

            # Build upon arrays of interpolated values
            mass_interp = np.concatenate((mass_interp, np.ones(len(Teff)) * mass))
            age_interp = np.concatenate((age_interp, age_arr))
            Teff_interp = np.concatenate((Teff_interp, Teff))
            logL_interp = np.concatenate((logL_interp, logL))
            logG_interp = np.concatenate((logG_interp, logG))

            print( 'Done {0} of {1}'.format(cnt, len(masses)))
            cnt+=1

        # Now, construct the iso_*.fits files for each age, write files
        # to iso subdirectory
        # First check to see if subdirectory exists
        if not os.path.exists('iso/'):
            os.mkdir('iso')

        # Now for the loop
        ages = np.unique(age_interp)
        print( 'Writing iso files')
        for age in ages:
            good = np.where( age_interp == age)

            t = Table( (mass_interp[good], Teff_interp[good], logL_interp[good],
                       logG_interp[good]), names=('Mass', 'Teff', 'logL', 'logG') )

            # Write out as fits table
            name = 'iso_{0:3.2f}.fits'.format(age)
            t.write('iso/'+name, format='fits', overwrite=True)

        return

    def test_age_interp(self, onlineIso, interpIso):
        r"""
        Compare one of our interpolated ischrones with one
        of the isochrones provided online by Baraffe+15.
        """
        true_iso = Table.read(onlineIso, format='ascii')
        our_iso = Table.read(interpIso, format='fits')

        # Compare the two isochrones using plots. Look at mass vs. Teff,
        # mass vs. logG, mass vs. logL. Ideally these isochrones should
        # be identical
        py.figure(1, figsize=(10,10))
        py.clf()
        py.plot(true_iso['col1'], true_iso['col2'], 'k.', ms = 10)
        py.plot(our_iso['Mass'], our_iso['Teff'], 'r.', ms = 10)
        py.xlabel('Mass')
        py.ylabel('Teff')
        py.savefig('interp_test1.png')

        py.figure(2, figsize=(10,10))
        py.clf()
        py.plot(true_iso['col1'], true_iso['col3'], 'k.', ms = 10)
        py.plot(our_iso['Mass'], our_iso['logL'], 'r.', ms = 10)
        py.xlabel('Mass')
        py.ylabel('logL')
        py.savefig('interp_test2.png')

        py.figure(3, figsize=(10,10))
        py.clf()
        py.plot(true_iso['col1'], true_iso['col4'], 'k.', ms = 10)
        py.plot(our_iso['Mass'], our_iso['logG'], 'r.', ms = 10)
        py.xlabel('Mass')
        py.ylabel('logG')
        py.savefig('interp_test3.png')

        # Look at the difference between values (assumes the masses are lined up)
        Teff_diff = np.mean(abs(true_iso['col2'][7:] - our_iso['Teff']))
        logL_diff = np.mean(abs(true_iso['col3'][7:] - our_iso['logL']))
        logG_diff = np.mean(abs(true_iso['col4'][7:] - our_iso['logG']))

        print( 'Average abs difference in Teff: {0}'.format(Teff_diff))
        print( 'Average abs difference in logL: {0}'.format(logL_diff))
        print( 'Average abs difference in logg: {0}'.format(logG_diff))

        return

def compare_Baraffe_Pisa(BaraffeIso, PisaIso):
    """
    Compare the Baraffe isochrones to the Pisa isochrones, since they overlap
    over a significant portion of mass space.
    """
    b = Table.read(BaraffeIso, format='fits')
    p = Table.read(PisaIso, format='ascii')

    name = BaraffeIso.split('_')
    age = name[1][:4]

    # Extract paramters we need
    b_mass = b['Mass']
    b_logT = np.log10(b['Teff'])
    b_logL = b['logL']
    b_logG = b['logG']

    p_mass = p['col3']
    p_logT = p['col2']
    p_logL = p['col1']
    p_logG = p['col4']

    m05_b = np.where( abs(b_mass - 0.5) == min(abs(b_mass - 0.5)) )
    m05_p = np.where( abs(p_mass - 0.5) == min(abs(p_mass - 0.5)) )

    # Comparison plots
    py.figure(1, figsize=(10,10))
    py.clf()
    py.plot(b_logT, b_logL, 'k-', linewidth=2, label='Baraffe+15')
    py.plot(b_logT[m05_b], b_logL[m05_b], 'k.', ms=10)
    py.plot(p_logT, p_logL, 'r', linewidth=2, label='Pisa')
    py.plot(p_logT[m05_p], p_logL[m05_p], 'r.', ms=10)
    py.xlabel('logT')
    py.ylabel('logL')
    py.title(age)
    py.axis([4.4, 3.4, -3, 4])
    #py.gca().invert_xaxis()
    py.legend()
    py.savefig('BaraffePisa_comp_{0}.png'.format(age))

    py.figure(2, figsize=(10,10))
    py.clf()
    py.plot(b_mass, b_logL, 'k-', linewidth=2, label='Baraffe+15')
    py.plot(b_mass[m05_b], b_logL[m05_b], 'k.', ms=10)
    py.plot(p_mass, p_logL, 'r', linewidth=2, label='Pisa')
    py.plot(p_mass[m05_p], p_logL[m05_p], 'r.', ms=10)
    py.xlabel('Mass')
    py.ylabel('logL')
    py.title(age)
    #py.axis([4.4, 3.4, -3, 4])
    #py.gca().invert_xaxis()
    py.legend()
    py.savefig('BaraffePisa_comp_mass_{0}.png'.format(age))

    return

#===============================#
# MIST v.1 (Choi+16)
#===============================#
class MISTv1(StellarEvolution):
    """
    Define intrinsic properties for the MIST v1 stellar
    models.

    Models originally downloaded from `online server <http://waps.cfa.harvard.edu/MIST/interp_isos.html>`_.

    Parameters
    ----------
    version: '1.0' or '1.2', optional
        Specify which version of MIST models you want. Version 1.0
        was downloaded from MIST website on 2/2017, while Version 1.2
        was downloaded on 8/2018 (solar metallicity)
        and 4/2019 (other metallicities). Default is 1.2.

    synthpop_extension: boolean (default True)
        If True, the isochrones are extended down to a minimum initial
        mass of 0.1Msun using grids interpolated via `SynthPop
        <https://ui.adsabs.harvard.edu/abs/2025AJ....169..317K/abstract>`_.
        If False, the web-downloaded MIST isochrones are used with their
        varying lower mass limits. True option is only valid for version=1.2.
    """
    def __init__(self, version=1.2, synthpop_extension=True):
        # define metallicity parameters for MIST models
        self.z_list = [0.0000014,   # [Fe/H] = -4.00
                       0.0000045,   # [Fe/H] = -3.50
                       0.000014,    # [Fe/H] = -3.00
                       0.000045,    # [Fe/H] = -2.50
                       0.00014, # [Fe/H] = -2.00
                       0.00025, # [Fe/H] = -1.75
                       0.00045, # [Fe/H] = -1.50
                       0.00080, # [Fe/H] = -1.25
                       0.0014,  # [Fe/H] = -1.00
                       0.0025,  # [Fe/H] = -0.75
                       0.0045,  # [Fe/H] = -0.50
                       0.0080,  # [Fe/H] = -0.25
                       0.014,   # [Fe/H] = 0.00
                       0.025,   # [Fe/H] = 0.25
                       0.045]   # [Fe/H] = 0.50

        # populate list of isochrone ages (log scale)
        self.age_list = np.arange(5.01, 10.30+0.005, 0.01)

        # Set version directory
        self.version = version
        self.synthpop_extension = synthpop_extension
        if self.version == 1.0:
            self.model_version_name = 'MISTv1.0'
            version_dir = 'v1.0/'
            if synthpop_extension:
                warnings.warn('SynthPop isochrone extension not supported for MISTv1.0 isochrones')
                self.synthpop_extension = False
        elif self.version == 1.2:
            self.model_version_name = 'MISTv1.2'
            version_dir = 'v1.2/'
        else:
            raise ValueError('Version {0} not supported for MIST isochrones'.format(version))

        # Specify location of model files
        self.model_dir = models_dir+'MISTv1/' + version_dir
        if self.synthpop_extension:
            self.model_version_name = self.model_version_name + '-synthpop'
            self.model_extension_dir = models_dir+'MISTv1/' + version_dir[:-1] + '-synthpop/'
            if not os.path.exists(self.model_extension_dir):
                raise ValueError(f'Missing {self.model_extension_dir}. Please download the latest SPISEA data at https://w.astro.berkeley.edu/~jlu/spisea/spisea_models.tar.gz.')

        else:
            self.model_extension_dir = None

        # Specifying metallicity
        self.z_solar = 0.0142
        self.z_file_map = {0.0000014: 'z0000014/',
                           0.0000045: 'z0000045/',
                           0.000014: 'z000014/',
                           0.000045: 'z000045/',
                           0.00014: 'z00014/',
                           0.00025: 'z00025/',
                           0.00045: 'z00045/',
                           0.00080: 'z00080/',
                           0.0014: 'z0014/',
                           0.0025: 'z0025/',
                           0.0045: 'z0045/',
                           0.0080: 'z0080/',
                           0.014: 'z014/',
                           0.025: 'z025/',
                           0.045: 'z045/'}

        # Define required evo_grid number (now 1.2 for synthpop extension)
        self.evo_grid_min = 1.2

    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the MISTv1
        collection.
        """
        # First, error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)

        # convert metallicity to mass fraction
        z_defined = self.z_solar * (10.**metallicity)

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if ((z_defined < np.min(self.z_list)-0.1) or
                (z_defined > np.max(self.z_list)+0.1)):
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + 'iso/' + z_dir + iso_file
        if self.synthpop_extension:
            addl_iso_file = self.model_extension_dir + 'iso/' + z_dir + iso_file

        # return isochrone data. Column locations depend on
        # version
        iso = Table.read(full_iso_file, format='fits')
        if self.synthpop_extension:
            addl_iso = Table.read(addl_iso_file, format='fits')
            for row in addl_iso:
                iso.add_row(row)
            iso.sort('col3')
            iso.meta['comments2'] = addl_iso.meta['comments']
        if self.version == 1.0:
            iso.rename_column('col7', 'Z')
            iso.rename_column('col2', 'logAge')
            iso.rename_column('col3', 'mass')
            iso.rename_column('col4', 'logT')
            iso.rename_column('col5', 'logg')
            iso.rename_column('col6', 'logL')
            iso.rename_column('col65', 'phase')
        elif self.version == 1.2:
            iso.rename_column('col2', 'logAge')
            iso.rename_column('col3', 'mass')
            iso.rename_column('col4', 'mass_current')
            iso.rename_column('col9', 'logL')
            iso.rename_column('col14', 'logT')
            iso.rename_column('col17', 'logg')
            iso.rename_column('col79', 'phase')

        # For MIST isochrones, anything with phase = 6 is a WD.
        # Following our IFMR convention, change the phase designation
        # to 101
        isWD = np.where(iso['phase'] == 6)[0]
        iso['phase'][isWD] = 101

        # Define "isWR" column based on phase info
        isWR = Column([False] * len(iso), name='isWR')
        idx_WR = np.where(iso['phase'] == 9)[0]
        isWR[idx_WR] = True
        iso.add_column(isWR)

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

    def format_isochrones(self):
        r"""
        Parse isochrone file downloaded from MIST web server,
        create individual isochrone files for the different ages.
        Assumes all files start with MIST_iso*

        Parameters:
        -----------
        input_iso_dir: path
            Points to MISTv1/<version>/iso directory.

        metallicity_list: array
            List of metallicity directories to check (i.e. z015 is solar)

        """
        # Get input iso dir, metallicity list from evo object
        input_iso_dir = '{0}/iso'.format(self.model_dir)
        metallicity_list = list(self.z_file_map.values())

        # Store current directory for later
        start_dir = os.getcwd()

        # Move into isochrone directory
        os.chdir(input_iso_dir)

        # Work on each metallicity isochrones individually
        for metal in metallicity_list:
            # More into metallicity directory, read isochrone file
            os.chdir(metal)

            # Read all available iso files, stack them together
            isoFile = glob.glob('MIST_iso*')
            print( 'Read Input: this is slow')
            iso_f = Table()
            for ii in isoFile:
                tmp = Table.read(ii, format='ascii')
                iso_f = vstack([iso_f, tmp])
            print( 'Done')

            # Extract the unique ages
            ages_all = iso_f['col2']
            age_arr = np.unique(ages_all)

            # For each unique age, extract the proper rows and make corresponding
            # table
            print( 'Making individual isochrone files')
            for age in age_arr:
                good = np.where(ages_all == age)
                tmp = iso_f[good]

                # Need to make sure the tables are unmasked...this causes
                # problems later
                tmp2 = Table(tmp, masked=False)

                #Write table
                tmp2.write('iso_{0:4.2f}.fits'.format(age))

            # Move back into iso directory
            os.chdir('..')

        # Return to starting directory
        os.chdir(start_dir)
        return
        
#===========================================#
# COSMIC Breivik+ 2020 - not normal evo model
#===========================================#
class COSMIC(StellarEvolution):
    """
    Evolve objects using COSMIC
    `Breivik et al. 2020 <hhttps://ui.adsabs.harvard.edu/abs/2019arXiv191100903B/abstract>`_.

    See code website here: <https://cosmic-popsynth.github.io/>`_.

    Parameters
    ----------
    BSEDict: dict or string, optional
        Binary Stellar Evolution dictionary for COSMIC evolution. 
        Default is 'default' which uses the dictionary from COSMIC docs with zsun = 0.02.

    keep_disrupted_companions: bool, optional
        When True, if the system is disrupted, the companions are added to the primary table.
        When False, the companion is deleted.
        Default is True.

    keep_COSMIC_tables : bool, optional
        Allows COSMIC tables to be accessible on the SPISEA evolution object.
        Default is False.
    """
    def __init__(self, BSEDict='default', keep_disrupted_companions=True, keep_COSMIC_tables=False): 
        if BSEDict == 'default':
            self.BSEDict = {
                                "pts1": 0.001, "pts2": 0.01, "pts3": 0.02, "zsun": 0.02, "windflag": 3,
                                "eddlimflag": 0, "neta": 0.5, "bwind": 0.0, "hewind": 0.5, "beta": 0.125,
                                "xi": 0.5, "acc2": 1.5, "LBV_flag": 1, "alpha1": [1.0, 1.0],
                                "lambdaf": 0.0, "ceflag": 1, "cekickflag": 2, "cemergeflag": 1,
                                "cehestarflag": 0, "qcflag": 5,
                                "qcrit_array": [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                                "kickflag": 5, "sigma": 265.0, "bhflag": 1, "bhsigmafrac": 1.0,
                                "sigmadiv": -20.0, "ecsn": 2.25, "ecsn_mlow": 1.6, "aic": 1, "ussn": 1,
                                "polar_kick_angle": 90.0,
                                "natal_kick_array": [[-100.0, -100.0, -100.0, -100.0, 0.0], [-100.0, -100.0, -100.0, -100.0, 0.0]],
                                "mm_mu_ns": 400.0, "mm_mu_bh": 200.0, "remnantflag": 4,
                                "fryer_mass_limit": 0, "mxns": 3.0, "fryer_fmix": 1.0,
                                "fryer_mcrit_nsbh": 5.75, "rembar_massloss": 0.5, "wd_mass_lim": 1,
                                "maltsev_mode": 0, "maltsev_fallback": 0.5, "maltsev_pf_prob": 0.1,
                                "pisn": -2, "ppi_co_shift": 0.0, "ppi_extra_ml": 0.0, "bhspinflag": 0,
                                "bhspinmag": 0.0, "grflag": 1, "eddfac": 10, "gamma": -2, "don_lim": -1,
                                "acc_lim": [-1, -1], "smt_periastron_check": 0, "tflag": 1, "ST_tide": 1,
                                "fprimc_array": [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0],
                                "ifflag": 1, "wdflag": 1, "epsnov": 0.001, "bdecayfac": 1,
                                "bconst": 3000, "ck": 1000, "rejuv_fac": 1.0, "rejuvflag": 0,
                                "bhms_coll_flag": 0, "htpmb": 1, "ST_cr": 1, "rtmsflag": 0
                            }

        else:
            self.BSEDict = BSEDict

        self.external_evol = True
        self.z_solar = 0.02 #0.014
        self.keep_disrupted_companions = keep_disrupted_companions
        self.keep_COSMIC_tables = keep_COSMIC_tables
        self.model_version_name = "COSMIC"


    def evolve(self, star_systems, companions, logAge, metallicity):
        from cosmic.utils import p_from_a, a_from_p
        from cosmic.sample.initialbinarytable import InitialBinaryTable
        from cosmic.evolve import Evolve
        
        companion_system_idxs = companions['system_idx']
        m1s = star_systems['mass']
        m2s = np.zeros(len(star_systems))
        m2s[companion_system_idxs] = companions['mass']
        
        a_Rsuns = (10**companions['log_a'])*u.AU.to('Rsun')
        porbs = np.zeros(len(star_systems))
        porbs[companion_system_idxs] = p_from_a(a_Rsuns, m1s[companion_system_idxs], m2s[companion_system_idxs])
        
        eccs = np.zeros(len(star_systems))
        eccs[companion_system_idxs] = companions['e']
        
        kstar1s = (m1s >= 0.7).astype(int) # 1 if MS above 0.7 and 0 if MS below 0.7
        kstar2s = (m2s >= 0.7).astype(int) # 1 if MS above 0.7 and 0 if MS below 0.7

        binary_pop = InitialBinaryTable.InitialBinaries(m1=m1s, m2=m2s, porb=porbs,
                                                   ecc=eccs, tphysf=[10**logAge/1e6]*len(m1s),
                                                   kstar1=kstar1s, kstar2=kstar2s, metallicity=[self.z_solar*10**metallicity]*len(m1s))
        
        bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary_pop, BSEDict=self.BSEDict)

        final_binaries = bcm[bcm['tphys'] > 0] #only gives the first and last idx, so this takes final one
        
        # Add number for system idx since we're about to manipuluate them a bunch
        star_systems['system_idx'] = np.arange(len(star_systems))
        
        # Remove systems that don't show up in bcm final (very few)
        if len(final_binaries) != len(star_systems):
            # Save or append to initC
            initC_fail_path = 'initC_fail.csv'
            exists = os.path.exists(initC_fail_path)
            initC.to_csv(
                    initC_fail_path,
                    mode='a' if exists else 'w',
                    header=not exists,
                    index=False
                    )

            # Save or append failing binary bcm values to table
            mask = ~bcm.loc[bcm['tphys'] == 0, 'bin_num'].isin(final_binaries['bin_num'])
            failing_binaries = bcm.loc[bcm['tphys'] == 0].loc[mask]
            
            print('missing binaries', bcm.loc[bcm['tphys'] == 0].loc[mask])
            missing_binary_csv_path = 'missing_cosmic_binaries_bcm.csv'
            exists = os.path.exists( missing_binary_csv_path)
            failing_binaries.to_csv(
                    missing_binary_csv_path,
                    mode='a' if exists else 'w',
                    header=not exists,
                    index=False
                    )

            # Remove failing binaries from all tables and redefine quantites
            bad_bin_nums = failing_binaries['bin_num']
            bad_mask_ss = np.isin(star_systems['system_idx'], bad_bin_nums)
            star_systems.remove_rows(np.where(bad_mask_ss)[0])
            
            bad_mask_comp = np.isin(companions['system_idx'], bad_bin_nums)
            companions.remove_rows(np.where(bad_mask_comp)[0])

            bpp = bpp[~bpp['bin_num'].isin(bad_bin_nums)]
            bcm = bcm[~bcm['bin_num'].isin(bad_bin_nums)]
            kick_info = kick_info[~kick_info['bin_num'].isin(bad_bin_nums)]

            # Redefine system idx since companions refer to the positions
            star_systems['old_system_idx'] = star_systems['system_idx']
            star_systems['system_idx'] = np.arange(len(star_systems))
            idx_map = dict(zip(star_systems['old_system_idx'], star_systems['system_idx']))
            companions['system_idx'] = np.array([idx_map[i] for i in companions['system_idx']])

            # redefine the bin_num to the new values too
            bpp['bin_num'] = bpp['bin_num'].map(idx_map)
            bcm['bin_num'] = bcm['bin_num'].map(idx_map)
            kick_info['bin_num'] = kick_info['bin_num'].map(idx_map)

            assert(set(companions['system_idx']) - set(star_systems['system_idx']) == set())
            assert(set(bpp['bin_num']) - set(star_systems['system_idx']) == set())
            assert(set(bcm['bin_num']) - set(star_systems['system_idx']) == set())
            assert(set(kick_info['bin_num']) - set(star_systems['system_idx']) == set())

            # reset the index to the bin_num column
            bpp = bpp.set_index('bin_num', drop=False)
            bcm = bcm.set_index('bin_num', drop=False)
            kick_info = kick_info.set_index('bin_num', drop=False)
            
            companion_system_idxs = companions['system_idx']
            m1s = star_systems['mass']
            m2s = np.zeros(len(star_systems))
            m2s[companion_system_idxs] = companions['mass']
            
            a_Rsuns = (10**companions['log_a'])*u.AU.to('Rsun')
            porbs = np.zeros(len(star_systems))
            porbs[companion_system_idxs] = p_from_a(a_Rsuns, m1s[companion_system_idxs], m2s[companion_system_idxs])
            
            eccs = np.zeros(len(star_systems))
            eccs[companion_system_idxs] = companions['e']
            
            kstar1s = (m1s >= 0.7).astype(int) # 1 if MS above 0.7 and 0 if MS below 0.7
            kstar2s = (m2s >= 0.7).astype(int) # 1 if MS above 0.7 and 0 if MS below 0.7

            final_binaries = bcm[bcm['tphys'] > 0] #only gives the first and last idx, so this takes final one
            
            print("WARNING: Some binaries didn't make it. Something went wrong with COSMIC. Saved to initC_fail.csv and missing_cosmic_binaries_bcm.csv")
        
        # initializes kick columns with zeros
        star_systems['kick_x'] = 0
        star_systems['kick_y'] = 0
        star_systems['kick_z'] = 0

        # rotates output kicks from cosmic to inclination
        # picks random direction for singles
        inclinations = np.arccos(2 * np.random.rand(len(star_systems)) - 1.0) # initializes random inclinations in radians
        existing_inclinations = np.deg2rad(companions['i'])
        inclinations[companion_system_idxs] = existing_inclinations 
        inclinations = np.repeat(inclinations, 2) # accounts for two rows per system in kick table
        rotated_kick_values_1 = self.get_kick_differential(kick_info['delta_vsysx_1'], kick_info['delta_vsysy_1'], 
                                                         kick_info['delta_vsysz_1'], inclination = inclinations)
        rotated_kick_values_2 = self.get_kick_differential(kick_info['delta_vsysx_2'], kick_info['delta_vsysy_2'], 
                                                         kick_info['delta_vsysz_2'], inclination = inclinations)
        kick_info['delta_vsysx_1_rot'] = rotated_kick_values_1.d_x
        kick_info['delta_vsysy_1_rot'] = rotated_kick_values_1.d_y
        kick_info['delta_vsysz_1_rot'] = rotated_kick_values_1.d_z
        
        kick_info['delta_vsysx_2_rot'] = rotated_kick_values_2.d_x
        kick_info['delta_vsysy_2_rot'] = rotated_kick_values_2.d_y
        kick_info['delta_vsysz_2_rot'] = rotated_kick_values_2.d_z
        
        
        star_systems['mass_current'] = final_binaries['mass_1']
        star_systems['Teff'] = final_binaries['teff_1']
        star_systems['L'] = final_binaries['lum_1']
        star_systems['logg'] = self.calc_logg(final_binaries['mass_1'], final_binaries['rad_1'])
        
        # Takes sum of the delta kicks in case there was a kick, no disruption, then second kick
        # Even for isolated stars, take sum since second row is blank
        primary_kick_sum = (
                kick_info
                .groupby(level=0)[["delta_vsysx_1_rot", "delta_vsysy_1_rot", "delta_vsysz_1_rot"]]
                .sum()
            )
        star_systems["kick_x"] = primary_kick_sum["delta_vsysx_1_rot"].reindex(star_systems["system_idx"], fill_value=0).to_numpy()
        star_systems["kick_y"] = primary_kick_sum["delta_vsysy_1_rot"].reindex(star_systems["system_idx"], fill_value=0).to_numpy()
        star_systems["kick_z"] = primary_kick_sum["delta_vsysz_1_rot"].reindex(star_systems["system_idx"], fill_value=0).to_numpy()
        
        companions['mass_current'] = final_binaries['mass_2'][companion_system_idxs]
        companions['Teff'] = final_binaries['teff_2'][companion_system_idxs]
        companions['L'] = final_binaries['lum_2'][companion_system_idxs]
        companions['logg'] = self.calc_logg(final_binaries['mass_2'][companion_system_idxs], final_binaries['rad_2'][companion_system_idxs])
        # Also take sum of companion kicks
        companion_kick_sum = (
                kick_info
                .groupby(level=0)[["delta_vsysx_2_rot", "delta_vsysy_2_rot", "delta_vsysz_2_rot"]]
                .sum()
            )
        companions['kick_x'] = companion_kick_sum["delta_vsysx_2_rot"].reindex(companion_system_idxs, fill_value=0).to_numpy()
        companions['kick_y'] = companion_kick_sum["delta_vsysy_2_rot"].reindex(companion_system_idxs, fill_value=0).to_numpy()
        companions['kick_z'] = companion_kick_sum["delta_vsysz_2_rot"].reindex(companion_system_idxs, fill_value=0).to_numpy()
        
        loga = np.log10(final_binaries['sep'][companion_system_idxs]*u.Rsun.to('AU'))
        companions['log_a'] = loga
        
        fixed_phases1 = final_binaries['kstar_1'].to_numpy(copy=True)
        fixed_phases1[np.where((final_binaries['kstar_1'] >= 10) & (final_binaries['kstar_1'] <= 12))[0]] = 101
        fixed_phases1[np.where(final_binaries['kstar_1'] == 13)[0]] = 102
        fixed_phases1[np.where(final_binaries['kstar_1'] == 14)[0]] = 103
        star_systems['phase'] = fixed_phases1
        
        fixed_phases2 = final_binaries['kstar_2'][companion_system_idxs].to_numpy(copy=True)
        fixed_phases2[np.where((final_binaries['kstar_2'][companion_system_idxs] >= 10) & (final_binaries['kstar_2'][companion_system_idxs] <= 12))[0]] = 101
        fixed_phases2[np.where(final_binaries['kstar_2'][companion_system_idxs] == 13)[0]] = 102
        fixed_phases2[np.where(final_binaries['kstar_2'][companion_system_idxs] == 14)[0]] = 103
        companions['phase'] = fixed_phases2
        
        # maybe add WR designation


        # Take the disrupted binaries and put the companions into the star_system table (if desired)
        # don't include massless remnant companions
        disrupted_binary_companion_idxs = np.where((final_binaries['bin_state'][companion_system_idxs] == 2) & (final_binaries['kstar_2'][companion_system_idxs] != 15))[0]
        disrupted_binary_companions_num = 0
        if self.keep_disrupted_companions and len(disrupted_binary_companion_idxs) > 0:
            disrupted_binary_companions = companions[disrupted_binary_companion_idxs]
            disrupted_binary_companions['systemMass'] = disrupted_binary_companions['mass_current']
            disrupted_binary_companions['isMultiple'] = [False]*len(disrupted_binary_companions)
            disrupted_binary_companions['N_companions'] = [0]*len(disrupted_binary_companions)
            disrupted_binary_companions.remove_columns(['system_idx', 'log_a', 'e', 'i', 'Omega', 'omega'])
            disrupted_binary_companions['system_idx'] = np.arange(len(disrupted_binary_companions)) + len(star_systems)
            star_systems = vstack([star_systems, disrupted_binary_companions])
            disrupted_binary_companions_num = len(disrupted_binary_companions)

        #Drop merged companions and totally disappeared systems
        # Also promote companions to primaries when the initial primary "merged" (if desired)
        mr_companion_only_idxs = np.where((final_binaries['kstar_1'][companion_system_idxs] != 15) & (final_binaries['kstar_2'][companion_system_idxs] == 15))[0] #mr for massless remnant
        disappeared_system_companion_idxs = np.where((final_binaries['kstar_1'][companion_system_idxs] == 15) & (final_binaries['kstar_2'][companion_system_idxs] == 15))[0]
        companions_to_mr_primaries_idxs = np.where((final_binaries['kstar_1'][companion_system_idxs] == 15) & (final_binaries['kstar_2'][companion_system_idxs] != 15))[0]
        
        mr_primary_only_idx = np.where((final_binaries['kstar_1'] == 15) & (final_binaries['kstar_2'] != 15))[0] #mr for massless remnant
        disappeared_system_primaries = np.where((final_binaries['kstar_1'] == 15) & (final_binaries['kstar_2'] == 15))[0]
        
        
        delete_primary_idxs = np.concatenate((mr_primary_only_idx, disappeared_system_primaries))
        delete_companion_idxs = np.concatenate((disrupted_binary_companion_idxs, mr_companion_only_idxs, disappeared_system_companion_idxs, companions_to_mr_primaries_idxs))
        primaries_to_deleted_companion_idxs = companions[delete_companion_idxs]['system_idx']
        
        #Fix binary specification of primaries that lost their companions
        #star_systems['isMultiple'][primaries_to_deleted_companion_idxs] = False
        #star_systems['N_companions'][primaries_to_deleted_companion_idxs] = 0

        mask = np.isin(star_systems['system_idx'],
               primaries_to_deleted_companion_idxs)

        star_systems['N_companions'][mask] = 0
        star_systems['isMultiple'][mask] = False
        
        # Promote the companions to merged primaries to primaries
        if self.keep_disrupted_companions and len(companions_to_mr_primaries_idxs) > 0:
            companions_to_mr_primaries = companions[companions_to_mr_primaries_idxs]
            companions_to_mr_primaries['systemMass'] = companions_to_mr_primaries['mass_current']
            companions_to_mr_primaries['isMultiple'] = [False]*len(companions_to_mr_primaries)
            companions_to_mr_primaries['N_companions'] = [0]*len(companions_to_mr_primaries)
            companions_to_mr_primaries.remove_columns(['system_idx', 'log_a', 'e', 'i', 'Omega', 'omega'])
            companions_to_mr_primaries['system_idx'] = np.arange(len(companions_to_mr_primaries)) + len(star_systems) + disrupted_binary_companions_num
            star_systems = vstack([star_systems, companions_to_mr_primaries])
        
        star_systems.remove_rows(delete_primary_idxs)
        companions.remove_rows(delete_companion_idxs) #if kstar 1 is 15 take the seocnd star and if kstar2 is 15 take the other

        #Reassign system_idx vals
        star_systems['system_idx_new'] = np.arange(len(star_systems))
        mapping = np.empty(star_systems['system_idx'].max() + 1, dtype=int)
        mapping[star_systems['system_idx']] = star_systems['system_idx_new']
        companions['system_idx'] = mapping[companions['system_idx']]
        star_systems.remove_columns(['system_idx', 'system_idx_new'])

        # Preserve a scalar kick magnitude alongside the vector components.
        for table in (star_systems, companions):
            table['kick'] = np.sqrt(table['kick_x']**2 + table['kick_y']**2 + table['kick_z']**2)

        # Make sure we didn't break anything by manipulating the number of companions
        assert star_systems['N_companions'].sum() == len(companions)
        
        if self.keep_COSMIC_tables:
            self.bpp = bpp
            self.bcm = bcm
            self.initC = initC
            self.kick_info = kick_info

        return star_systems, companions

        
    def calc_logg(self, masses, radii):
        """
        Inputs
        ------
        masses : array-like 
            Masses of objects in Msun
            
        radii : array-like
            Radii of stars in Rsun
    
        Returns
        -------
        logg : array-like
            Log10 surface gravity in cgs
        """
        return np.log10(((np.array(c.G.to('Rsun^3/(Msun*s^2)').value*masses/((radii)**2))*u.Rsun/u.s**2).to('cm/s^2')).value)

    def get_kick_differential(self, delta_v_sys_x, delta_v_sys_y, delta_v_sys_z, phase=None, inclination=None):
        """Calculate the :class:`~astropy.coordinates.CylindricalDifferential` from a combination of the natal
        kick, Blauuw kick and orbital motion.
    
        via cogsworth https://github.com/TomWagg/cogsworth/blob/main/cogsworth/kicks.py
    
        Parameters
        ----------
        delta_v_sys_x : :class:`~astropy.units.Quantity` [velocity]
            Change in systemic velocity due to natal and Blauuw kicks in BSE :math:`(v_x, v_y, v_z)` frame
            (see Fig A1 of `Hurley+02 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_)
        delta_v_sys_y : :class:`~astropy.units.Quantity` [velocity]
            Change in systemic velocity due to natal and Blauuw kicks in BSE :math:`(v_x, v_y, v_z)` frame
            (see Fig A1 of `Hurley+02 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_)
        delta_v_sys_z : :class:`~astropy.units.Quantity` [velocity]
            Change in systemic velocity due to natal and Blauuw kicks in BSE :math:`(v_x, v_y, v_z)` frame
            (see Fig A1 of `Hurley+02 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_)
        phase : np.array
            Orbital phase angle in radians
        inclination : np.array
            Inclination to the Galactic plane in radians
    
        Returns
        -------
        kick_differential : :class:`~astropy.coordinates.CylindricalDifferential`
            Kick differential
        """
        # orbital phase angle and inclination to Galactic plane
        thetas = np.random.uniform(0, 2 * np.pi, size = len(delta_v_sys_x)) if phase is None else phase
        phis = np.arccos(2 * np.random.rand(len(delta_v_sys_x)) - 1.0) if inclination is None else inclination
    
        # rotate BSE (v_x, v_y, v_z) into Galactocentric (v_X, v_Y, v_Z)
        v_X = delta_v_sys_x * np.cos(thetas) - delta_v_sys_y * np.sin(thetas) * np.cos(phis)\
            + delta_v_sys_z * np.sin(thetas) * np.sin(phis)
        v_Y = delta_v_sys_x * np.sin(thetas) + delta_v_sys_y * np.cos(thetas) * np.cos(phis)\
            - delta_v_sys_z * np.cos(thetas) * np.sin(phis)
        v_Z = delta_v_sys_y * np.sin(phis) + delta_v_sys_z * np.cos(phis)
        kick_differential = coords.CartesianDifferential(v_X, v_Y, v_Z)
    
        return kick_differential


#==============================#
# Merged model classes
#==============================#
class MergedFlameSkimmerRedDiamondbackMIST(StellarEvolution):
    """
    Load the merged Sonora Flame Skimmer, Sonora Red Diamondback, and
    MIST v1.2 evolution grid.

    Flame Skimmer supplies the low-mass evolution through 9 Gyr. At older
    ages Red Diamondback is used only at its native metallicities (-0.5, 0.0,
    and +0.5); the other merged-grid metallicities retain explicitly flagged
    Flame age extrapolations. MIST v1.2 supplies the stellar sequence.

    The merged grid covers [M/H] = -1.0 to +0.5 in 0.25-dex steps. For a
    request below -1.0, the best available [M/H]=-1.0 brown-dwarf evolution
    rows are reused as an explicitly warned proxy and joined with MIST v1.2 at
    the requested metallicity. Photometry can remain unavailable for old or
    low-metallicity brown dwarfs when the atmosphere grid has no matching
    spectrum; the evolution loader does not invent fluxes.

    Parameters
    ----------
    model_dir : str or path-like, optional
        Root directory containing the seven metallicity subdirectories from
        ``iso/z0014`` through ``iso/z045``. When omitted, use the installed SPISEA model location
        ``evolution/merged/flame_red_diamondback_mist``.
    """

    def __init__(self, model_dir=None):
        self.z_solar = 0.0142
        self.metallicity_list = np.array(
            [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5], dtype=float
        )
        z_list = list(self.z_solar * 10.0**self.metallicity_list)

        regular_ages = np.round(np.arange(5.0, 10.15 + 0.025, 0.05), 2)
        exact_ages = np.array([
            np.log10(9.0e9),
            np.log10(15.0e9),
        ])
        age_list = np.array(sorted(set(regular_ages) | set(exact_ages)))

        # The mass sampling varies with age, so there is no single exact list.
        mass_list = []

        if model_dir is None:
            model_dir = os.path.join(
                models_dir,
                'merged',
                'flame_red_diamondback_mist',
            )
            self._uses_installed_grid = True
        else:
            model_dir = os.path.abspath(os.fspath(model_dir))
            self._uses_installed_grid = False

        StellarEvolution.__init__(
            self,
            model_dir,
            age_list,
            mass_list,
            z_list,
        )
        self.model_version_name = (
            'MergedFlameSkimmerRedDiamondbackMIST'
        )
        self.z_solar = 0.0142
        self.metallicity_list = np.array(
            [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5], dtype=float
        )
        self.z_file_map = {
            -1.0: 'z0014',
            -0.75: 'z0025',
            -0.5: 'z0045',
            -0.25: 'z0080',
            0.0: 'z014',
            0.25: 'z025',
            0.5: 'z045',
        }
        self._mist_only_evolution = None
        self._mist_only_min_metallicity = -4.0
        self._low_metallicity_warning_emitted = False

        # Minimum version of the currently installed SPISEA evolution grid.
        self.evo_grid_min = 3.0

    @staticmethod
    def _isochrone_filename(log_age):
        if np.isclose(log_age, np.log10(9.0e9), rtol=0.0, atol=1.0e-8):
            return 'iso_9.954243.fits'
        if np.isclose(log_age, np.log10(15.0e9), rtol=0.0, atol=1.0e-8):
            return 'iso_10.176091.fits'
        return 'iso_{0:.2f}.fits'.format(log_age)

    def _read_static_merged_isochrone(self, log_age, metallicity):
        iso_file = self._isochrone_filename(log_age)
        full_iso_file = os.path.join(
            self.model_dir,
            'iso',
            self.z_file_map[float(metallicity)],
            iso_file,
        )
        if not os.path.isfile(full_iso_file):
            raise FileNotFoundError(
                'Isochrone file not found: {0}'.format(full_iso_file)
            )
        return Table.read(full_iso_file, format='fits')

    @staticmethod
    def _standardize_proxy_component(frame, model_ref):
        output = frame.copy()
        required_input = {'mass', 'mass_current', 'logL', 'logT'}
        missing_input = sorted(required_input - set(output.columns))
        if missing_input:
            raise ValueError(
                'Low-metallicity proxy component is missing columns: {0}'
                .format(', '.join(missing_input))
            )
        # SPISEA's ordinary MIST loader does not always expose logR. Derive
        # it from the Stefan-Boltzmann identity rather than rejecting an
        # otherwise valid requested-metallicity MIST isochrone. Likewise,
        # derive logg only when the delegated loader did not provide it.
        if 'logR' not in output:
            output['logR'] = (
                0.5 * pd.to_numeric(output['logL'], errors='raise')
                - 2.0 * (
                    pd.to_numeric(output['logT'], errors='raise')
                    - np.log10(5772.0)
                )
            )
        if 'logg' not in output:
            output['logg'] = (
                4.438
                + np.log10(pd.to_numeric(output['mass_current'], errors='raise'))
                - 2.0 * pd.to_numeric(output['logR'], errors='raise')
            )
        required = {'mass', 'mass_current', 'logL', 'logT', 'logR', 'logg'}
        missing = sorted(required - set(output.columns))
        if missing:
            raise ValueError(
                'Low-metallicity proxy component is missing columns: {0}'
                .format(', '.join(missing))
            )
        if 'logT_WR' not in output:
            output['logT_WR'] = output['logT']
        if 'phase' not in output:
            output['phase'] = 0
        if 'isWR' not in output:
            output['isWR'] = False
        if 'model_ref' not in output:
            output['model_ref'] = model_ref
        if 'is_interpolated' not in output:
            output['is_interpolated'] = False
        if 'is_extrapolated' not in output:
            output['is_extrapolated'] = False
        if 'interp_axis' not in output:
            output['interp_axis'] = 'none'
        if 'extrap_axis' not in output:
            output['extrap_axis'] = 'none'
        return output

    def _low_metallicity_proxy_isochrone(
            self, requested_age, actual_log_age, requested_metallicity):
        """Combine -1 brown dwarfs with MIST at the requested [M/H]."""
        if not self._low_metallicity_warning_emitted:
            warnings.warn(
                'Requested [M/H]={0:+.2f} is below the Flame Skimmer '
                'limit. Brown-dwarf evolution is approximated with the '
                '[M/H]=-1.00 merged grid, while MIST v1.2 uses the requested '
                'metallicity. Brown-dwarf photometry may be unavailable '
                'where the atmosphere grid has no matching spectrum.'
                .format(requested_metallicity),
                UserWarning,
                stacklevel=3,
            )
            self._low_metallicity_warning_emitted = True

        minus_one = self._read_static_merged_isochrone(
            actual_log_age, -1.0
        ).to_pandas()
        minus_one = self._standardize_proxy_component(
            minus_one, 'Sonora [M/H]=-1 brown-dwarf proxy'
        )
        brown_dwarfs = minus_one.loc[
            pd.to_numeric(minus_one['mass'], errors='coerce') <= 0.080000001
        ].copy()
        if len(brown_dwarfs) < 4:
            raise ValueError(
                'The [M/H]=-1 merged isochrone has too few brown-dwarf rows '
                'for the low-metallicity proxy.'
            )
        brown_dwarfs['model_ref'] = 'Sonora [M/H]=-1 brown-dwarf proxy'

        if self._mist_only_evolution is None:
            self._mist_only_evolution = MISTv1(
                version=1.2, synthpop_extension=False
            )
        delegated_age = float(10.0**actual_log_age)
        mist_table = self._mist_only_evolution.isochrone(
            age=delegated_age,
            metallicity=requested_metallicity,
        )
        mist = self._standardize_proxy_component(
            mist_table.to_pandas(), 'original MIST v1.2'
        )
        mist = mist.loc[
            pd.to_numeric(mist['mass'], errors='coerce') > 0.080000001
        ].sort_values('mass', kind='stable').reset_index(drop=True)
        brown_dwarfs = brown_dwarfs.sort_values(
            'mass', kind='stable'
        ).reset_index(drop=True)
        if len(mist) < 4:
            raise ValueError(
                'MIST v1.2 has too few rows above 0.08 Msun for the '
                'low-metallicity proxy junction.'
            )

        low_edge = float(brown_dwarfs['mass'].iloc[-1])
        high_edge = float(mist['mass'].iloc[0])
        if high_edge <= low_edge:
            raise ValueError(
                'The low-metallicity brown-dwarf and MIST components overlap '
                'or are not ordered in initial mass.'
            )
        step = 0.002
        first = math.ceil((low_edge + 1.0e-12) / step) * step
        bridge_mass = np.round(
            np.arange(first, high_edge - 1.0e-12, step), 12
        )
        if len(bridge_mass) == 0:
            raise ValueError(
                'No PCHIP bridge rows fit inside the low-metallicity mass gap.'
            )

        low_anchor = brown_dwarfs.tail(4)
        high_anchor = mist.head(4)
        anchor_mass = np.concatenate([
            low_anchor['mass'].to_numpy(dtype=float),
            high_anchor['mass'].to_numpy(dtype=float),
        ])
        if not np.all(np.diff(anchor_mass) > 0):
            raise ValueError(
                'Low-metallicity PCHIP anchors are not strictly increasing.'
            )
        bridge = pd.DataFrame({'mass': bridge_mass})
        for quantity in ('logL', 'logT'):
            anchor_values = np.concatenate([
                low_anchor[quantity].to_numpy(dtype=float),
                high_anchor[quantity].to_numpy(dtype=float),
            ])
            function = interpolate.PchipInterpolator(
                anchor_mass, anchor_values, extrapolate=False
            )
            bridge[quantity] = function(bridge_mass)
            lower = min(float(anchor_values[3]), float(anchor_values[4]))
            upper = max(float(anchor_values[3]), float(anchor_values[4]))
            values = bridge[quantity].to_numpy(dtype=float)
            if values.min() < lower - 1.0e-10 or values.max() > upper + 1.0e-10:
                raise ValueError(
                    'Low-metallicity PCHIP {0} left its endpoint range.'
                    .format(quantity)
                )

        solar_teff = 5772.0
        solar_logg = 4.438
        bridge['mass_current'] = bridge['mass']
        bridge['logR'] = (
            0.5 * bridge['logL']
            - 2.0 * (bridge['logT'] - np.log10(solar_teff))
        )
        bridge['logg'] = (
            solar_logg + np.log10(bridge['mass_current'])
            - 2.0 * bridge['logR']
        )
        bridge['logT_WR'] = bridge['logT']
        bridge['phase'] = 99
        bridge['isWR'] = False
        bridge['model_ref'] = 'PCHIP [M/H]=-1 proxy-to-MIST bridge'
        bridge['is_interpolated'] = True
        bridge['is_extrapolated'] = bool(
            brown_dwarfs['is_extrapolated'].astype(bool).any()
        )
        bridge['interp_axis'] = 'mass'
        bridge['extrap_axis'] = np.where(
            bridge['is_extrapolated'], 'age', 'none'
        )

        columns = [
            'mass', 'mass_current', 'logL', 'logT', 'logR', 'logg',
            'logT_WR', 'phase', 'isWR', 'model_ref', 'is_interpolated',
            'is_extrapolated', 'interp_axis', 'extrap_axis',
        ]
        combined = pd.concat(
            [
                brown_dwarfs.reindex(columns=columns),
                bridge.reindex(columns=columns),
                mist.reindex(columns=columns),
            ],
            ignore_index=True,
        ).sort_values('mass', kind='stable').reset_index(drop=True)
        numeric = combined[
            ['mass', 'mass_current', 'logL', 'logT', 'logR', 'logg']
        ].to_numpy(dtype=float)
        if not np.isfinite(numeric).all():
            raise ValueError(
                'Low-metallicity proxy isochrone contains non-finite '
                'evolution values.'
            )
        if not np.all(np.diff(combined['mass'].to_numpy(dtype=float)) > 0):
            raise ValueError(
                'Low-metallicity proxy isochrone masses are not strictly '
                'increasing.'
            )

        combined['phase'] = combined['phase'].astype(int)
        for name in ('isWR', 'is_interpolated', 'is_extrapolated'):
            combined[name] = combined[name].fillna(False).astype(bool)
        for name in ('model_ref', 'interp_axis', 'extrap_axis'):
            combined[name] = combined[name].fillna('none').astype(str)

        iso = Table.from_pandas(combined)
        iso.meta['log_age'] = math.log10(requested_age)
        iso.meta['log_age_act'] = actual_log_age
        iso.meta['metallicity_in'] = requested_metallicity
        iso.meta['metallicity_act'] = requested_metallicity
        iso.meta['brown_dwarf_source_metallicity'] = -1.0
        iso.meta['metallicity_policy'] = (
            '[M/H]=-1 brown-dwarf evolution proxy plus MIST v1.2 at the '
            'requested metallicity'
        )
        iso.meta['brown_dwarf_photometry'] = (
            'may be unavailable where atmosphere coverage is missing'
        )
        return iso

    def isochrone(self, age=1.e8, metallicity=0.0):
        """Return the nearest available isochrone for age in years."""
        if age <= 0:
            raise ValueError('Age must be positive and in years.')

        if self._uses_installed_grid:
            self.evo_grid_num = check_evo_grid_number(
                self.evo_grid_min,
                models_dir,
            )

        log_age = math.log10(age)
        age_array = np.asarray(self.age_list, dtype=float)
        if log_age < age_array.min() or log_age > age_array.max():
            raise ValueError(
                'Requested logAge {0} is out of bounds between {1} and {2}.'
                .format(log_age, age_array.min(), age_array.max())
            )

        age_index = np.argmin(np.abs(age_array - log_age))
        actual_log_age = float(age_array[age_index])

        if metallicity < self.metallicity_list.min():
            if metallicity < self._mist_only_min_metallicity:
                raise ValueError(
                    'Requested metallicity {0} is below the MIST v1.2 '
                    'minimum of {1}.'.format(
                        metallicity, self._mist_only_min_metallicity
                    )
                )
            return self._low_metallicity_proxy_isochrone(
                requested_age=float(age),
                actual_log_age=actual_log_age,
                requested_metallicity=float(metallicity),
            )

        if metallicity > self.metallicity_list.max():
            raise ValueError(
                'Requested metallicity {0} is above the merged-grid maximum '
                'of {1}.'.format(
                    metallicity,
                    self.metallicity_list.max(),
                )
            )

        metallicity_index = np.argmin(
            np.abs(self.metallicity_list - metallicity)
        )
        actual_metallicity = float(
            self.metallicity_list[metallicity_index]
        )

        iso = self._read_static_merged_isochrone(
            actual_log_age, actual_metallicity
        )
        required_columns = {
            'mass', 'mass_current', 'logL', 'logT', 'logR', 'logg',
            'logT_WR', 'phase', 'isWR',
        }
        missing = sorted(required_columns - set(iso.colnames))
        if missing:
            raise ValueError(
                'Merged isochrone is missing columns: {0}'.format(
                    ', '.join(missing)
                )
            )

        iso.meta['log_age'] = log_age
        iso.meta['log_age_act'] = actual_log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = actual_metallicity
        return iso


class MergedPhillipsBaraffePisaEkstromParsec(StellarEvolution):
    """
    This is a combination of several different evolution models:

    * Phillips (`Phillips et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020A%26A...637A..38P/abstract>`_)
    * Baraffe (`Baraffe et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract>`_)
    * Pisa (`Tognelli et al. 2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...533A.109T/abstract>`_)
    * Geneva (`Ekstrom et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.146E/abstract>`_)
    * Parsec (version 1.2s, `Bressan+12 <https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract>`_)

    The model used depends on the age of the population and what stellar masses
    are being modeled:
    

    For logAge < 7.4:

    * Phillips: 0.01 - 0.07 M_sun
    * Phillips/Baraffe transition: 0.070 - 0.075 M_sun
    * Baraffe: 0.075 - 0.4 M_sun
    * Baraffe/Pisa transition: 0.4 - 0.5 M_sun 
    * Pisa: 0.5 M_sun to the highest mass in Pisa isochrone (typically 5 - 7 Msun)
    * Geneva: Highest mass of Pisa models to 120 M_sun

    For logAge > 7.4:

    * Phillips: 0.01 - 0.075 M_sun
    * Phillips/Parsec v1.2s transition: 0.075 - 0.2 M_sun
    * Parsec v1.2s: full mass range above 0.2 M_sun
    
    Parameters
    ----------
    rot: boolean, optional
        If true, then use rotating Ekstrom models. Default is true.
    """
    def __init__(self, rot=True):
        # populate list of model masses (in solar masses)
        mass_list = [(0.01 + i*0.005) for i in range(181)] # generates masses from 0.01 - 1 M_sun
        
        # define metallicity parameters for Geneva models
        z_list = [0.015]
        
        # populate list of isochrone ages (log scale)
        age_list = np.arange(6.0, 10.0, 0.01).tolist()
        
        # specify location of model files
        model_dir = models_dir + 'merged/phillips_baraffe_pisa_ekstrom_parsec/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.015
        
        # Switch to specify rotating/non-rotating models
        if rot:
            self.z_file_map = {0.015: 'z015_rot/'}
        else:
            self.z_file_map = {0.015: 'z015_norot/'}

        # Define required evo_grid number
        self.evo_grid_min = 3.0
        
    
    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Baraffe-Pisa-Ekstrom-Parsec 
        collection
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)
        
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)
        
        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            logger.error('Requested age {0} is out of bounds.'.format(log_age))
            
        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            logger.error('Requested metallicity {0} is out of bounds.'.format(z_defined))

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.dat'.format(self.age_list[age_idx])
        
        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + z_dir + iso_file

        # return isochrone data
        iso = Table.read(full_iso_file, format='ascii')
        iso.rename_column('col1', 'mass')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'logL')
        iso.rename_column('col4', 'logg')
        iso.rename_column('col5', 'logT_WR')
        iso.rename_column('col6', 'mass_current')
        iso.rename_column('col7', 'phase')
        iso.rename_column('col8', 'model_ref')

        # Define "isWR" column based on phase info
        isWR = Column([False] * len(iso), name='isWR')
        idx_WR = np.where(iso['logT'] != iso['logT_WR'])
        isWR[idx_WR] = True
        iso.add_column(isWR)

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)
        
        # Assume mass of brown dwarfs does not change over their lifetime
        #bd_idx = iso['mass'] < 0.08
        #iso['mass_current'][bd_idx] = iso['mass'][bd_idx]

        # Handling NaN effective temperatures
        #nan_teff_idx = np.isnan(iso['logT'])
        #if np.any(nan_teff_idx):
        #    iso['logT'][nan_teff_idx] = self.estimate_teff(iso['mass'][nan_teff_idx])
            
        return iso


class MergedBaraffePisaEkstromParsec(StellarEvolution):
    """
    This is a combination of several different evolution models:

    * Baraffe (`Baraffe et al. 2015 <https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/abstract>`_)
    * Pisa (`Tognelli et al. 2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...533A.109T/abstract>`_)
    * Geneva (`Ekstrom et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.146E/abstract>`_)
    * Parsec (version 1.2s, `Bressan+12 <https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract>`_)

    The model used depends on the age of the population and what stellar masses
    are being modeled:


    For logAge < 7.4:

    * Baraffe: 0.08 - 0.4 M_sun
    * Baraffe/Pisa transition: 0.4 - 0.5 M_sun
    * Pisa: 0.5 M_sun to the highest mass in Pisa isochrone (typically 5 - 7 Msun)
    * Geneva: Highest mass of Pisa models to 120 M_sun

    For logAge > 7.4:

    * Parsec v1.2s: full mass range

    Parameters
    ----------
    rot: boolean, optional
        If true, then use rotating Ekstrom models. Default is true.
    """
    def __init__(self, rot=True):
        if rot:
            self.model_version_name = "MergedBaraffePisaEkstromParsec-rot"
        else:
            self.model_version_name = "MergedBaraffePisaEkstromParsec-norot"
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]

        # define metallicity parameters for Geneva models
        z_list = [0.015]

        # populate list of isochrone ages (log scale)
        age_list = np.arange(6.0, 10.091, 0.01).tolist()

        # specify location of model files
        model_dir = models_dir + 'merged/baraffe_pisa_ekstrom_parsec/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.015

        # Switch to specify rotating/non-rotating models
        if rot:
            self.z_file_map = {0.015: 'z015_rot/'}
        else:
            self.z_file_map = {0.015: 'z015_norot/'}

        # Define required evo_grid number
        self.evo_grid_min = 1.0


    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Baraffe-Pisa-Ekstrom-Parsec
        collection
        """
        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)

        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if ((log_age < np.min(self.age_list)) or (log_age > np.max(self.age_list))):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if ((z_defined < np.min(self.z_list)) or
                (z_defined > np.max(self.z_list))):
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])


        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + z_dir + iso_file

        # return isochrone data
        try:
            # FITS version of file (older model evolution grids)
            iso = Table.read(full_iso_file, format='fits')
        except:
            # ASCII version of files (newer model evo grids
            iso_file = 'iso_{0:.2f}.dat'.format(self.age_list[age_idx])
            full_iso_file = self.model_dir + z_dir + iso_file

            iso = Table.read(full_iso_file, format='ascii')

        iso.rename_column('col1', 'mass')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'logL')
        iso.rename_column('col4', 'logg')
        iso.rename_column('col5', 'logT_WR')
        iso.rename_column('col6', 'mass_current')
        iso.rename_column('col7', 'phase')
        iso.rename_column('col8', 'model_ref')

        # Define "isWR" column based on phase info
        isWR = Column([False] * len(iso), name='isWR')
        idx_WR = np.where(iso['logT'] != iso['logT_WR'])
        isWR[idx_WR] = True
        iso.add_column(isWR)

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso


class MergedPisaEkstromParsec(StellarEvolution):
    """
    Same as MergedBaraffePisaEkstromParsec, but without
    the Baraffe models.

    Parameters
    ----------
    rot: boolean, optional
        If true, then use rotating Ekstrom models. Default is true.
    """
    def __init__(self, rot=True):
        if rot:
            self.model_version_name = "MergedPisaEkstromParsec-rot"
        else:
            self.model_version_name = "MergedPisaEkstromParsec-norot"
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]

        # define metallicity parameters for Geneva models
        z_list = [0.015]

        # populate list of isochrone ages (log scale)
        age_list = np.arange(6.0, 8.001, 0.01).tolist()

        # specify location of model files
        model_dir = models_dir + 'merged/pisa_ekstrom_parsec/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.015

        #Switch to specify rot/notot
        if rot:
            self.z_file_map = {0.015: 'z015_rot/'}
        else:
            self.z_file_map = {0.015: 'z015_norot/'}

        # Define required evo_grid number
        self.evo_grid_min = 1.0

        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)


    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Pisa-Ekstrom-Parsec collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if (log_age < self.age_list[0]) or (log_age > self.age_list[-1]):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if not z_defined in self.z_list:
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + z_dir + iso_file

        # return isochrone data
        iso = Table.read(full_iso_file, format='fits')
        iso.rename_column('col1', 'mass')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'logL')
        iso.rename_column('col4', 'logg')
        iso.rename_column('col5', 'logT_WR')
        iso.rename_column('col6', 'model_ref')

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

class MergedSiessGenevaPadova(StellarEvolution):
    """
    This is a combination of several different evolution models.

    The model used depends on the age of the population and what stellar masses
    are being modeled:

    * Siess (`Siess et al. 2000 <https://ui.adsabs.harvard.edu/abs/2000A%26A...358..593S/abstractt>`_)
    * Geneva (`Meynet & Maeder 2003 <https://ui.adsabs.harvard.edu/abs/2003A%26A...404..975M/abstract>`_)
    * Padova (`Marigo et al. 2008 <https://ui.adsabs.harvard.edu/abs/2008A%26A...482..883M/abstract>`_)

    For logAge < 7.4:

    * Siess: 0.1 - 7 M_sun
    * Siess/Geneva transition: 7 - 9 M_sun
    * Geneva: > 9 M_sun

    For logAge > 7.4:

    * Padova: full mass range
    """
    def __init__(self):
        """
        Define intrinsic properties for merged Siess-meynetMaeder-Padova
        stellar models.
        """
        self.model_version_name = "MergedSiessGenevaPadova"
        # populate list of model masses (in solar masses)
        mass_list = [(0.1 + i*0.005) for i in range(181)]

        # define metallicity parameters for Geneva models
        z_list = [0.02]

        # populate list of isochrone ages (log scale)
        age_list = np.arange(5.5, 7.41, 0.01).tolist()
        age_list.append(7.48)
        idx = np.arange(7.50, 8.01, 0.05)
        for ii in idx:
            age_list.append(ii)
        age_list.append(8.30)
        age_list.append(8.48)
        age_list.append(8.60)
        age_list.append(8.70)
        age_list.append(8.78)
        age_list.append(8.85)
        age_list.append(8.90)
        age_list.append(8.95)
        age_list.append(9.00)
        age_list.append(9.30)
        age_list.append(9.60)
        age_list.append(9.70)
        age_list.append(9.78)

        # specify location of model files
        model_dir = models_dir + 'merged/siess_meynetMaeder_padova/'
        StellarEvolution.__init__(self, model_dir, age_list, mass_list, z_list)
        self.z_solar = 0.02

        # Metallicity map
        self.z_file_map = {0.02: 'z02/'}

        # Define required evo_grid number
        self.evo_grid_min = 1.0

        # Error check to see if installed evolution model
        # grid is compatible with code version. Also return
        # current grid num
        self.evo_grid_num = check_evo_grid_number(self.evo_grid_min, models_dir)


    def isochrone(self, age=1.e8, metallicity=0.0):
        r"""
        Extract an individual isochrone from the Siess-Geneva-Padova collection.
        """
        # convert metallicity to mass fraction
        z_defined = self.z_solar*10.**metallicity

        log_age = math.log10(age)

        # check age and metallicity are within bounds
        if (log_age < self.age_list[0]) or (log_age > self.age_list[-1]):
            raise ValueError(f'Requested age {log_age} is out of bounds between {np.min(self.age_list)} and {np.max(self.age_list)}.')

        if not z_defined in self.z_list:
            raise ValueError(f'Requested metallicity {z_defined} is out of bounds between {np.min(self.z_list)} and {np.max(self.z_list)}.')

        # Find nearest age in grid to input grid
        age_idx = np.where(abs(np.array(self.age_list) - log_age) == min(abs(np.array(self.age_list) - log_age)) )[0][0]
        iso_file = 'iso_{0:.2f}.fits'.format(self.age_list[age_idx])

        # find closest metallicity value
        z_idx = np.where(abs(np.array(self.z_list) - z_defined) == min(abs(np.array(self.z_list) - z_defined)) )[0][0]
        z_dir = self.z_file_map[self.z_list[z_idx]]

        # generate isochrone file string
        full_iso_file = self.model_dir + z_dir + iso_file

        # return isochrone data
        iso = Table.read(full_iso_file, format='ascii')
        iso.rename_column('col1', 'mass')
        iso.rename_column('col2', 'logT')
        iso.rename_column('col3', 'logL')
        iso.rename_column('col4', 'logg')
        iso.rename_column('col5', 'logT_WR')
        iso.rename_column('col6', 'model_ref')

        iso.meta['log_age'] = log_age
        iso.meta['metallicity_in'] = metallicity
        iso.meta['metallicity_act'] = np.log10(self.z_list[z_idx] / self.z_solar)

        return iso

#================================================#

def make_isochrone_pisa_interp(log_age, metallicity=0.015,
                         tracks=None, test=False):
    """
    Read in a set of isochrones and generate an isochrone at log_age
    that is well sampled at the full range of masses.

    Puts isochrones is Pisa2011/iso/<metal>/
    """
    # If logage > 8.0, quit immediately...grid doesn't go that high
    if log_age > 8.0:
        print( 'Age too high for Pisa grid (max logAge = 8.0)')
        return

    # Directory with where the isochrones will go (both downloaded and interpolated)
    rootDir = models_dir + '/Pisa2011/iso/'
    metSuffix = 'z' + str(metallicity).split('.')[-1]
    rootDir += metSuffix + '/'

    # Can we find the isochrone directory?
    if not os.path.exists(rootDir):
        print( 'Failed to find Pisa PMS isochrones for metallicity = ' + metSuffix)
        return

    # Check to see if isochrone at given age already exists. If so, quit
    if os.path.exists(rootDir+'iso_{0:3.2f}.fits'.format(log_age)):
        print( 'Isochrone at logAge = {0:3.2f} already exists'.format(log_age))
        return

    # Name/directory for interpolated isochrone
    isoFile = rootDir+'iso_%3.2f.fits' % log_age
    outSuffix = '_%.2f' % (log_age)

    print( '*** Generating Pisa isochrone for log t = %3.2f and Z = %.3f' % \
        (log_age, metallicity))
    import time
    print( time.asctime(), 'Getting original Pisa isochrones.')
    iso = get_orig_pisa_isochrones(metallicity=metallicity)

    # First thing is to find the isochrones immediately above and below desired
    # age
    iso_log_ages = iso.log_ages
    tmp = np.append(iso_log_ages, log_age)

    # Find desired age in ordered sequence; isolate model younger and older
    tmp.sort()
    good = np.where(tmp == log_age)
    young_model_logage = tmp[good[0]-1]
    old_model_logage = tmp[good[0]+1]

    # Isolate younger/older isochrones
    young_ind = np.where(iso.log_ages == young_model_logage)
    old_ind = np.where(iso.log_ages == old_model_logage)

    young_iso = iso.isochrones[young_ind[0]]
    old_iso = iso.isochrones[old_ind[0]]

    # Need both younger and older model on same temperature grid for time
    # interpolation. Will adopt mass grid of whichever model is closer in time
    if abs(young_model_logage - log_age) <= abs(old_model_logage - log_age):
        # Use young model mass grid
        young_iso, old_iso = interpolate_iso_tempgrid(young_iso, old_iso)

    else:
        # Use old model mass grid
        old_iso, young_iso = interpolate_iso_tempgrid(old_iso, young_iso)

    # Now, can interpolate in time over the two models. Do this star by star.
    # Work in linear time here!!
    numStars = len(young_iso.M)

    interp_iso = Isochrone(log_age)
    interp_iso.log_Teff = np.zeros(numStars, dtype=float)
    interp_iso.log_L = np.zeros(numStars, dtype=float)
    interp_iso.log_g = np.zeros(numStars, dtype=float)
    interp_iso.M = young_iso.M # Since mass grids should already be matched

    for i in range(numStars):
        # Do interpolations in linear space
        model_ages = [10**young_model_logage[0], 10**old_model_logage[0]]
        target_age = 10**log_age
        #model_ages = [young_model_logage[0], old_model_logage[0]]
        #target_age = log_age

        # Build interpolation functions
        Teff_arr = [10**young_iso.log_Teff[i], 10**old_iso.log_Teff[i]]
        logL_arr = [10**young_iso.log_L[i], 10**old_iso.log_L[i]]
        logg_arr = [10**young_iso.log_g[i], 10**old_iso.log_g[i]]

        f_log_Teff = interpolate.interp1d(model_ages, Teff_arr, kind='linear')
        f_log_L = interpolate.interp1d(model_ages, logL_arr, kind='linear')
        f_log_g = interpolate.interp1d(model_ages, logg_arr, kind='linear')

        interp_iso.log_Teff[i] = np.log10(f_log_Teff(target_age))
        interp_iso.log_L[i] = np.log10(f_log_L(target_age))
        interp_iso.log_g[i] = np.log10(f_log_g(target_age))

    # If indicated, plot new isochrone along with originals it was interpolated
    # from
    if test:
        py.figure(1)
        py.clf()
        py.plot(interp_iso.log_Teff, interp_iso.log_L, 'k-', label = 'Interp')
        py.plot(young_iso.log_Teff, young_iso.log_L, 'b-',
                label = 'log Age = {0:3.2f}'.format(young_model_logage[0]))
        py.plot(old_iso.log_Teff, old_iso.log_L, 'r-',
                label = 'log Age = {0:3.2f}'.format(old_model_logage[0]))
        rng = py.axis()
        py.xlim(rng[1], rng[0])
        py.xlabel('log Teff')
        py.ylabel('log L')
        py.legend()
        py.title('Pisa 2011 Isochrone at log t = %.2f' % log_age)
        py.savefig(rootDir + 'plots/interp_isochrone_at' + outSuffix + '.png')

    print( time.asctime(), 'Finished.')

    # Write output to file, MUST BE IN SAME ORDER AS ORIG FILES
    _out = open(isoFile, 'w')

    _out.write('%10s  %10s  %10s  %10s\n' %
               ('# log L', 'log Teff', 'Mass', 'log g'))
    _out.write('%10s  %10s  %10s  %10s\n' %
               ('# (Lsun)', '(Kelvin)', '(Msun)', '(cgs)'))

    for ii in range(len(interp_iso.M)):
        _out.write('%10.4f  %10.4f  %10.4f  %10.4f\n' %
                   (interp_iso.log_L[ii], interp_iso.log_Teff[ii], interp_iso.M[ii],
                    interp_iso.log_g[ii]))

    _out.close()

    return

def get_orig_pisa_isochrones(metallicity=0.015):
    """
    Helper code to get the original pisa isochrones at given metallicity.
    These are downloaded online
    """
    pms_dir = models_dir + '/Pisa2011/iso/iso_orig/'
    metSuffix = 'z' + str(metallicity).split('.')[-1]
    pms_dir += metSuffix + '/'

    if not os.path.exists(pms_dir):
        print( 'Failed to find Siess PMS isochrones for metallicity = ' + metSuffix)
        return

    # Collect the isochrones
    files = glob.glob(pms_dir + '*.dat')
    count = len(files)

    data = objects.DataHolder()

    data.isochrones = []
    data.log_ages = []

    # Extract useful params from isochrones
    for ff in range(len(files)):
        d = Table.read(files[ff], format='ascii')

        # Extract logAge from filename
        log_age = float(files[ff].split('_')[2][:-4])

        # Create an isochrone object
        iso = Isochrone(log_age)
        iso.M = d['col3']
        iso.log_Teff = d['col2']
        iso.log_L = d['col1']

        # If a log g column exist, extract it. Otherwise, calculate
        # log g from T and L and add column at end
        if len(d.keys()) == 3:

            # Calculate log g from T and L
            L_sun = 3.8 * 10**33 #cgs
            SB_sig = 5.67 * 10**-5 #cgs
            M_sun = 2. * 10**33 #cgs
            G_const = 6.67 * 10**-8 #cgs

            radius = np.sqrt( (10**d['col1'] * L_sun) /
                          (4 * np.pi * SB_sig *  (10**d['col2'])**4) )
            g = (G_const * d['col3'] * M_sun) / radius**2


            iso.log_g = np.log10(g.astype(np.float))
        else:
            iso.log_g = d['col4']

        data.isochrones.append(iso)
        data.log_ages.append(log_age)

        # If it doesn't already exist, add a column with logg vals. This will
        # be appended at the end
        if len(d.keys()) == 3:
            logg_col = Column(iso.log_g, name = 'col4')
            d.add_column(logg_col, index=3)
            d.write(files[ff],format='ascii')
    data.log_ages = np.array(data.log_ages)

    # Resort so that everything is in order of increasing age
    sdx = data.log_ages.argsort()
    data.masses = data.log_ages[sdx]
    data.isochrones = [data.isochrones[ss] for ss in sdx]

    return data

class Isochrone(object):
    def __init__(self, log_age):
        self.log_age = log_age
