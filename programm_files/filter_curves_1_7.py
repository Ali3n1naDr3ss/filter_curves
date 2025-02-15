# loops

import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.gridspec import GridSpec

nipy_spectral = plt.get_cmap('nipy_spectral')
PuRd = plt.get_cmap('PuRd')
YlOrRd = plt.get_cmap('YlGn')
oranges = plt.get_cmap('Oranges')

LBG_spectrum = '/Users/user/Documents/filter_curves/data_files/REBELS-05_1dspec_wRMS.dat'
LBG_wavelength = []
LBG_flux = []
with open(LBG_spectrum,'r', encoding='utf-8-sig') as file:
    lines = file.readlines()[2:]
    for line in lines:
        columns = line.split()
        LBG_wavelength.append(float(columns[0])) # microns
        LBG_flux.append(float(columns[1])) # nJy

z = 6.496
c = 3e18 # Angstrom/s
def LBG_spectrum_prep(LBG_wavelength, LBG_flux):

    def flux_per_wavelength(flux, wavelength):
        '''
        Convets nJy to erg s^-1 cm^-2 Å^-1
        '''
        flux_per_wavelength = []
        for flux, wavelength in zip(flux, wavelength):
            flux_per_A = (flux * c / (wavelength**2) * 1e-21)
            flux_per_wavelength.append(flux_per_A)
        return flux_per_wavelength

    LBG_flux = flux_per_wavelength(LBG_flux, LBG_wavelength)
    LBG_flux_scaled = []
    for i in range(len(LBG_flux)):
        LBG_flux[i] = LBG_flux[i] * 40 # scaling factor
        LBG_flux_scaled.append(LBG_flux[i])

    return LBG_flux_scaled
LBG_flux = LBG_spectrum_prep(LBG_wavelength, LBG_flux)

def Lyman_break_position(z):
    '''
    INPUT(s)
    z(float): redshift of LBG spectrum

    OUTPUT(s)
    lyman_break_observed(float): observed wavelength of Lyman break in microns
    y(np.array): y-axis values for the Lyman break plot
    '''
    y = np.linspace(-100, 1500, 1000)
    lyman_break_rest = 1216 * 1e-4 # Angstrom to microns
    lyman_break_observed = lyman_break_rest*(z+1) 
    lyman_break_observed = [lyman_break_observed] * len(y)
    
    return lyman_break_observed, y
lyman_break_observed, y = Lyman_break_position(z)

telescopes = ['VISTA', 'Subaru', 'Euclid']
vista_filters = ['Y', 'J', 'H', 'Ks']
subaru_filters = ['B', 'V', 'gplus', 'rplus', 'iplus']
euclid_filters = ['H', 'Y', 'J', 'VIS']

def get_data_files(telescopes, vista_filters, subaru_filters, euclid_filters):
    '''
    Loops over each filter name in the list of filters for each telescope and 
    creates a dictionary of filter names and their respective files.

    INPUT(s)
    telescopes(list): list of telescopes
    vista_filters(list): list of VISTA filters
    subaru_filters(list): list of Subaru filters
    euclid_filters(list): list of Euclid filters
    
    OUTPUT(s)
    vista_files(dict): dictionary of VISTA filters and their respective files
    subaru_files(dict): dictionary of Subaru filters and their respective files
    euclid_files(dict): dictionary of Euclid filters and their respective files
    '''
    for telescope in telescopes:
        if telescope == 'VISTA':
            vista_files = {}
            for vista_filter in vista_filters:
                vista_files[vista_filter] = f'/Users/user/Documents/filter_curves/data_files/VISTA_VIRCAM_{vista_filter}.csv'
        elif telescope == 'Subaru':
            subaru_files = {}
            for subaru_filter in subaru_filters:
                subaru_files[subaru_filter] = f'/Users/user/Documents/filter_curves/data_files/subaru_{subaru_filter}.csv'
        elif telescope == 'Euclid':
            euclid_files = {}
            for euclid_filter in euclid_filters:
                euclid_files[euclid_filter] = f'/Users/user/Documents/filter_curves/data_files/Euclid_{euclid_filter}.csv'
        else:
            print('No such telescope-filter combination')
    return vista_files, subaru_files, euclid_files

vista_files, subaru_files, euclid_files = get_data_files(telescopes, vista_filters, subaru_filters, euclid_filters)

# extract data from csv files
def get_wl_transmission(file):
    '''
    Extracts wavelength and transmission data from csv files converting Angstrom to microns.

    INPUT(s)
    file(str): path to csv file


    OUTPUT(s)
    wavelength_list(list): list of wavelengths
    transmission_list(list): list of transmission values
    '''
    wavelength_list = []
    transmission_list = []
    with open(file,'r', encoding='utf-8-sig') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            wavelength_list.append(float(row[0])* 1e-4) # Angstrom to microns
            transmission_list.append(100*float(row[1])) # convert to percentage
    return wavelength_list, transmission_list

for telescope in telescopes:
    if telescope == 'VISTA':
        for vista_filter in vista_filters:
            vista_files[vista_filter] = get_wl_transmission(vista_files[vista_filter])
    elif telescope == 'Subaru':
        for subaru_filter in subaru_filters:
            subaru_files[subaru_filter] = get_wl_transmission(subaru_files[subaru_filter])
    elif telescope == 'Euclid':
        for euclid_filter in euclid_filters:
            euclid_files[euclid_filter] = get_wl_transmission(euclid_files[euclid_filter])
base=0.5
for telescope in telescopes:
    if telescope == 'VISTA':
        for vista_filter in vista_filters:
            plt.fill_between(vista_files[vista_filter][0], vista_files[vista_filter][1],
                              color=PuRd(base+0.75*vista_filters.index(vista_filter)), alpha=0.2, label=f'Vista {vista_filter}')
    elif telescope == 'Subaru':
        for subaru_filter in subaru_filters:
            plt.fill_between(subaru_files[subaru_filter][0], subaru_files[subaru_filter][1],
                              color=YlOrRd(base+0.75*subaru_filters.index(subaru_filter)), alpha=0.2, label=f'Subaru {subaru_filter}')
    elif telescope == 'Euclid':
        for euclid_filter in euclid_filters:
            plt.fill_between(euclid_files[euclid_filter][0], euclid_files[euclid_filter][1],
                              color=oranges(base+0.75*euclid_filters.index(euclid_filter)), alpha=0.2, label=f'Euclid {euclid_filter}')
plt.plot(LBG_wavelength, LBG_flux, label='LBG Spectrum', color='gray')
plt.plot(lyman_break_observed, y, label='Lyman Break', color='magenta', linestyle='--')
plt.title('Filter Curves')
plt.xlabel(r'Wavelength $(\mu)$')
plt.ylabel('Transmission (%)')
plt.xlim(0.3, 2.5)
plt.ylim(0,110)
plt.legend()
plt.show()


"""

for telescope in telescopes:
    if telescope == 'VISTA':
        for vista_filter in vista_filters:
            plt.fill_between(vista_files[vista_filter][0], vista_files[vista_filter][1],
                              color=nipy_spectral(0.3-((vista_filters.index(vista_filter)+30)/100)), alpha=0.2, label=f'Vista {vista_filter}')
    elif telescope == 'Subaru':
        for subaru_filter in subaru_filters:
            plt.fill_between(subaru_files[subaru_filter][0], subaru_files[subaru_filter][1],
                              color=nipy_spectral(0.85-((subaru_filters.index(subaru_filter)+30)/100)), alpha=0.2, label=f'Subaru {subaru_filter}')
    elif telescope == 'Euclid':
        for euclid_filter in euclid_filters:
            plt.fill_between(euclid_files[euclid_filter][0], euclid_files[euclid_filter][1],
                              color=nipy_spectral(0.1-((euclid_filters.index(euclid_filter)+30)/100)), alpha=0.2, label=f'Euclid {euclid_filter}')
plt.title('Filter Curves')
"""