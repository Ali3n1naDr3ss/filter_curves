# works nice! colours not very clear

# imports and colour spectra from 
# https://matplotlib.org/stable/gallery/color/colormap_reference.html
import numpy as np
import matplotlib.pyplot as plt
import csv
from labellines import labelLines

# Define constants 
z_REBELS = 6.496 # rredshift
c = 3e18 # speed of light (Angstrom/s)
base = 0.5 # base for shaded plot areas - ensures no plots are white

telescopes = ['VISTA', 'Subaru', 'Euclid']
vista_filters = ['Y', 'J', 'H', 'Ks']
subaru_filters = ['B', 'V', 'gplus', 'rplus', 'iplus']
euclid_filters = ['H', 'Y', 'J', 'VIS']

def LBG_spectrum_prep():
    LBG_spectrum = '/Users/user/Documents/filter_curves/data_files/REBELS-05_1dspec_wRMS.dat'
    LBG_wavelength = []
    LBG_flux = []
    with open(LBG_spectrum,'r', encoding='utf-8-sig') as file:
        lines = file.readlines()[2:]
        for line in lines:
            columns = line.split()
            LBG_wavelength.append(float(columns[0])) # microns
            LBG_flux.append(float(columns[1])) # nJy
    
    def flux_per_wavelength(LBG_flux, LBG_wavelength):
        '''
        Convets nJy to erg s^-1 cm^-2 Å^-1
        '''
        flux_per_wavelength = []
        for LBG_flux, LBG_wavelength in zip(LBG_flux, LBG_wavelength):
            flux_per_A = (LBG_flux * c / (LBG_wavelength**2) * 1e-21)
            flux_per_wavelength.append(flux_per_A)
        return flux_per_wavelength

    LBG_flux = flux_per_wavelength(LBG_flux, LBG_wavelength)
    LBG_flux_scaled = []
    for i in range(len(LBG_flux)):
        LBG_flux[i] = LBG_flux[i] * 40 # scaling factor
        LBG_flux_scaled.append(LBG_flux[i])
    LBG_flux = LBG_flux_scaled

    return LBG_wavelength, LBG_flux

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

def remove_spectrum_noise(lyman_break_observed, LBG_wavelength, LBG_flux):
    '''
    Removes the noise from the short-wavelength regime of the provided spectrum.

    INPUT(s)
    lyman_break_observed(float): observed wavelength of Lyman break in microns
    LBG_wavelength(list): list of wavelengths measurements for each sample point on spectrum
    LBG_flux(list): list of flux measurements for each sample point on spectrum

    OUTPUT(s)
    noise_removed_wavelength(list): list of wavelengths with noise removed
    noise_removed_flux(list): list of flux with noise removed
    '''
    noise_removed_wavelength = []
    noise_removed_flux = []
    
    for i in range(len(LBG_wavelength)):
        if LBG_wavelength[i] >= lyman_break_observed[0]:
            noise_removed_wavelength.append(LBG_wavelength[i])
            noise_removed_flux.append(LBG_flux[i])
    return noise_removed_wavelength, noise_removed_flux

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

def sort_filter_data(vista_files, subaru_files, euclid_files):
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
    return vista_files, subaru_files, euclid_files

def filter_curve_plotter(vista_files, subaru_files, euclid_files, LBG_wavelength, LBG_flux, lyman_break_observed, y):
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
            labelled_lines = []
            ###x#vals = []
            for euclid_filter in euclid_filters:
                line, = plt.plot(euclid_files[euclid_filter][0], euclid_files[euclid_filter][1], label=f'Euclid {euclid_filter}',
                          color=oranges(base+0.75*euclid_filters.index(euclid_filter)))
                labelled_lines.append(line)
                ##x#data = line.get_xdata()
                #x#midpoint = (xdata[0] + xdata[-1]) / 2
                #xvals.append(xmidpoint)
                plt.fill_between(euclid_files[euclid_filter][0], euclid_files[euclid_filter][1],
                                color=oranges(base+0.75*euclid_filters.index(euclid_filter)), alpha=0.2, label=f'Euclid {euclid_filter}')
    
    plt.plot(LBG_wavelength, LBG_flux, label=f'z = {z_REBELS} LBG Spectrum', color='gray')
    plt.plot(lyman_break_observed, y, label=f'Lyman Break at $z =$ {z_REBELS}', color='magenta', linestyle='--')
    #labelled_lines = labelled_lines.append(plt.plot(LBG_wavelength, LBG_flux, label=f'$z =${z_REBELS} LBG Spectrum', color='gray'))

    plt.title('Filter Curves')
    plt.xlabel(r'Observed Wavelength $(\mu)$')
    plt.ylabel('Transmission (%)')
    plt.xlim(0.3, 3.0)
    plt.ylim(0,110)
    #plt.legend(loc='best')
    labelLines(labelled_lines, align=True) #xvals=xmidpoint, zorder=2.5)
    plt.show()

PuRd = plt.get_cmap('PuRd')
YlOrRd = plt.get_cmap('PuBuGn')
oranges = plt.get_cmap('Oranges')


LBG_wavelength, LBG_flux = LBG_spectrum_prep()
lyman_break_observed, y = Lyman_break_position(z_REBELS)
LBG_wavelength, LBG_flux = remove_spectrum_noise(lyman_break_observed, LBG_wavelength, LBG_flux)
vista_files, subaru_files, euclid_files = get_data_files(telescopes, vista_filters, subaru_filters, euclid_filters)
vista_files, subaru_files, euclid_files = sort_filter_data(vista_files, subaru_files, euclid_files)
filter_curve_plotter(vista_files, subaru_files, euclid_files, LBG_wavelength, LBG_flux, lyman_break_observed, y)

