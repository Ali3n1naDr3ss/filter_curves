# based on version 1.7 
# 
# TODO: 
# RECALCULATE BREAK POSITION INSTEP FUNCTION GENERATOR
    # - use the same method as in the LBG_spectrum_prep function

import numpy as np
import matplotlib.pyplot as plt
import csv
from labellines import labelLines


# LBGs at other redshifts - use step functions for now
def LBG_step_generator():
    """
    Generates step functions to represent the Lyman break at different redshifts
    """
    redshifts = [7, 8, 9, 10, 11]
    observed_wavelengths = []
    rest_frame_wavelength = 912*1e-4 # microns

    # calculate observed wl
    def observed_wavelength_calculator(rest_frame_wavelength, z):
        '''
        INPUT(s)
        rest_frame_wavelength(float): rest wavelength in microns
        z(float): redshift

        OUTPUT(s)
        observed_wavelength(float): observed wavelength in microns
        '''
        observed_wavelength = rest_frame_wavelength * (1 + z)
        return observed_wavelength

    for z in redshifts:
        observed_wavelength = observed_wavelength_calculator(rest_frame_wavelength, z)
        observed_wavelengths.append(observed_wavelength)

    # Step function
    step_colour = plt.get_cmap('Reds')
    step_colour_base=0.25  
    y_values = np.linspace(0, 100, 100)
    step_size = 15
    labelled_lines_LBGs = []

    for i in range(len(redshifts)):
        step_label = f'z = {redshifts[i]}'
        break_values = [observed_wavelengths[i]] * len(y_values)
        x_values = np.linspace(observed_wavelengths[i], observed_wavelengths[i] + 10, len(y_values))
        y_const = np.linspace(max(y_values) - i * step_size, max(y_values) - i * step_size, len(y_values))  # -i*step_size to move the step down

        if i == 0:
            plt.plot(break_values, y_values, color=step_colour(step_colour_base + 0.2 * i))
            plt.plot(x_values, y_const, color=step_colour(step_colour_base + 0.2 * i), label=step_label)
        else:
            plt.plot(break_values[:-i * step_size], y_values[:-i * step_size], color=step_colour(step_colour_base + 0.2 * i))  # vertical
            plt.plot(x_values, y_const, color=step_colour(step_colour_base + 0.2 * i), label=step_label)  # horizontal

    labelLines(align=True, xvals=[1.75,1.75,1.75,1.75,1.75], zorder=2.5)

# Define constants 
z_REBELS05 = 6.496 # actual sprectrum redshift
c = 3e18 # speed of light (Angstrom/s)
base = 0.2 # base for shaded plot areas - ensures no plots are white

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

def filter_curve_plotter(vista_files, subaru_files, euclid_files, LBG_wavelength, LBG_flux):
    for telescope in telescopes:
        if telescope == 'VISTA':
            for vista_filter in vista_filters:
                plt.fill_between(vista_files[vista_filter][0], vista_files[vista_filter][1],
                                color=PuRd(base+0.5*vista_filters.index(vista_filter)), alpha=0.1, label=f'Vista {vista_filter}')
                plt.annotate('VISTA filters', (1.65, 90), color=PuRd(0.5))
        elif telescope == 'Subaru':
            for subaru_filter in subaru_filters:
                plt.fill_between(subaru_files[subaru_filter][0], subaru_files[subaru_filter][1],
                                color=YlOrRd(base+0.5*subaru_filters.index(subaru_filter)), alpha=0.1, label=f'Subaru {subaru_filter}')
                plt.annotate('Subaru filters', (0.5, 90), color=YlOrRd(0.5))
        elif telescope == 'Euclid':
            labelled_lines_euclid = []
            for euclid_filter in euclid_filters:
                line, = plt.plot(euclid_files[euclid_filter][0], euclid_files[euclid_filter][1], label=f'Euclid {euclid_filter}',
                          color=oranges(base+0.5*euclid_filters.index(euclid_filter)))
                labelled_lines_euclid.append(line)
                plt.fill_between(euclid_files[euclid_filter][0], euclid_files[euclid_filter][1],
                                color=oranges(base+0.5*euclid_filters.index(euclid_filter)), alpha=0.4, label=f'Euclid {euclid_filter}')
    line, = plt.plot(LBG_wavelength, LBG_flux, label=f'z = {z_REBELS05} LBG Spectrum', color='blue')
    labelled_lines_euclid.append(line)
    plt.title('Filter Curves')
    plt.xlabel(r'Observed Wavelength $(\mu)$')
    plt.ylabel('Transmission (%)')
    plt.xlim(0.3, 2.5)
    plt.ylim(0,110)
    labelLines(labelled_lines_euclid, align=False, xvals=[1.75,1.05,1.4,0.7,2.0], zorder=2.5, outline_width=7)
    plt.show()

PuRd = plt.get_cmap('PuRd')
YlOrRd = plt.get_cmap('PuBuGn')
oranges = plt.get_cmap('Oranges')

LBG_wavelength, LBG_flux = LBG_spectrum_prep()
lyman_break_observed, y = Lyman_break_position(z_REBELS05)
LBG_wavelength, LBG_flux = remove_spectrum_noise(lyman_break_observed, LBG_wavelength, LBG_flux)
vista_files, subaru_files, euclid_files = get_data_files(telescopes, vista_filters, subaru_filters, euclid_filters)
vista_files, subaru_files, euclid_files = sort_filter_data(vista_files, subaru_files, euclid_files)
LBG_step_generator()
filter_curve_plotter(vista_files, subaru_files, euclid_files, LBG_wavelength, LBG_flux)

