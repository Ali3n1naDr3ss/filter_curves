import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.gridspec import GridSpec
nipy_spectral = plt.get_cmap('nipy_spectral')

LBG_spectrum = '/Users/user/Documents/filter_curves/data_files/REBELS-05_1dspec_wRMS.dat'
# CISTA-VIRCAM
vista_y = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_Y.dat') 
vista_j = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_J.dat')
vista_h = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_H.dat')
vista_Ks = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_Ks.dat')

# Subaru - Optical Filters
subaru_B = '/Users/user/Documents/filter_curves/data_files/subaru_B.csv'
subaru_V = '/Users/user/Documents/filter_curves/data_files/subaru_V.csv'
subaru_gplus = '/Users/user/Documents/filter_curves/data_files/subaru_gplus.csv'
subaru_rplus = '/Users/user/Documents/filter_curves/data_files/subaru_rplus.csv'
subaru_iplus = '/Users/user/Documents/filter_curves/data_files/subaru_iplus.csv'

# Euclid Filters
euclid_H = '/Users/user/Documents/filter_curves/data_files/Euclid_H.csv'
euclid_Y = '/Users/user/Documents/filter_curves/data_files/Euclid_Y.csv'
euclid_J = '/Users/user/Documents/filter_curves/data_files/Euclid_J.csv'
euclid_VIS = '/Users/user/Documents/filter_curves/data_files/Euclid_VIS.csv'

files = [LBG_spectrum, vista_y, vista_j, vista_h, vista_Ks, subaru_B, subaru_V, subaru_gplus, subaru_rplus, subaru_iplus, euclid_H, euclid_Y, euclid_J, euclid_VIS]

# VISTA Wavelengths
def vista_wavelength_converter(file):
    ''' Converts Angstrom to microns fior raw VISTA data'''
    vista_wavelengths = []
    vista_wavelength = file[:,0] * 1e-4
    vista_wavelengths.append(vista_wavelength)
    return vista_wavelengths


# Prepare Subaru and Euclid data
def csv_to_list(file):
    ''' Converts csv files to lists for Subaru and Euclid data'''
    wavelength_list = []
    transmission_list = []
    with open(file,'r', encoding='utf-8-sig') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            wavelength_list.append(float(row[0])* 1e-4) # Angstrom to microns)
            transmission_list.append(100*float(row[1])) # convert to percentage
    return wavelength_list, transmission_list

vista_filters = [vista_y, vista_j, vista_h, vista_Ks]
for vista_filter in vista_filters:
    vista_wavelength = vista_wavelength_converter(vista_filter)
    #print(vista_wavelength)
    vista_flux = vista_filter[:,1]

subaru_filters = [subaru_B, subaru_V, subaru_gplus, subaru_rplus, subaru_iplus]
for subaru_filter in subaru_filters:
    subaru_wavelength, subaru_transmission = csv_to_list(subaru_filter)

euclid_filters = [euclid_H, euclid_Y, euclid_J, euclid_VIS]
for euclid_filter in euclid_filters:
    euclid_wavelength, euclid_transmission = csv_to_list(euclid_filter,)

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

y = np.linspace(-100, 1500, 1000)
lyman_break_rest = 1216 * 1e-4 # Angstrom to microns
lyman_break_observed = lyman_break_rest*(z+1) 
lyman_break_observed = [lyman_break_observed] * len(y)


plt.figure(figsize=(15,10))
gs = GridSpec(3, 3)

ax1 = plt.subplot(gs[0, 0])
for vista_filter in vista_filters:
    ax1.plot(vista_wavelength, vista_flux, label='VISTA-?', color=nipy_spectral(vista_filter + 0.1))
ax1.set_title('VISTA-VIRCAM Filters')
ax1.set_ylabel(r'Transmission ($\%$)')
ax1.set_xlabel(r'Wavelength ($\mu$)')
ax1.set_ylim(0,110)
ax1.legend(loc='lower right')

ax2 = plt.subplot(gs[0, 1])
for subaru_filter in subaru_filters:
    ax2.plot(subaru_wavelength, subaru_transmission, label='Subaru ?', color=nipy_spectral(subaru_filter + 0.05))
ax2.set_xlabel(r'Wavelength ($\mu$)')
ax2.set_title('Subaru Optical Filters')
ax2.set_ylim(0,110)
ax2.legend(loc='lower right')

ax3 = plt.subplot(gs[0, 2])
for euclid_filter in euclid_filters:
    ax3.plot(euclid_wavelength, euclid_transmission, label='Euclid ?', color=nipy_spectral(euclid_filter + 0.05))
ax3.set_xlabel(r'Wavelength ($\mu$)')
ax3.set_ylim(0,110)
ax3.legend(loc='lower right')
ax3.set_title('Euclid Filters')

ax4 = plt.subplot(gs[1:, :])    
ax4.plot(LBG_wavelength, LBG_flux, label='LBG Spectrum', color=nipy_spectral(1))
for vista_filter in vista_filters:
    ax4.plot(vista_wavelength[vista_filter], vista_flux, label='VISTA-?', color=nipy_spectral(vista_filter + 0.1))
for subaru_filter in subaru_filters:
    ax4.plot(subaru_wavelength[subaru_filter], subaru_transmission, label='Subaru ?', color=nipy_spectral(subaru_filter + 0.05))
for euclid_filter in euclid_filters:
    ax4.plot(euclid_wavelength[euclid_filter], euclid_transmission, label='Euclid ?', color=nipy_spectral(euclid_filter + 0.05))
ax4.plot(lyman_break_observed ,y, label='Lyman Limit', color='purple')

ax4.plot(lyman_break_observed ,y, label='Lyman Limit', color='purple')
ax4.set_ylim(0,2.5)
ax4.set_xlim(0.65, 5.5)
ax4.set_xlabel(r'Wavelength ($\mu$)')
ax4.set_ylabel(r'Flux (10^{-19} $erg s^{-1} cm^{-2} Å^{-1}$)')
ax4.legend(loc='upper right')
ax4.set_title('All Filters')

plt.tight_layout()
plt.show()