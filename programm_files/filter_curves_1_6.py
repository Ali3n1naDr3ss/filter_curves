# works well :)
# LBG spectrm and LB plotted
# Filter curves for VISTA< Subaru, Euclid shaded

import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.gridspec import GridSpec
from labellines import labelLines

nipy_spectral = plt.get_cmap('nipy_spectral')

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

LBG_flux = LBG_flux_scaled

y = np.linspace(-100, 1500, 1000)
lyman_break_rest = 1216 * 1e-4 # Angstrom to microns
lyman_break_observed = lyman_break_rest*(z+1) 
lyman_break_observed = [lyman_break_observed] * len(y)

# open csv files
def csv_to_list(file, wavelength_list, transmission_list):
    with open(file,'r', encoding='utf-8-sig') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            wavelength_list.append(float(row[0])* 1e-4) # Angstrom to microns)
            transmission_list.append(100*float(row[1])) # convert to percentage
    return wavelength_list, transmission_list

# VISTA-VIRCAM
vista_Y = '/Users/user/Documents/filter_curves/data_files/VISTA_VIRCAM_Y.csv'
vista_Y_wavelength = []
vista_Y_transmission = []
vista_Y_wavelength, vista_Y_transmission = csv_to_list(vista_Y, vista_Y_wavelength, vista_Y_transmission)

vista_J = '/Users/user/Documents/filter_curves/data_files/VISTA_VIRCAM_J.csv'
vista_J_wavelength = []
vista_J_transmission = []
vista_J_wavelength, vista_J_transmission = csv_to_list(vista_J, vista_J_wavelength, vista_J_transmission)

vista_H = '/Users/user/Documents/filter_curves/data_files/VISTA_VIRCAM_H.csv'
vista_H_wavelength = []
vista_H_transmission = []
vista_H_wavelength, vista_H_transmission = csv_to_list(vista_H, vista_H_wavelength, vista_H_transmission)

vista_Ks = '/Users/user/Documents/filter_curves/data_files/VISTA_VIRCAM_Ks.csv'
vista_Ks_wavelength = []
vista_Ks_transmission = []
vista_Ks_wavelength, vista_Ks_transmission = csv_to_list(vista_Ks, vista_Ks_wavelength, vista_Ks_transmission)

subaru_B = '/Users/user/Documents/filter_curves/data_files/subaru_B.csv'
subaru_B_wavelength = []
subaru_B_transmission = []
subaru_B_wavelength, subaru_B_transmission = csv_to_list(subaru_B, subaru_B_wavelength, subaru_B_transmission)

subaru_V = '/Users/user/Documents/filter_curves/data_files/subaru_V.csv'
subaru_V_wavelength = []
subaru_V_transmission = []
subaru_V_wavelength, subaru_V_transmission = csv_to_list(subaru_V, subaru_V_wavelength, subaru_V_transmission)

subaru_gplus = '/Users/user/Documents/filter_curves/data_files/subaru_gplus.csv'
subaru_gplus_wavelength = []
subaru_gplus_transmission = []
subaru_gplus_wavelength, subaru_gplus_transmission = csv_to_list(subaru_gplus, subaru_gplus_wavelength, subaru_gplus_transmission)

subaru_rplus = '/Users/user/Documents/filter_curves/data_files/subaru_rplus.csv'
subaru_rplus_wavelength = []
subaru_rplus_transmission = []
subaru_rplus_wavelength, subaru_rplus_transmission = csv_to_list(subaru_rplus, subaru_rplus_wavelength , subaru_rplus_transmission)

subaru_iplus = '/Users/user/Documents/filter_curves/data_files/subaru_iplus.csv'
subaru_iplus_wavelength = []
subaru_iplus_transmission = []
subaru_iplus_wavelength, subaru_iplus_transmission = csv_to_list(subaru_iplus, subaru_iplus_wavelength, subaru_iplus_transmission)

euclid_H = '/Users/user/Documents/filter_curves/data_files/Euclid_H.csv'
euclid_H_wavelength = []
euclid_H_transmission = []
euclid_H_wavelength, euclid_H_transmission = csv_to_list(euclid_H, euclid_H_wavelength, euclid_H_transmission)  

euclid_Y = '/Users/user/Documents/filter_curves/data_files/Euclid_Y.csv'
euclid_Y_wavelength = []
euclid_Y_transmission = []
euclid_Y_wavelength, euclid_Y_transmission = csv_to_list(euclid_Y, euclid_Y_wavelength, euclid_Y_transmission)  

euclid_J = '/Users/user/Documents/filter_curves/data_files/Euclid_J.csv'
euclid_J_wavelength = []
euclid_J_transmission = []
euclid_J_wavelength, euclid_J_transmission = csv_to_list(euclid_J, euclid_J_wavelength, euclid_J_transmission)  

euclid_VIS = '/Users/user/Documents/filter_curves/data_files/Euclid_VIS.csv'
euclid_VIS_wavelength = []
euclid_VIS_transmission = []
euclid_VIS_wavelength, euclid_VIS_transmission = csv_to_list(euclid_VIS, euclid_VIS_wavelength, euclid_VIS_transmission)

plt.plot(LBG_wavelength, LBG_flux, label='(REBELS-05) LBG Spectrum', color=nipy_spectral(1))
plt.plot(lyman_break_observed ,y, label='Lyman Limit', color='purple', linestyle='--')
plt.plot(vista_Y_wavelength, vista_Y_transmission, label='VISTA-Y', color=nipy_spectral(0.1))
plt.fill_between(vista_Y_wavelength, vista_Y_transmission, color=nipy_spectral(0.1), alpha=0.2)
plt.plot(vista_J_wavelength, vista_J_transmission, label='VISTA-J', color=nipy_spectral(0.15))
plt.fill_between(vista_J_wavelength, vista_J_transmission, color=nipy_spectral(0.15), alpha=0.2)
plt.plot(vista_J_wavelength, vista_J_transmission, label='VISTA-J', color=nipy_spectral(0.15))
plt.fill_between(vista_J_wavelength, vista_J_transmission, color=nipy_spectral(0.15), alpha=0.2)
plt.plot(vista_H_wavelength, vista_H_transmission, label='VISTA-H', color=nipy_spectral(0.2))
plt.fill_between(vista_H_wavelength, vista_H_transmission, color=nipy_spectral(0.2), alpha=0.2)
plt.plot(vista_Ks_wavelength, vista_Ks_transmission, label='VISTA-Ks', color=nipy_spectral(0.25))
plt.fill_between(vista_Ks_wavelength, vista_Ks_transmission, color=nipy_spectral(0.25), alpha=0.2)
plt.plot(subaru_B_wavelength, subaru_B_transmission, label='Subaru B', color=nipy_spectral(0.4))
plt.fill_between(subaru_B_wavelength, subaru_B_transmission, color=nipy_spectral(0.4), alpha=0.2)
plt.plot(subaru_V_wavelength, subaru_V_transmission, label='Subaru V', color=nipy_spectral(0.45))
plt.fill_between(subaru_V_wavelength, subaru_V_transmission, color=nipy_spectral(0.45), alpha=0.2)
plt.plot(subaru_gplus_wavelength, subaru_gplus_transmission, label='Subaru $g^+$', color=nipy_spectral(0.5))
plt.fill_between(subaru_gplus_wavelength, subaru_gplus_transmission, color=nipy_spectral(0.5), alpha=0.2)
plt.plot(subaru_rplus_wavelength, subaru_rplus_transmission, label='Subaru $r^+$', color=nipy_spectral(0.55))
plt.fill_between(subaru_rplus_wavelength, subaru_rplus_transmission, color=nipy_spectral(0.55), alpha=0.2)
plt.plot(subaru_iplus_wavelength, subaru_iplus_transmission, label='Subaru $i^+$', color=nipy_spectral(0.6))
plt.fill_between(subaru_iplus_wavelength, subaru_iplus_transmission, color=nipy_spectral(0.6), alpha=0.2)
plt.plot(euclid_VIS_wavelength, euclid_VIS_transmission, label='Euclid VIS', color=nipy_spectral(0.7))
plt.fill_between(euclid_VIS_wavelength, euclid_VIS_transmission, color=nipy_spectral(0.7), alpha=0.2)
plt.plot(euclid_Y_wavelength, euclid_Y_transmission, label='Euclid Y', color=nipy_spectral(0.75))
plt.fill_between(euclid_Y_wavelength, euclid_Y_transmission, color=nipy_spectral(0.75), alpha=0.2)
plt.plot(euclid_J_wavelength, euclid_J_transmission, label='Euclid J', color=nipy_spectral(0.8))
plt.fill_between(euclid_J_wavelength, euclid_J_transmission, color=nipy_spectral(0.8), alpha=0.2)
plt.plot(euclid_H_wavelength, euclid_H_transmission, label='Euclid H', color=nipy_spectral(0.85))
plt.fill_between(euclid_H_wavelength, euclid_H_transmission, color=nipy_spectral(0.85), alpha=0.2)
plt.xlim(0.3, 3.0)
plt.ylim(0,110)
#plt.legend(loc='upper right')
plt.ylabel(r"Transmission $(\%)$")
plt.xlabel(r"Wavelength $(\mu)$")
plt.title('All Filters')
labelLines(plt.gca().get_lines(), zorder=5)
plt.tight_layout()
plt.show()