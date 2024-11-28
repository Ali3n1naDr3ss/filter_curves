# individual filter curve plots look good
# Filter plots + LBG spectrum: need finishing

import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.gridspec import GridSpec

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

y = np.linspace(-100, 1500, 1000)
lyman_break_rest = 1216 * 1e-4 # Angstrom to microns
lyman_break_observed = lyman_break_rest*(z+1) 
lyman_break_observed = [lyman_break_observed] * len(y)

# VISTA-VIRCAM
vista_y = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_Y.dat') 
vista_wavelength_y = vista_y[:,0] * 1e-4
vista_flux_y = vista_y[:,1]
vista_j = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_J.dat') 
vista_wavelength_j = vista_j[:,0] * 1e-4
vista_flux_j = vista_j[:,1]
vista_h = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_H.dat')
vista_wavelength_h = vista_h[:,0] * 1e-4
vista_flux_h = vista_h[:,1]
vista_Ks = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_Ks.dat')
vista_wavelenght_Ks = vista_Ks[:,0] * 1e-4
vista_flux_Ks = vista_Ks[:,1]

# Subaru and Euclid
def csv_to_list(file, wavelength_list, transmission_list):
    with open(file,'r', encoding='utf-8-sig') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            wavelength_list.append(float(row[0])* 1e-4) # Angstrom to microns)
            transmission_list.append(100*float(row[1])) # convert to percentage
    return wavelength_list, transmission_list

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

plt.figure(figsize=(15,10))
plt.subplot(2,1,1)
"""
ax1 = plt.subplot(gs[0, 0])
ax1.plot(vista_wavelength_y, vista_flux_y, label='VISTA-Y', color=nipy_spectral(0.1))
ax1.plot(vista_wavelength_j, vista_flux_j, label='VISTA-J', color=nipy_spectral(0.2))
ax1.plot(vista_wavelength_h, vista_flux_h, label='VISTA-H', color=nipy_spectral(0.3))
ax1.plot(vista_wavelenght_Ks, vista_flux_Ks, label='VISTA-Ks', color=nipy_spectral(0.4))
ax1.set_title('VISTA-VIRCAM Filters')
ax1.set_ylabel(r'Transmission ($\%$)')
ax1.set_xlabel(r'Wavelength ($\mu$)')
ax1.set_ylim(0,110)
ax1.legend(loc='lower right')

ax2 = plt.subplot(gs[0, 1])
ax2.plot(subaru_B_wavelength, subaru_B_transmission, label='Subaru B', color=nipy_spectral(0.4))
ax2.plot(subaru_V_wavelength, subaru_V_transmission, label='Subaru V', color=nipy_spectral(0.45))
ax2.plot(subaru_gplus_wavelength, subaru_gplus_transmission, label='Subaru $g^+$', color=nipy_spectral(0.5))
ax2.plot(subaru_rplus_wavelength, subaru_rplus_transmission, label='Subaru $r^+$', color=nipy_spectral(0.6))
ax2.plot(subaru_iplus_wavelength, subaru_iplus_transmission, label='Subaru $i^+$', color=nipy_spectral(0.65))
ax2.set_xlabel(r'Wavelength ($\mu$)')
ax2.set_title('Subaru Optical Filters')
ax2.set_ylim(0,110)
ax2.legend(loc='lower right')

ax3 = plt.subplot(gs[0, 2])
ax3.plot(euclid_VIS_wavelength, euclid_VIS_transmission, label='Euclid VIS', color=nipy_spectral(0.7))
ax3.plot(euclid_Y_wavelength, euclid_Y_transmission, label='Euclid Y', color=nipy_spectral(0.75))
ax3.plot(euclid_J_wavelength, euclid_J_transmission, label='Euclid J', color=nipy_spectral(0.8))
ax3.plot(euclid_H_wavelength, euclid_H_transmission, label='Euclid H', color=nipy_spectral(0.85))
ax3.set_xlabel(r'Wavelength ($\mu$)')
ax3.set_ylim(0,110)
ax3.legend(loc='lower right')
ax3.set_title('Euclid Filters')


vista_flux_y = flux_per_wavelength(vista_flux_y, vista_wavelength_y)
vista_flux_j = flux_per_wavelength(vista_flux_j, vista_wavelength_j)
vista_flux_h = flux_per_wavelength(vista_flux_h, vista_wavelength_h)
vista_flux_Ks = flux_per_wavelength(vista_flux_Ks, vista_wavelenght_Ks)

subaru_B_wavelength = flux_per_wavelength(subaru_B_transmission, subaru_B_wavelength)
subaru_V_wavelength = flux_per_wavelength(subaru_V_transmission, subaru_V_wavelength)
subaru_gplus_wavelength = flux_per_wavelength(subaru_gplus_transmission, subaru_gplus_wavelength)
subaru_rplus_wavelength = flux_per_wavelength(subaru_rplus_transmission, subaru_rplus_wavelength)
subaru_iplus_wavelength = flux_per_wavelength(subaru_iplus_transmission, subaru_iplus_wavelength)

euclid_VIS_wavelength = flux_per_wavelength(euclid_VIS_transmission, euclid_VIS_wavelength)
euclid_Y_wavelength = flux_per_wavelength(euclid_Y_transmission, euclid_Y_wavelength)
euclid_J_wavelength = flux_per_wavelength(euclid_J_transmission, euclid_J_wavelength)
euclid_H_wavelength = flux_per_wavelength(euclid_H_transmission, euclid_H_wavelength)
"""
 
#ax4.plot(LBG_wavelength, LBG_flux, label='LBG Spectrum', color=nipy_spectral(1))
plt.plot(vista_wavelength_y, vista_flux_y, label='VISTA-Y', color=nipy_spectral(0.1))
plt.plot(vista_wavelength_j, vista_flux_j, label='VISTA-J', color=nipy_spectral(0.15))
plt.plot(vista_wavelength_h, vista_flux_h, label='VISTA-H', color=nipy_spectral(0.2))
plt.plot(vista_wavelenght_Ks, vista_flux_Ks, label='VISTA-Ks', color=nipy_spectral(0.25))
plt.plot(subaru_B_wavelength, subaru_B_transmission, label='Subaru B', color=nipy_spectral(0.4))
plt.plot(subaru_V_wavelength, subaru_V_transmission, label='Subaru V', color=nipy_spectral(0.45))
plt.plot(subaru_gplus_wavelength, subaru_gplus_transmission, label='Subaru $g^+$', color=nipy_spectral(0.5))
plt.plot(subaru_rplus_wavelength, subaru_rplus_transmission, label='Subaru $r^+$', color=nipy_spectral(0.55))
plt.plot(subaru_iplus_wavelength, subaru_iplus_transmission, label='Subaru $i^+$', color=nipy_spectral(0.6))
plt.plot(euclid_VIS_wavelength, euclid_VIS_transmission, label='Euclid VIS', color=nipy_spectral(0.7))
plt.plot(euclid_Y_wavelength, euclid_Y_transmission, label='Euclid Y', color=nipy_spectral(0.75))
plt.plot(euclid_J_wavelength, euclid_J_transmission, label='Euclid J', color=nipy_spectral(0.8))
plt.plot(euclid_H_wavelength, euclid_H_transmission, label='Euclid H', color=nipy_spectral(0.85))
plt.xlim(0.0, 2.0)
plt.ylim(0,110)
plt.legend(loc='upper right')
plt.ylabel(r'Transmission ($\%$)')
plt.xlabel(r'Wavelength ($\mu$)')
plt.title('All Filters')

plt.subplot(2,1,2)
plt.plot(LBG_wavelength, LBG_flux, label='(REBELS-05) LBG Spectrum', color=nipy_spectral(1))
plt.plot(lyman_break_observed ,y, label='Lyman Limit', color='purple')
plt.ylim(0,2.5)
plt.xlabel(r'Wavelength ($\mu$)')
plt.ylabel(r'Flux ($10^{-19} erg s^{-1} cm^{-2} Å^{-1}$)')
plt.legend(loc='upper right')
plt.title(f'(LBG Spectrum and Lyman Limit at z = {z})')

plt.legend()
plt.tight_layout()
plt.show()