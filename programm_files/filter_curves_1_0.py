import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.gridspec import GridSpec

nipy_spectral = plt.get_cmap('nipy_spectral')

# VISTA-VIRCAM
vista_y = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_Y.dat')
vista_j = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_J.dat')
vista_h = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_H.dat')
vistta_Ks = np.loadtxt('/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_Ks.dat')

# Subaru - Optical Filters
def csv_to_list(file, wavelength_list, transmission_list):
    with open(file,'r', encoding='utf-8-sig') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            wavelength_list.append(float(row[0]))
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
gs = GridSpec(3, 3)

ax1 = plt.subplot(gs[0, 0])
ax1.plot(vista_y[:,0], vista_y[:,1], label='VISTA-Y', color=nipy_spectral(0.1))
ax1.plot(vista_j[:,0], vista_j[:,1], label='VISTA-J', color=nipy_spectral(0.2))
ax1.plot(vista_h[:,0], vista_h[:,1], label='VISTA-H', color=nipy_spectral(0.3))
ax1.plot(vistta_Ks[:,0], vistta_Ks[:,1], label='VISTA-Ks', color=nipy_spectral(0.4))
ax1.set_title('VISTA-VIRCAM Filters')
ax1.set_ylabel('Transmission (%)')
ax1.set_xlabel('Wavelength (Angstrom)')
ax1.set_ylim(0,110)
ax1.legend(loc='lower right')

ax2 = plt.subplot(gs[0, 1])
ax2.plot(subaru_B_wavelength, subaru_B_transmission, label='Subaru B', color=nipy_spectral(0.4))
ax2.plot(subaru_V_wavelength, subaru_V_transmission, label='Subaru V', color=nipy_spectral(0.45))
ax2.plot(subaru_gplus_wavelength, subaru_gplus_transmission, label='Subaru $g^+$', color=nipy_spectral(0.5))
ax2.plot(subaru_rplus_wavelength, subaru_rplus_transmission, label='Subaru $r^+$', color=nipy_spectral(0.6))
ax2.plot(subaru_iplus_wavelength, subaru_iplus_transmission, label='Subaru $i^+$', color=nipy_spectral(0.65))
ax2.set_xlabel('Wavelength (Angstrom)')
ax2.set_title('Subaru Optical Filters')
ax2.set_ylim(0,110)
ax2.legend(loc='lower right')

ax3 = plt.subplot(gs[0, 2])
ax3.plot(euclid_VIS_wavelength, euclid_VIS_transmission, label='Euclid VIS', color=nipy_spectral(0.7))
ax3.plot(euclid_Y_wavelength, euclid_Y_transmission, label='Euclid Y', color=nipy_spectral(0.75))
ax3.plot(euclid_J_wavelength, euclid_J_transmission, label='Euclid J', color=nipy_spectral(0.8))
ax3.plot(euclid_H_wavelength, euclid_H_transmission, label='Euclid H', color=nipy_spectral(0.85))
ax3.set_xlabel('Wavelength (Angstrom)')
ax3.set_ylim(0,110)
ax3.legend(loc='lower right')
ax3.set_title('Euclid Filters')

ax4 = plt.subplot(gs[1:, :])
ax4.plot(vista_y[:,0], vista_y[:,1], label='VISTA-Y', color=nipy_spectral(0.1))
ax4.plot(vista_j[:,0], vista_j[:,1], label='VISTA-J', color=nipy_spectral(0.15))
ax4.plot(vista_h[:,0], vista_h[:,1], label='VISTA-H', color=nipy_spectral(0.2))
ax4.plot(vistta_Ks[:,0], vistta_Ks[:,1], label='VISTA-Ks', color=nipy_spectral(0.25))
ax4.plot(subaru_B_wavelength, subaru_B_transmission, label='Subaru B', color=nipy_spectral(0.4))
ax4.plot(subaru_V_wavelength, subaru_V_transmission, label='Subaru V', color=nipy_spectral(0.45))
ax4.plot(subaru_gplus_wavelength, subaru_gplus_transmission, label='Subaru $g^+$', color=nipy_spectral(0.5))
ax4.plot(subaru_rplus_wavelength, subaru_rplus_transmission, label='Subaru $r^+$', color=nipy_spectral(0.55))
ax4.plot(subaru_iplus_wavelength, subaru_iplus_transmission, label='Subaru $i^+$', color=nipy_spectral(0.6))
ax4.plot(euclid_VIS_wavelength, euclid_VIS_transmission, label='Euclid VIS', color=nipy_spectral(0.7))
ax4.plot(euclid_Y_wavelength, euclid_Y_transmission, label='Euclid Y', color=nipy_spectral(0.75))
ax4.plot(euclid_J_wavelength, euclid_J_transmission, label='Euclid J', color=nipy_spectral(0.8))
ax4.plot(euclid_H_wavelength, euclid_H_transmission, label='Euclid H', color=nipy_spectral(0.85))
ax4.set_xlabel('Wavelength (Angstrom)')
ax4.set_ylabel('Transmission (%)')
ax4.set_title('All Filters')
ax4.set_ylim(0,110)
ax4.legend(loc='lower right')

plt.tight_layout()
plt.show()