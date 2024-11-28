import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.gridspec import GridSpec
nipy_spectral = plt.get_cmap('nipy_spectral')

vista_filter_names = ['Y', 'J', 'H', 'Ks']
subaru_filter_names = ['B', 'V', 'gplus', 'rplus', 'iplus']
euclid_filter_names = ['H', 'Y', 'J', 'VIS']
vista_filters = []
subaru_filters = []
euclid_filters = []
for filter_name in vista_filter_names:
    vista_filter = f'/Users/user/Documents/filter_curves/data_files/VISTA_Filters_at80K_forETC_{filter_name}.dat'
    vista_filters.append(vista_filter)

for filter_name in subaru_filter_names:
    subaru_filter = f'/Users/user/Documents/filter_curves/data_files/subaru_{filter_name}.csv'
    subaru_filters.append(subaru_filter)

for filter_name in euclid_filter_names:
    euclid_filter = f'/Users/user/Documents/filter_curves/data_files/Euclid_{filter_name}.csv'
    euclid_filters.append(euclid_filter)


all_filters = vista_filters + subaru_filters + euclid_filters
# retrieved all telescope data
# now convert wavelengths to microns

#all_vista_filters = all_filters[0]
#vista_y = all_vista_filters[0:]
#vista_y_wl = vista_y[:,0]
#print(all_filters[0][0:][:,0]) -> all_filters[telescope index][telescope filter index][wavelength]

for filter in all_filters:
    if 'VISTA' in filter:
        vista_filters = []
        
        vista_filter = np.loadtxt(filter)
        vista_filter_wavelength = vista_filter[:,0] * 1e-4
        vista_filter_flux = vista_filter[:,1]
    elif 'subaru' in filter:
        subaru_filter_wavelength = []
        subaru_filter_flux = []
        with open(filter,'r', encoding='utf-8-sig') as filter_file:
            reader = csv.reader(filter_file)
            for filter_data in reader:
                filter_wavelength = filter_data[0]
    elif 'Euclid' in filter:
        euclid_filter_wavelength = []
        euclid_filter_flux = []
        with open(filter,'r', encoding='utf-8-sig') as filter_file:
            reader = csv.reader(filter_file)
            for filter_data in reader:
                filter_wavelength = filter_data[0]
print(vista_filter_wavelength)   

