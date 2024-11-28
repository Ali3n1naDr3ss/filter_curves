from filter_curves_1_8 import *

LBG_wavelength, LBG_flux = LBG_spectrum_prep()
lyman_break_observed, y = Lyman_break_position(z)
vista_files, subaru_files, euclid_files = get_data_files(telescopes, vista_filters, subaru_filters, euclid_filters)
vista_files, subaru_files, euclid_files = sort_filter_data(vista_files, subaru_files, euclid_files)
filter_curve_plotter(vista_files, subaru_files, euclid_files, LBG_wavelength, LBG_flux, lyman_break_observed, y)