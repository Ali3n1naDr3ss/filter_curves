import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines

# LBGs at other redshifts - use step functions for now
def LBG_step_generator():
    """
    Generates step functions to represent the Lyman break at different redshifts
    """
    redshifts = [7, 8, 9, 10, 11]
    observed_wavelengths = []
    rest_frame_wavelength = 912 * 1e-4  # microns

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
    step_colour_base = 0.5
    y_values = np.linspace(0, 100, 100)
    step_size = 15
    labelled_lines_LBGs = []

    for i in range(len(redshifts)):
        step_label = f'z = {redshifts[i]}'
        break_values = [observed_wavelengths[i]] * len(y_values)
        x_values = np.linspace(observed_wavelengths[i], observed_wavelengths[i] + 10, len(y_values))
        y_const = np.linspace(max(y_values) - i * step_size, max(y_values) - i * step_size, len(y_values))  # -i*step_size to move the step down

        if i == 0:
            line1, = plt.plot(break_values, y_values, color=step_colour(step_colour_base + 0.2 * i))
            line2, = plt.plot(x_values, y_const, color=step_colour(step_colour_base + 0.2 * i))
        else:
            line1, = plt.plot(break_values[:-i * step_size], y_values[:-i * step_size], color=step_colour(step_colour_base + 0.2 * i))  # vertical
            line2, = plt.plot(x_values, y_const, color=step_colour(step_colour_base + 0.2 * i))  # horizontal

        # Append the lines to the list for labeling
        labelled_lines_LBGs.append(line1)
        labelled_lines_LBGs.append(line2)
    # Example x-values for labeling
    xvals = [observed_wavelengths[i] + 5 for i in range(len(redshifts))]

    # Label the lines at specified x-values
    labelLines(labelled_lines_LBGs, xvals=xvals, align=True, zorder=2.5)    
    
    plt.show()

LBG_step_generator()