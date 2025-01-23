import numpy as np
from numpy import e, pi, sin, cos
import matplotlib.pyplot as plt


def sum_normal_modes(
        string_len_m=1., velocity=1., number_modes=40, 
        source_position_m=0.2, rcv_position_m=0.7, seismogram_duration_s=1.25, 
        nstep=500, tau=.02, plot=True, plot_lines=None
        ):
    """
    Sum normal modes to create a traveling wave and resulting seismogram.
    Transcribed from Stein and Wysession A.8.1 following SW section 2.2

    :param string_len_m: Total length of string L in unit meters
    :param velocity: Speed of the wave in meters per second (m/s)
    :param number_modes: Total number of modes to sum up from n=0
    :param source_position: Where is the input source in unit meters
    :param rc_position_m: Where is the receiver in unit meters
    :param seismogram_duration_s: Duration of the recorded waveform in unit sec.
    :param nstep: Number of time steps
    :param dt: (Optional) length of one time step in seconds
    :param tau: Term controlling the shape of the source pulse
    """
    # Initialize Empty Displacement field for t=0
    time = np.linspace(0, seismogram_duration_s, nstep)
    displacement = np.zeros(nstep)
    modes = []
    
    # Each mode needs to be defined as a term in the summation of n
    for mode in range(1, number_modes + 1):
        # Spatial terms
        source_term = sin(mode * pi * source_position_m / string_len_m)
        rcv_term = sin(mode * pi * rcv_position_m / string_len_m)

        # Time Independent Terms
        omega_n = pi * velocity * mode / string_len_m  # Mode frequency 
        f_omega_n = e ** (-1 * (omega_n * tau)**2 / 4)  # Weighting Term

        # Generate the time-dependent contribution of this mode, m
        mode_disp = source_term * f_omega_n * rcv_term * cos(omega_n * time)
        modes.append(mode_disp)
        displacement += mode_disp


    if plot:
        f, axs = plt.subplots(figsize=(10, 8), nrows=2, height_ratios=[15,1], 
                              sharex=True)
        for m, mode in enumerate(modes):
            # Normalize modes between -0.5 and 0.5 so they fit in a space of one
            # then offset modes by one so they stack on top of each other, 
            # starting at the top with m=0 then plotting downwards
            sep = 1.25  # separation factor between adjacent modes
            scaled_mode = mode / (2 * mode.max()) + (len(modes) - sep * m)
            axs[0].plot(time, scaled_mode, c="k")
            axs[0].text(-0.5, len(modes) - sep * m, f"{m+1:>2}", c="r")  
        
        # Marker for the source, plot a line so we don't need to figure out
        # the corresponding displacement
        if plot_lines:
            for line in plot_lines:
                axs[0].axvline(line, c="r", ls="--")
                axs[1].axvline(line, c="r", ls="--")

        # We don't care about the y axis values, just the representations
        axs[0].get_yaxis().set_ticks([])
        axs[1].get_yaxis().set_ticks([])

        axs[0].set_ylabel("Modes")
        axs[1].plot(time, displacement, c="k")
        axs[1].set_xlabel("Time")
        axs[1].set_ylabel("Sum")

        # Remove empty whitespace between figures
        plt.subplots_adjust(wspace=0, hspace=0)

        plt.show()

    return time, displacement, modes


def mode_sum_over_x(
        time, string_len_m=20., velocity=3., number_modes=40, 
        source_position_m=8, tau=.02, plot=True
        ):
    """
    Sum normal modes to create a traveling wave and resulting seismogram.
    Transcribed from Stein and Wysession A.8.1 following SW section 2.2

    :param string_len_m: Total length of string L in unit meters
    :param velocity: Speed of the wave in meters per second (m/s)
    :param number_modes: Total number of modes to sum up from n=0
    :param source_position: Where is the input source in unit meters
    :param rc_position_m: Where is the receiver in unit meters
    :param seismogram_duration_s: Duration of the recorded waveform in unit sec.
    :param nstep: Number of time steps
    :param dt: (Optional) length of one time step in seconds
    :param tau: Term controlling the shape of the source pulse
    """
    nx = 1000
    # Initialize Empty Displacement field for t=0
    dx = string_len_m / nx
    distance = np.linspace(0, string_len_m, nx)
    displacement = distance * 0
    modes = []
    
    # Each mode needs to be defined as a term in the summation of n
    for mode in range(1, number_modes + 1):
        # Spatial terms
        source_term = sin(mode * pi * source_position_m / string_len_m)
        rcv_term = sin(mode * pi * distance / string_len_m)

        # Time Independent Terms
        omega_n = pi * velocity * mode / string_len_m  # Mode frequency 
        f_omega_n = e ** (-1 * (omega_n * tau)**2 / 4)  # Weighting Term

        # Generate the time-dependent contribution of this mode, m
        mode_disp = source_term * f_omega_n * rcv_term * cos(omega_n * time)
        modes.append(mode_disp)
        displacement += mode_disp

    if plot:
        f, axs = plt.subplots(figsize=(10, 8), nrows=2, height_ratios=[15,1], 
                              sharex=True)
        for m, mode in enumerate(modes):
            # Normalize modes between -0.5 and 0.5 so they fit in a space of one
            # then offset modes by one so they stack on top of each other, 
            # starting at the top with m=0 then plotting downwards
            sep = 1.25  # separation factor between adjacent modes
            scaled_mode = mode / (2 * mode.max()) + (len(modes) - sep * m)
            axs[0].plot(distance, scaled_mode, c="k")
            axs[0].text(-0.5, len(modes) - sep * m, f"{m+1:>2}", c="r")  
        
        # Marker for the source, plot a line so we don't need to figure out
        # the corresponding displacement
        axs[0].axvline(source_position_m, c="r", ls="--")
        axs[1].axvline(source_position_m, c="r", ls="--")

        # We don't care about the y axis values, just the representations
        axs[0].get_yaxis().set_ticks([])
        axs[1].get_yaxis().set_ticks([])

        axs[0].set_ylabel("Modes")
        axs[1].plot(distance, displacement, c="k")
        axs[1].set_xlabel("Distance")
        axs[1].set_ylabel("Sum")

        # Remove empty whitespace between figures
        plt.subplots_adjust(wspace=0, hspace=0)

        plt.show()

    return distance, displacement, modes



if __name__ == "__main__":
    # distance, displacement, modes = mode_sum_over_x(time=1.5, plot=True)
    time, displacement, modes = sum_normal_modes(plot=True, plot_lines=[0.5, 0.9, 1.1])
