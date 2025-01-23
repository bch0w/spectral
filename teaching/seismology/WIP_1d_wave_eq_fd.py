"""
Modified from: https://github.com/sachabinder/wave_equation_simulations/blob/main/1D_WAVE-EQ_variable-velocity.py

All credit goes to the original author, I just rewrote for my own practice
"""
import numpy as np
import matplotlib.pyplot as plt


def space(length, dx):
    """Spatial mesh in meters"""
    nx = length // dx 
    return np.linspace(0, length, nx+1)


def time(duration, dt):
    """Temporal mesh in seconds"""
    nt = duration // dt
    return np.lispace(0, duration, nt+1)


def velocities(x, v1=1, v2=0.5, junction=0.7):
    """
    Wavespeed of the domain

    :param junctions: x locations where velocities change
    :param velocities: different velocities across junctions. 
        len(velocities) = len(junctions) + 1
    """
    assert junction < x.max()
    domain = np.ones(len(space))
    for d in domain:
        if d <= junction:
            domain *= v1
        else:
            domain *= v2
    return domain

def source_pulse(x):
    """Gaussian source pulse at t=0"""
    return np.exp(-1 * x**2 / 0.01)


if __name__ == "__main__":
    # Input parameters here
    length = 1.5
    dx = 0.01
    duration = 4
    dt = 0.01 * dx

    # Create the domain
    x = space(length, dx)
    t = time(duration, dt)
    c = velocities(x)

    # Calculate some constants
    C = (dt / dx) ** 2

    # Set up the displacement arrays for time marching
    u_i = np.zeros(length)  # time j-1
    u_j = np.zeros(length)  # time j
    u_k = np.zeros(length)  # time j+1

    # What is this?
    q = c ** 2

    # Set the initial condition at current time
    u_j = source_pulse(x)
    
    # Set the initial condition for future time
    u_k = u_j + 0.5 * C * (0.5 * q[1:] + q[2:






    



    
