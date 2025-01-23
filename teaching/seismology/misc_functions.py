"""
Simple things I code up that may be useful to keep
"""
import numpy as np

def reflection_coefficient(rho_1, v_1, rho_2, v_2):
    """
    SW2.2.15 1D wave of a string ref. coeff. for
    density (rho) and velocity (v)
    """
    num = rho_1 * v_1 - rho_2 * v_2
    den = rho_1 * v_1 + rho_2 * v_2
    return num/den


def transmission_coefficient(rho_1, v_1, rho_2, v_2):
    """
    SW2.2.16 1D wave of a string tran. coeff. for
    density (rho) and velocity (v)
    """
    num = 2 * rho_1 * v_1
    den = rho_1 * v_1 + rho_2 * v_2
    return num/den

