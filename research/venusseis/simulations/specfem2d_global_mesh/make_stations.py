"""
Make a STATIONS file for stations on the surface at various angles 
"""
import numpy as np

def surface_xy(r, angle):
    x = r * 1E3 * np.cos(np.radians(angle))
    y = r * 1E3 * np.sin(np.radians(angle))
    return x, y

r_venus = 6051.8

for angle in [90, 75, 60, 45, 30, 15, 0, -15, -30, -45, -60, -75, -90]:
    x, y = surface_xy(r=r_venus, angle=angle)
    if angle < 0:
        name = str(angle).replace("-", "N")
    else:
        name = f"{angle:0>3}"
    print(f"S{name}\tXX\t{x:10.2f}\t{y:10.2f}\t0.0\t0.0")


