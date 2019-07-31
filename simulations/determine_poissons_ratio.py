import os
import glob
import numpy as np
from pyatoa.utils.operations.conversions import read_fortran_binary


def read_slices(path, qty):
    """
    Read all individual proc files for a given path and quantity
    :type path: str
    :param path: path to bin files
    :type qty: str
    :param qty: quantity to read, e.g. 'vp', 'vs', 'x' etc.
    :rtype: np.array
    :return: data array
    """
    data = np.array([])
    binfiles = glob.glob(os.path.join(path, "*{}.bin".format(qty)))
    binfiles.sort()
    print("reading {} to {}".format(os.path.basename(binfiles[0]),
                                    os.path.basename(binfiles[-1]))
          )
    for proc in binfiles:
        temp = read_fortran_binary(proc)
        data = np.concatenate((data, temp))

    return data


def poissons_ratio(vp, vs):
    """
    Calculate Poisson's ratio based on Specfem3d Cartesian
    shared/check_mesh_resolution.f90

    Theoretical limits of Piisson's ratio -1 < sigma < 0.5



    :type vp: np.array or float:q
    :param vp: p-wave velocity
    :type vs: np.array or float
    :param vs: s-wave velocity
    :rtype: float
    :return: Poissons Ratio
    """
    return 0.5 * (vp ** 2 - (2 * vs ** 2)) / ((vp ** 2) - (vs ** 2))


path_to_xyz = "./xyz"
path_to_vel = "./model_0009"

# Read in Vp values
vp = read_slices(path_to_vel, "vp")
vs = read_slices(path_to_vel, "vs")
# x = read_slices(path_to_xyz, "x")
# y = read_slices(path_to_xyz, "y")
# z = read_slices(path_to_xyz, "z")

poissons = poissons_ratio(vp, vs)

