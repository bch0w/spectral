"""
Generate a checkerboard tomography model from an existing tomography .xyz file
For use with Specfem3D Cartesian.
"""
import os
import numpy as np
from scipy import signal
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rcParams['image.cmap'] = 'seismic'


def xyz_reader(xyz_fid, save=True):
    """
    read .xyz file used as input for specfem and return information
    header parsing specifically related to specfem
    data is delineated as :
    x, y, z, vp[m/s], vs[m/s], rho[kg/m**3], Qp, Qs
    -origin and endpoints are the coordinates in meters
    -spacing given in units of meters for UTM projections
    -nx, ny, nz = [(end_x - orig_x)/spacing_x] + 1
    :type xyz_fid: str
    :param xyz_fid: file to read in
    :type save: bool
    :param save: to save the read files into numpy arrays
    :rtype header: dict
    :return header: a dictionary with relevant header informadtion
    :rtype data: np.ndarray
    :return data: all the data contained in the xyz file
    """
    if os.path.exists(xyz_fid + ".npy") and os.path.exists(xyz_fid + ".npz"):
        header = np.load(xyz_fid + ".npz")
        data = np.load(xyz_fid + ".npy")
        return header, data
    else:
        with open(xyz_fid) as f:
            lines = f.readlines()

            # parse the header and make sure the values are floats
            orig_x, orig_y, orig_z, end_x, end_y, end_z = \
                [float(_) for _ in lines[0].strip().split()]
            spacing_x, spacing_y, spacing_z = \
                [float(_) for _ in lines[1].strip().split()]
            nx, ny, nz = [float(_) for _ in lines[2].strip().split()]
            vp_min, vp_max, vs_min, vs_max, rho_min, rho_max = \
                [float(_) for _ in lines[3].strip().split()]

            header = {"orig_x": orig_x, "orig_y": orig_y, "orig_z": orig_z,
                      "end_x": end_x, "end_y": end_y, "end_z": end_z,
                      "spacing_x": spacing_x, "spacing_y": spacing_y,
                      "spacing_z": spacing_z, "nx": nx, "ny": ny, "nz": nz,
                      "vp_min": vp_min, "vp_max": vp_max, "vs_min": vs_min,
                      "vs_max": vs_max, "rho_min": rho_min, "rho_max": rho_max,
                      }

            data = np.array([_.strip().split() for _ in lines[4:]], dtype=float)
            if save:
                np.save(xyz_fid, data)  # save data in .npy
                np.savez(xyz_fid, **header)  # save header in .npz

        return header, data


def checkerboardiphy(xyz_fid, spacing_m, perturbation=0.02, taper_signal=None,
                     plot=True):
    """
    Read in the data from an XYZ tomography file and create a checkerboard
    overlay which has +/- {perturbation} checkers overlain on the tomography
    file. The checkers also include tapering so that there are no sharp
    transitions inside the checkerboard

    :type xyz_fid: str
    :param xyz_fid: path to file to read
    :type spacing_m: float
    :param spacing_m: spacing in meters
    :type perturbation: float
    :param perturbation: perturbation to give to checkers
    :type taper_signal: scipy.signal.function
    :param taper_signal: scipy signal to taper checkers by
    :type plot: bool
    :param plot: plot the overlay for confirmation
    :return:
    """
    # Read in the data
    header, data = xyz_reader(xyz_fid=xyz_fid, save=True)

    # Initialize starting values
    checker_overlay = np.zeros(len(data))
    x = 1

    print(xyz_fid)
    # Loop through the x-axis, setting left and right boundaries
    for x_left in np.arange(header["orig_x"], header["end_x"], spacing_m):
        y = 1
        print("x_left: {}/{}".format(x_left, header["end_x"]))
        x_right = x_left + spacing_m

        # Create the tapered checker for a given checker defined by bounds
        x_checker = np.arange(x_left, x_right, header["spacing_x"])
        x_window = taper_signal(len(x_checker))

        # Loop through the y-axis, setting lower and upper boundaries
        for y_bot in np.arange(header["orig_y"], header["end_y"], spacing_m):
            print("\ty_bot: {}/{}".format(y_bot, header["end_y"]))
            y_top = y_bot + spacing_m

            # Create the tapered checker for the given checker defined by bounds
            y_checker = np.arange(y_bot, y_top, header["spacing_y"])
            y_window = taper_signal(len(y_checker))

            # Determine the sign of the overall checker, which is set by the
            # alternating sign of each inner checker axis
            checker = x * y

            # Determine where the data falls within this checker's bounds
            checker_indices = np.where(
                (data[:, 0] >= x_left) & (data[:, 0] < x_right) &
                (data[:, 1] >= y_bot) & (data[:, 1] < y_top)
            )
            # For each of the given indices, figure out the resulting overlay
            # Which is a multiplication of the overall sign, and the combination
            # of the x and y tapers within the checker
            for ind in checker_indices[0]:
                checker_overlay[ind] = checker * (
                        x_window[np.where(x_checker == data[ind, 0])[0]] *
                        y_window[np.where(y_checker == data[ind, 1])[0]]
                )

            # Flip the sign of the y-axis checker
            y *= -1
        # Flip the sign of the x-axis checker
        x *= -1

    # Apply the checker overlay
    data_out = np.copy(data)
    for i in range(3, np.shape(data)[1]):
        data_out[:, i] = data_out[:, i] + (perturbation * data_out[:, i])

    # Generate a quick plot to show a representation of the checkerboard
    if plot is not None:
        z_ind = np.where(data[:, 2] > 0)[0]
        plt.scatter(x=data[z_ind, 0], y=data[z_ind, 1],
                    c=checker_overlay[z_ind]
                    )
        plt.xlabel("UTM-60 EASTING (m)")
        plt.ylabel("UTM-60 NORTHING (m)")
        plt.title("{f}\n +/- {p}, {t} taper, {s}m spacing".format(
            f=xyz_fid, p=perturbation, t=chosen_signal.__name__, s=spacing_m)
        )
        plt.savefig("{}_checker_{}km.png".format(xyz_fid, int(spacing*1E-3)))

    return checker_overlay, data_out


if __name__ == "__main__":
    path = "./"
    fid_template = "nz_north_eberhart2015_{}.xyz"
    spacing = 80000.
    chosen_signal = signal.hanning

    for section in ["shallow", "crust", "mantle"]:
        fid = fid_template.format(section)
        overlay, checkerboard = checkerboardiphy(
            xyz_fid=os.path.join(path, fid), spacing_m=spacing.,
            perturbation=0.02, taper_signal=chosen_signal, plot=True
        )
        # save the data with a new tag
        fid_out = (fid.split('.')[0] +
                   "_checker_{}km.".format(int(spacing*1E-3)) +
                   fid.split('.')[1]
                   )
        np.save(fid_out, checkerboard)








