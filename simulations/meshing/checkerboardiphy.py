"""
Generate a checkerboard tomography model from an existing tomography .xyz file
For use with Specfem3D Cartesian.
"""
import os
import numpy as np
from scipy import signal
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
                     plot_fid=None):
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
    :type plot_fid: bool
    :param plot_fid: plot the overlay for confirmation
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

    # Apply the checker overlay, only to Vp and Vs values, not to rho or Q
    data_out = np.copy(data)
    data_out[:, 3] = data_out[:, 3] + (perturbation * data_out[:, 3])  # Vs
    data_out[:, 4] = data_out[:, 4] + (perturbation * data_out[:, 4])  # Vp

    # Generate a quick plot to show a representation of the checkerboard
    if plot_fid is not None:
        z_ind = np.where(data[:, 2] > 0)[0]
        plt.scatter(x=data[z_ind, 0], y=data[z_ind, 1],
                    c=checker_overlay[z_ind]
                    )
        plt.xlabel("UTM-60 EASTING (m)")
        plt.ylabel("UTM-60 NORTHING (m)")
        plt.title("{f}\n +/- {p}, {t} taper, {s}m spacing".format(
            f=xyz_fid, p=perturbation, t=taper_signal.__name__, s=spacing_m)
        )
        # plot coastline
        coastline = np.load("./nz_resf_coast_mod_utm60H_xyz.npy")
        coastline = coastline[
                np.where((coastline[:, 0] > header["orig_x"]) &
                         (coastline[:, 0] < header["end_x"]) & 
                         (coastline[:, 1] > header["orig_y"]) &
                         (coastline[:, 1] < header["end_y"])
                         )[0]]
        plt.scatter(coastline[:, 0], coastline[:, 1], c='k', marker='.')
        plt.savefig("{}.png".format(plot_fid))

    return checker_overlay, data_out


def parse_data_to_header(data):
    """
    Make a new header dictionary with adjusted data
    :param data:
    :return:
    """
    x_values = np.unique(data[:, 0])
    y_values = np.unique(data[:, 1])
    z_values = np.unique(data[:, 2])
    parsed_header = {
        "orig_x": x_values.min(), "orig_y": y_values.min(),
        "orig_z": z_values.min(), "end_x": x_values.max(),
        "end_y": y_values.max(), "end_z": z_values.max(),
        "spacing_x": abs(x_values[1] - x_values[0]),
        "spacing_y": abs(y_values[1] - y_values[0]),
        "spacing_z": abs(z_values[1] - z_values[0]),
        "nx": len(x_values), "ny": len(y_values), "nz": len(z_values),
        "vp_min": data[:, 4].min(), "vp_max": data[:, 3].max(),
        "vs_min": data[:, 4].min(), "vs_max": data[:, 4].max(),
        "rho_min": data[:, 5].max(), "rho_max": data[:, 5].max(),
    }
    return parsed_header


def write_xyz(header, data, fid_out):
    """
    Write out a new xyz file with proper header and data
    """
    with open(fid_out, "w") as f:
        # write header
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
            header["orig_x"], header["orig_y"], header["orig_z"],
            header["end_x"], header["end_y"], header["end_z"])
        )
        f.write("{:.1f} {:.1f} {:.1f}\n".format(
            header["spacing_x"], header["spacing_y"], header["spacing_z"])
                )
        f.write("{:d} {:d} {:d}\n".format(
            header["nx"], header["ny"], header["nz"])
        )
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
            header["vp_min"], header["vp_max"], header["vs_min"],
            header["vs_max"], header["rho_min"], header["rho_max"])
        )
        # write data by line
        for line in data:
            x, y, z, vp, vs, rho, qp, qs = line.tolist()
            f.write(
             "{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
                   x, y, z, vp, vs, rho, qp, qs)
            )


def call_checkerboardiphy():
    """
    Call script for the checkerboard function
    :return:
    """
    path = "./"
    fid_template = "nz_north_eberhart2015_{}.xyz"
    spacing = 80000.
    chosen_signal = signal.hanning

    # Create checkers with varying levels of perturbation
    for perturbation in [0, 0.02, 0.05, 0.1]:
        for section in ["shallow", "crust", "mantle"]:
            fid = fid_template.format(section)

            # Save the outputs with a new tag
            tag = "_checker_{space}km_{win}_{prt}pct.".format(
                space=int(spacing * 1E-3), win=chosen_signal.__name__,
                prt=int(perturbation * 1E2)
            )
            fid_out = (fid.split('.')[0] + tag + fid.split('.')[1])

            # Only plot the shallow section
            if section == "shallow":
                plot_fid = fid_out
            else:
                plot_fid = None

            # Create the checkerboard data
            fid = fid_template.format(section)
            overlay, checkerboard_data = checkerboardiphy(
                xyz_fid=os.path.join(path, fid), spacing_m=spacing,
                perturbation=perturbation, taper_signal=chosen_signal, plot_fid=plot_fid
            )
            checkerboard_header = parse_data_to_header(checkerboard_data)

            # Save the new data to a numpy array for easy i/o
            np.save(fid_out, checkerboard_data)

            # Write the new data to a new xyz file
            write_xyz(header=checkerboard_header, data=checkerboard_data,
                      fid_out=fid_out)


if __name__ == "__main__":
    call_checkerboardiphy()







