"""
Generate a checkerboard tomography model from an existing tomography .xyz file
For use with Specfem3D Cartesian.
"""
import os
import numpy as np
from scipy import signal


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


def determine_checkers(data, header, spacing_m=50000):
    """
    read files in, define bounds, return
    :type data: np.ndarray
    :param data: data read in using xyz_reader
    :type header: dict
    :param header: the header from xyz_reader
    :type spacing_m: int
    :param spacing_m: spacing of checkers in meters
    :rtype checkerboard: np.array
    :return checkerboard: the checkerboard the same size as the input data
    """
    def define_checkers(point, min_bound, max_bound, checker_size, grid_size):
        """
        1-D operation
        Given a point along an axis, first determine whether this point exists
        in a positive or negative checker, and then determine where along a
        tapering inside this checker the point exists.

        :type point: float
        :param point: coordinate point to be evaluated
        :type min_bound: float
        :param min_bound: minimum domain boundary
        :type max_bound: float
        :param max_bound: maximum domain boundary
        :type checker_size: float
        :param checker_size: the spacing of the overall domain
        :type grid_size: float
        :param grid_size: the grid spacing along the given axis
        :rtype: int
        :return: value of the checker for the given information
        """
        sign = 1
        # check if the value is between
        for inner_bound in np.arange(min_bound, max_bound, checker_size):
            outer_bound = inner_bound + checker_size
            if inner_bound <= point < outer_bound:
                # Define the extent of the checker grid
                checker = np.arange(inner_bound, outer_bound, grid_size)

                # Taper the checker (NOTE: Different signals can be used here)
                taper_window = signal.tukey(len(checker))
                taper_amount = taper_window[np.where(point == checker)[0][0]]

                # Assign value based on the sign of the checker and tapering
                checker_value = sign * taper_amount

                return checker_value
            # If the value is not within these bounds, switch sign and move on
            sign *= -1
        else:
            # Takes the last value of the inner and outer bound
            print("this shouldn't have happened")
            return

    # Go through each line in the file and determine the checker value
    x_grid, y_grid, checkerboard_out = [], [], []
    for gridpoint in data:
        x_checker = define_checkers(point=gridpoint[0],  # this is the x value
                                    min_bound=header["orig_x"],
                                    max_bound=header["end_x"],
                                    checker_size=spacing_m,
                                    grid_size=header["spacing_x"]
                                    )
        y_checker = define_checkers(point=gridpoint[1],  # this is the y value
                                    min_bound=header["orig_y"],
                                    max_bound=header["end_y"],
                                    checker_size=spacing_m,
                                    grid_size=header["spacing_y"]
                                    )
        checkerboard_out.append(x_checker * y_checker)

    return np.array(checkerboard_out)


def convert_to_checkers(data, checkerboard, perturbation=0.1):
    """
    take checker boolean data and perturb base model with given perturbation

    x, y, z, vp[m/s], vs[m/s], rho[kg/m**3], Qp, Qs

    :return:
    """
    data_out = np.copy(data)
    for i, gridpoint in enumerate(data):
        # Set the value of the perturbation to add to the gridpoint
        perturb = checkerboard[i] * perturbation
        for j, value in gridpoint:
            # Make sure we skip over the coordinates, which occupy first 3
            if j in np.arange(0, 3):
                continue
            # New data is old data plus some perturbation of old data
            data_out[i] = data_in[i] + data_in[i] * perturb

    return data_out


# def checkerboardiphy(xyz_fid, spacing_m, perturbation=0.02):
#     """
#     Amalgamation of the above functions
#
#     This is too slow, doesn't work
#
#     :type xyz_fid: str
#     :param xyz_fid: path to file to read
#     :type spacing_m: float
#     :param spacing_m: spacing in meters
#     :type perturbation: float
#     :param perturbation: perturbation to give to checkers
#     :return:
#     """
#     # Read in the data
#     header, data = xyz_reader(xyz_fid=xyz_fid, save=True)
#
#     # Create the checkerboard based on the given data
#     checkerboard = determine_checkers(data, header, spacing_m)
#
#     # Convert the data into checkerboard
#     data_out = convert_to_checkers(data, checkerboard, perturbation)
#
#     return data_out


def checkerboardiphy(xyz_fid, spacing_m, perturbation=0.02):
    """
    Amalgamation of the above functions

    :type xyz_fid: str
    :param xyz_fid: path to file to read
    :type spacing_m: float
    :param spacing_m: spacing in meters
    :type perturbation: float
    :param perturbation: perturbation to give to checkers
    :return:
    """
    # Read in the data
    header, data = xyz_reader(xyz_fid=xyz_fid, save=True)

    # Initialize starting values
    checker_overlay = np.zeros(len(data))
    x, y = 1, 1

    # Loop through the x-axis, setting left and right boundaries
    for x_left in np.arange(header["orig_x"], header["end_x"], spacing_m):
        x_right = x_left + spacing_m

        # Create the tapered checker for a given checker defined by bounds
        x_checker = np.arange(x_left, x_right, header["spacing_x"])
        x_window = signal.tukey(len(x_checker))

        # Loop through the y-axis, setting lower and upper boundaries
        for y_bot in np.arange(header["orig_y"], header["end_y"], spacing_m):
            y_top = y_bot + spacing_m

            # Create the tapered checker for the given checker defined by bounds
            y_checker = np.arange(y_bot, y_top, header["spacing_y"])
            y_window = signal.tukey(len(y_checker))

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

    return checker_overlay, data_out


def plot_overlay():
    """
    Quick visualizations of the checkerboard overlay to make sure things are ok
    :return:
    """


if __name__ == "__main__":
    path = "./"
    fid_template = "nz_north_eberhart2015_{}.xyz"
    for section in ["shallow", "crust", "header"]:
        fid = fid_template.format(section)
        overlay, checkerboard = checkerboardiphy(xyz_fid=os.path.join(path, fid),
                                                 spacing_m=50000.,
                                                 perturbation=0.02
                                                 )
        # save the data with a new tag
        fid_out = fid.split('.')[0] + "_checker." + fid.split('.')[1]
        np.save(checkerboard, fid_out)








