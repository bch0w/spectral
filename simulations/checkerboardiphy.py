"""
Generate a checkerboard tomography model from an existing tomography .xyz file
For use with Specfem3D Cartesian
"""
import os
import numpy as np


def xyz_reader(filepath):
    """
    read .xyz file used as input for specfem and return information
    header parsing specifically related to specfem
    data is delineated as :
    x, y, z, vp[m/s], vs[m/s], rho[kg/m**3], Qp, Qs
    -origin and endpoints are the coordinates in meters
    -spacing given in units of meters for UTM projections
    -nx, ny, nz = [(end_x - orig_x)/spacing_x] + 1
    :param filepath:
    :return:
    """
    with open(filepath) as f:
        lines = f.readlines()

    # parse the header and make sure the values are floats
    orig_x, orig_y, orig_z, end_x, end_y, end_z = \
        [float(_) for _ in lines[0].strip().split()]
    spacing_x, spacing_y, spacing_z = \
        [float(_) for _ in lines[1].strip().split()]
    nx, ny, nz = [float(_) for _ in lines[2].strip().split()]
    vp_min, vp_max, vs_min, vs_max, rho_min, rho_max = \
        [float(_) for _ in lines[3].strip().split()]

    parsed_header = {"orig_x": orig_x, "orig_y": orig_y, "orig_z": orig_z,
                     "end_x": end_x, "end_y": end_y, "end_z": end_z,
                     "spacing_x": spacing_x, "spacing_y": spacing_y,
                     "spacing_z": spacing_z, "nx": nx, "ny": ny, "nz": nz,
                     "vp_min": vp_min, "vp_max": vp_max, "vs_min": vs_min,
                     "vs_max": vs_max, "rho_min": rho_min, "rho_max": rho_max,
                      }
    return lines[:4], lines[4:], parsed_header




def determine_checkers(data_list, bounds, spacing_m=50000, include_depth=False):
    """
    read files in, define bounds, return
    :type data_list: list
    :param data_list: data read in using xyz_reader
    :type bounds: dict
    :param bounds: the bounds for each coordinate, e.g. {"x": [x_min, xmax]...}
    :type data_out: list of str
    :param data_out: if path should be writen
    :type spacing_m: int
    :param spacing_m: spacing of checkers in meters
    :type include_depth: bool
    :param include_depth: checker in the z/depth coordinate
    :rtype ischeckers: list
    :return ischeckers: bool telling if a line is + or - checker, 0 or 1
    """
    def define_bounds(value, var_bounds):
        """
        x and y bounds are shared by each of the layers: shallow, crust, mantle
        :type value: float
        :param value: value to be evaluated
        :type var_bounds: list of floats
        :param var_bounds: [min_bound, max_bound]
        :rtype: bool
        :return: bool to say if the value falls within a checker or not
        """
        # defining z bounds requires extra information that x and y don't need
        min_value = var_bounds[0]
        max_value = var_bounds[1]

        checker_bool = True
        # check if the value is between
        for bound in np.arange(min_value, max_value, spacing_m):
            if bound <= value < bound + spacing_m:
                return checker_bool
            checker_bool = not checker_bool

        return checker_bool

    # go through each data file
    ischecker = [[] for _ in range(len(data_list))]
    for i, data in enumerate(data_list):
        for line in data:
            x_bool = define_bounds(float(line.strip().split()[0]), bounds["x"])
            y_bool = define_bounds(float(line.strip().split()[1]), bounds["y"])

            # we don't always want to include checkering to depth because
            # that affects turning of rays
            if include_depth:
                z_bool = define_bounds(
                    float(line.strip().split()[2]), bounds["z"])
                # checker only if all coordinates are true
                # OR if all components are false
                if x_bool == y_bool == z_bool:
                    checker = 1
                else:
                    checker = -1
            else:
                if x_bool == y_bool:
                    checker = 1
                else:
                    checker = -1
            ischecker[i].append(checker)

    return ischecker


def checkerboardiphy(data, ischeckers, perturbation=0.1):
    """
    take checker boolean data and perturb base model with given perturbation

    x, y, z, vp[m/s], vs[m/s], rho[kg/m**3], Qp, Qs

    :return:
    """
    for i, line in enumerate(data):
        (x, y, z, vp, vs, rho, qp, qs) = line.split()
        perturb = ischeckers[i] * perturbation
        vp = float(vp) + (float(vp) * perturb)
        vs = float(vs) + (float(vs) * perturb)
        rho = float(rho) + (float(rho) * perturb)
        qp = float(qp) + (float(qp) * perturb)
        qs = float(qs) + (float(qs) * perturb)
        data[i] = (
            "{} {} {} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
                x, y, z, vp, vs, rho, qp, qs)
        )

    return data


def write_out(header, data, file_id):
    """
    simple writer to create new tomography files
    :param header:
    :param data:
    :param file_id:
    :return:
    """
    with open(file_id, "w") as f:
        for h in header:
            f.write(h)
        for d in data:
            f.write(d)


if __name__ == "__main__":
    # read in data
    # path = "/Users/chowbr/Documents/subduction/data/KUPEDATA/tomo_files"
    path = "/seis/prj/fwi/bchow/data/KUPEDATA/tomo_files" 
    name_template = "nz_x1200_y600_eberhart2015_{}.xyz"
    fullpath = os.path.join(path, name_template)

    shallow_header, shallow_data, shallow_parsed = xyz_reader(
        fullpath.format("shallow"))
    crust_header, crust_data, crust_parsed = xyz_reader(
        fullpath.format("crust"))
    mantle_header, mantle_data, mantle_parsed = xyz_reader(
        fullpath.format("mantle"))
    data_list = [shallow_data, crust_data, mantle_data]

    # user defined paramters
    spacing_m = 50000.
    perturbation = 0.1

    # workflow
    bounds_dict = {"x": [shallow_parsed["orig_x"], shallow_parsed["end_x"]],
                   "y": [shallow_parsed["orig_y"], shallow_parsed["end_y"]],
                   "z": [mantle_parsed["orig_z"], shallow_parsed["end_z"]]
                   }
    ischeckers = determine_checkers(
        data_list, bounds_dict, spacing_m=spacing_m, include_depth=False
    )

    headers = [shallow_header, crust_header, mantle_header]
    data_out = [fullpath.format("shallow_checker"),
                fullpath.format("crust_checker"),
                fullpath.format("mantle_checker")
                ]
    for i, data in enumerate(data_list):
        new_data = checkerboardiphy(data, ischeckers[i], perturbation)
        write_out(headers[i], new_data, data_out[i])





