"""
Simple script to difference VTK files for the same model.
There is no real checking so user needs to be sure that the grid data is same!
"""
import os
import sys


def read_file(pathname):
    """
    read the file and return header info and data
    :type pathname: str
    :param pathname: full path to the .vtk file to read
    :rtype lines: list
    :return lines: data line by line from readlines()
    :rtype header_dict: dic
    :return header_dict: dictionary with all relevant header information
        from an unstructured_grid vtk file
    """
    with open(pathname, "r") as f:
        lines = f.readlines()

    # determine important line numbers, headers are specific to specfem3d output
    for i, line in enumerate(lines):
        if "POINTS" in line:
            points_n = int(line.split()[1])
            points_line = i + 1
        elif "CELL" in line and not "CELL_TYPE" in line:
            cells_n = int(line.strip().split()[1])
            cells_size = int(line.strip().split()[2])
            cells_line = i
        elif "CELL_TYPES" in line:
            cell_types_n = int(line.strip().split()[1])
            cell_types_line = i
        elif "POINT_DATA" in line:
            point_data_n = int(line.split()[1])
            point_data_line = i
        elif "SCALARS" in line:
            scalars = line.strip().split()[1]
            data_line = i+2

    # easier returns in a dictionary
    header_dict = {"points_n": points_n, "points_line": points_line,
                   "cells_n": cells_n, "cells_size": cells_size,
                   "cells_line": cells_line, "cell_types_n": cell_types_n,
                   "cell_types_line": cell_types_line,
                   "points_data_n": point_data_n,
                   "points_data_line": point_data_line, "scalars": scalars,
                   "data_line": data_line
                   }

    return lines, header_dict


def difference_vtk(model_a, model_b, path="./", reverse=1, write=None):
    """
    read each model and scan line by line, difference all necessary values

    the difference is defined as: c = a - b
    to reverse, set reverse == -1

    :type model_?: str
    :param model_?: name of the model file
    :type path: str
    :param path: path holding the models to be read in
    :type reverse: int
    :param reverse: 1 for a-b, -1 for b-a
    :type write: str
    :param write: name of the output file to be written
    :rtype differences: list
    :return differences: list of the differences in values between a and b
    """
    # read files
    model_a, header_dict_a = read_file(os.path.join(path, model_a))
    model_b, header_dict_b = read_file(os.path.join(path, model_b))

    # check that the files have the same characteristics before parsing
    for key in header_dict_a.keys():
        if header_dict_a[key] != header_dict_b[key]:
            sys.exit("{} not equal".format(key))

    # parse through models together and separate by len of line, skip header
    error_count = 1
    differences = []
    start = header_dict_a["data_line"]
    for x, y in zip(model_a[start:-1], model_b[start:-1]):
        try:
            difference = reverse * (float(x.strip()) - float(y.strip()))
            differences.append(difference)
        except ValueError:
            print("value error")

    if write:
        with open(os.path.join(path, write), "w") as f:
            f.writelines(model_a[:header_dict_a["data_line"]])
            for diff in differences:
                if diff == 0:
                    f.write("{:13.5f}    \n".format(float(diff)))
                elif abs(diff) > 1:
                    f.write("{:13.5f}    \n".format(float(diff)))
                else:
                    f.write("{:13.10f}    \n".format(float(diff)))
            print("\n")

    return differences


if __name__ == "__main__":
    path = "/Users/chowbr/Documents/subduction/tomo/testing/" \
           "seisflows/vtk_files//hikurangi_trial"
    model_b = "vp_twoevent_0001.vtk"
    model_a = "vp_nz_init.vtk"

    model_out = "{}_diff_{}_{}.vtk".format(model_a.split("_")[0],
                                           model_a.split("_")[2].split(".")[0],
                                           model_b.split("_")[2].split(".")[0]
                                           )

    differences = difference_vtk(model_a, model_b, path, write=model_out)




