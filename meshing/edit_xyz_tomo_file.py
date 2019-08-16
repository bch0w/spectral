import os
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


def lonlat_utm(lon_or_x, lat_or_y, utm_zone=60, inverse=False):
    """convert latitude and longitude coordinates to UTM projection
    from mesh_gen_helper.py (personal code)

    :type lon_or_x: float or int
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
    :type lat_or_y: float or int
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
    :type utm_zone: int
    :param utm_zone: UTM zone for conversion from WGS84
    :type inverse: bool
    :param inverse: if inverse == False, latlon => UTM, vice versa.
    :rtype x_or_lon: float
    :return x_or_lon: x coordinate in UTM or longitude in WGS84
    :rtype y_or_lat: float
    :return y_or_lat: y coordinate in UTM or latitude in WGS84
    """
    from pyproj import Proj
    projstr = (f"+proj=utm +zone={utm_zone}, +south +ellps=WGS84 +datum=WGS84"
               " +units=m +no_defs")
    myProj = Proj(projstr)
    x_or_lon, y_or_lat = myProj(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat


def extrapolate_outward(xyz_npy_fid, extend_x=0, extend_y=500000):
    """
    Extrapolate the boundaries of a mesh for each depth slice, and for given 
    amounts in the x and y hat directions

    :type xyz_npy_fid: str
    :param xyz_npy_fid: numpy .npy file containing the tomo file data
    :type extend_x: int
    :param extend_x: distance (meters) to extend x in either direction
    :type extend_y: int
    :param extend_y: distance (meters) to extend y in either direction
    :type xyz_new: np.array
    :return: the newly extrapolated array
    """
    # Load and distribute the data
    xyz = np.load(xyz_npy_fid)
    x_values = np.unique(xyz[:, 0])
    y_values = np.unique(xyz[:, 1])
    z_values = np.unique(xyz[:, 2])

    # Set some useful values to be used during extrapolation
    dx = x_values[1] - x_values[0]
    dy = y_values[1] - y_values[0]
    x_min = x_values.min()
    x_max = x_values.max()
    y_min = y_values.min()
    y_max = y_values.max()

    # Determine how much to extend the file, based on the original grid spacing,
    # if the desired extension is less than one grid space, just set as grid
    if extend_x:
        extend_x = int(extend_x//dx) or 1
    if extend_y:
        extend_y = int(extend_y//dy) or 1

    # Extrapolate for each depth slice
    xyz_new = np.array([])
    for z in z_values:
        depth_slice = xyz[np.where(xyz[:, 2] == z)[0]]
        # Extend the tomo file on the x-axis
        if extend_x:
            # Extend the tomo file on the x max side
            depth_slice_last_column = depth_slice[
                np.where(depth_slice[:, 0] == x_max)[0]]
            for i in range(1, extend_x + 1):
                depth_slice_add = depth_slice_last_column.copy()
                depth_slice_add[:, 0] += i * dx
                depth_slice = np.concatenate((depth_slice, depth_slice_add))

            # Extend the tomo file on the x min side
            depth_slice_first_column = depth_slice[
                np.where(depth_slice[:, 0] == x_min)[0]]
            for i in range(-extend_x, 0):
                depth_slice_add = depth_slice_first_column.copy()
                depth_slice_add[:, 0] += i * dx
                depth_slice = np.concatenate((depth_slice, depth_slice_add))

        # Do the same for the y extent
        if extend_y:
            # Extend the tomo file on the y max side
            depth_slice_top_row = depth_slice[
                np.where(depth_slice[:, 1] == y_max)[0]]
            for i in range(1, extend_y + 1):
                depth_slice_add = depth_slice_top_row.copy()
                depth_slice_add[:, 1] += i * dy
                depth_slice = np.concatenate((depth_slice, depth_slice_add))

            # Extend the tomo file on the y min side
            depth_slice_bottom_row = depth_slice[
                np.where(depth_slice[:, 1] == y_min)[0]]
            for i in range(-extend_y, 0):
                depth_slice_add = depth_slice_bottom_row.copy()
                depth_slice_add[:, 1] += i * dy
                depth_slice = np.concatenate((depth_slice, depth_slice_add))

        # Add the depth slices to the output array
        if len(xyz_new) == 0:
            xyz_new = depth_slice
        else:
            xyz_new = np.concatenate((xyz_new, depth_slice))

    return xyz_new


def convert_interp(latlon_fid, data_fid, spacing, plot=False, tag=''):
    """
    Matlab spits out some converted lat lon coordinates but they need to be in 
    UTM60 and also need to be interpolated because sampling in Carl's
    custom coordinate system is not the same as sampling in lat/lon or utm

    :type latlon_fid: str
    :param latlon_fid: the output.txt file given by matlab
    :type data_fid: str
    :param data_fid: path to the extrapolated data, saved as a .npy file
    :type spacing: float
    :param spacing: the desired grid spacing in the x and y directions
    :type plot: bool
    :param plot: plot a slice to check things are working
    :type tag: str
    :param tag: tag of the file, e.g. mantle
    """
    # Load in the all the data and distribute
    data = np.load(data_fid)
    latlon = np.loadtxt(latlon_fid)
    lat = latlon[:, 0]
    lon = latlon[:, 1]
    x_irreg, y_irreg = lonlat_utm(lon, lat)
    
    # These are taken from converting my chosen map coordinates to UTM-60 and
    # putting a small buffer around the dimensions and making sure they are 
    # in line with the nz_north tomography file
    # lat_min, lon_min = 172.9998, -42.5007;
    # lat_max, lon_max = 179.5077, -36.9488
    x_min, x_max = 170000., 725000.
    y_min, y_max = 5271000., 5910000.

    # For the coastline plot
    if plot:
        coastline = np.load("./nz_resf_coast_mod_utm60H_xyz.npy")
        coastline = coastline[
            np.where((coastline[:, 0] > x_min) &
                     (coastline[:, 0] < x_max) &
                     (coastline[:, 1] > y_min) &
                     (coastline[:, 1] < y_max)
                     )[0]]

    # Create the regular grid to be interpolated onto
    x_reg = np.arange(x_min, x_max, spacing)
    y_reg = np.arange(y_min, y_max, spacing)
    x_grid, y_grid = np.meshgrid(x_reg, y_reg)
    x_out = x_grid.flatten()
    y_out = y_grid.flatten()
    
    # Interpolate for each depth slice
    list_interp = []
    data_out = None
    for z_value in np.unique(data[:, 2]):
        print(z_value)
        data_ind = np.where(data[:, 2] == z_value)[0]

        # Loop over the values to be interpolated
        for val in range(3, 8):
            interp_vals = griddata(
                points=(x_irreg[data_ind], y_irreg[data_ind]),
                values=data[data_ind, val],
                xi=(x_grid, y_grid)
            )
            list_interp.append(interp_vals.flatten())
            
            # If plot, only plot one depth slice
            if plot:
                plt.contourf(x_grid, y_grid, interp_vals)
                plt.scatter(coastline[:, 0], coastline[:, 1], c='k',
                            marker='.', s=1.5
                            )
                plt.title(f"{tag} {z_value}")
                plt.xlabel("Easting (m)")
                plt.ylabel("Northing (m)")
                plt.show()
                plot = False

        # Create the ndarray by creating column vectors and mushing em together
        z_out = np.ones(len(x_grid.flatten())) * z_value
        vp_out = list_interp[0].flatten()
        vs_out = list_interp[1].flatten()
        rho_out = list_interp[2].flatten()
        qp_out = list_interp[3].flatten()
        qs_out = list_interp[4].flatten()

        # Arrange the data in the desired format
        depth_slice = np.column_stack(
            (x_out, y_out, z_out, vp_out, vs_out, rho_out, qp_out, qs_out)
        )
        # Add the data to all the other data
        if data_out is None:
            data_out = depth_slice
        data_out = np.concatenate((data_out, depth_slice))

    return data_out


def parse_data_to_header(data):
    """
    Parse data to header to write into .xyz file
    """
    parsed_header = {"orig_x": data[:, 0].min(), "orig_y": data[:, 1].min(),
                     "orig_z": data[:, 2].min(), "end_x": data[:, 0].max(),
                     "end_y": data[:, 1].max(), "end_z": data[:, 2].max(),
                     "spacing_x": 2000., "spacing_y": 2000.,
                     "spacing_z": 1000.,
                     "nx": len(np.unique(data[:, 0])),
                     "ny": len(np.unique(data[:, 1])),
                     "nz": len(np.unique(data[:, 2])),
                     "vp_min": data[:, 3].min(), "vp_max": data[:, 3].max(),
                     "vs_min": data[:, 4].min(), "vs_max": data[:, 4].max(),
                     "rho_min": data[:, 5].max(),
                     "rho_max": data[:, 5].max(),
                     }

    head_len = parsed_header["nx"] * parsed_header["ny"] * parsed_header["nz"] 
    assert head_len == len(data)

    return parsed_header


def trim_xyz_file(data, bounds):
    """
    Sometime the xyz file is too large, trim it down to new dimensions

    :type data: np.array
    :param data: data to trim
    :type bounds: list of floats
    :param bounds: bounds to trim to
    :rtype parsed_header: dict
    :return parsed_header: values of the newly trimmed xyz file
    :rtype data: np.array
    :return data: trimmed data
    """
    # convert data into numpy array for easier working
    for i, bound in enumerate(bounds):
        too_small = np.where(data[:, i] < bound[0])[0]
        too_large = np.where(data[:, i] > bound[1])[0]
    
    to_remove = np.unique(np.concatenate((too_small, too_large), 0))
    
    # delete in place
    data = np.delete(data, to_remove, 0)
    
    # make new header
    header = parse_data_to_header(data)

    return header, data


def write_new_xyz(data, fidout, write_header=False):
    """
    Write out a new xyz file with proper header and data

    :type data: np.array
    :param data: data to write into xyz file
    :type fidout: str
    :param fidout: file name to save to
    :type header: dict 
    :param header: header information to prepend to data, if not given, no head
    """
    with open(fidout, "w") as f:
        if write_header:
            # write header
            header = parse_data_to_header(data)
            
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


def plot_slice(array, slice_index=0):
    """
    Simple slice plot just to check, I was tired of rewriting this all the time

    :type array: np.array
    :param array: array to plot
    :type slice_index: int
    :param slice_index: of the given unique slices, which one?
    """
    depth_slice = array[np.where(array[:,2] == np.unique(
                                                    array[:,2])[slice_index])
                        ]
    plt.scatter(depth_slice[:,0], depth_slice[:,1], c=depth_slice[:,3])
    plt.show()

if __name__ == "__main__":
    tag = "shallow"

    if tag == "mantle":
        spacing = 8000.
    else:
        spacing = 2000.

    extrapo_fidout =  f'./nz_full_extrapolate_{tag}.xyz.npy'
    if not os.path.exists(extrapo_fidout): 
        # Extrapolate the nz_full file so that when it's rotated, no nans
        xyz_extrapolated = extrapolate_outward(
                xyz_npy_fid=f'./nz_full_eberhart2015_{tag}.xyz.npy'
                )
        
        # Save the new extrapolated file
        np.save(extrapo_fidout, xyz_extrapolated)
        
        # Plot a slice just for visual confirmation
        plot_slice(xyz_extrapolated)

        # Write the new data to be opened by matlab
        write_new_xyz(data=xyz_extrapolated, 
                      fidout=f'./nz_full_extrapolate_{tag}_noheader.xyz')
        
    # In between these two steps, matlab needs to be run
    else:
        # Convert coordinates from lat lon and interpolate along a regular grid
        data_out = convert_interp(latlon_fid=f'{tag}_latlon.txt',
                                  data_fid=extrapo_fidout, spacing=spacing, 
                                  plot=True)
       
        # Save the new xyz file
        np.save(f'./nz_full_utm60_extrapolate_{tag}.xyz.npy', data_out)



