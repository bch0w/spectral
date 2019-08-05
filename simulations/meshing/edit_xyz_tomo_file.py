import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

from checkerboardiphy import xyz_reader


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
    projstr = ("+proj=utm +zone={}, +south +ellps=WGS84 +datum=WGS84"
               " +units=m +no_defs".format(utm_zone))
    myProj = Proj(projstr)
    x_or_lon, y_or_lat = myProj(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat


def convert_interp(latlon_fid, xyz_npy_fid, spacing, plot=False):
    """
    Matlab spat out some converted lat lon coordinates but they need to be in 
    UTM60 and also need to be interpolated because sampling in Carl's
    custom coordinate system is not the same as sampling in lat/lon or utm
    """
    latlon = np.loadtxt(latlon_fid)
    lat = latlon[:, 0]
    lon = latlon[:, 1]
    x_irreg, y_irreg = lonlat_utm(lon, lat)
    
    # Read in the data file to replace the coordinates
    data = np.load(xyz_npy_fid)

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
        data_ind = np.where(data[:, 2] == z_value)[0]

        # Loop over the values to be interpolated
        for val in range(3, 8):
            interp_vals = griddata(
                points=(x_irreg[data_ind], y_irreg[data_ind]),
                values=data[data_ind, val],
                xi=(x_grid, y_grid)
            )
            list_interp.append(interp_vals.flatten())

            if plot:
                plt.contourf(x_grid, y_grid, interp_vals)
                plt.scatter(coastline[:, 0], coastline[:, 1], c='k',
                            marker='.', s=1.5
                            )
                plt.title("{fid} {z}".format(fid=xyz_npy_fid, z=z_value))
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

        depth_slice = np.column_stack(
            (x_out, y_out, z_out, vp_out, vs_out, rho_out, qp_out, qs_out)
        )
        if data_out is None:
            data_out = depth_slice
        data_out = np.concatenate((data_out, depth_slice))

    np.save('utm60_{}'.format(xyz_npy_fid), data_out)
    return data_out


def fill_nans(data_with_nans, data_to_fill):
    """
    Combining two tomography files

    :param fid_to_fill:
    :param fid_fill_with:
    :return:
    """
    # Determine the bounds of the good data, to avoid searching outside
    x_min_dtf = data_to_fill[:,0].min()
    x_max_dtf = data_to_fill[:,0].max()
    y_min_dtf = data_to_fill[:,1].min()
    y_max_dtf = data_to_fill[:,1].max()
        
    # Find rows that contain NaN's
    nan_ind = np.where((np.isnan(data_with_nans[:, 3]) == True) &
                       (data_with_nans[:,0] > x_min_dtf) &
                       (data_with_nans[:,0] < x_max_dtf) &
                       (data_with_nans[:,1] > y_min_dtf) &
                       (data_with_nans[:,1] < y_max_dtf)
                       )[0]

    for i in nan_ind:
        print(i)
        to_fill_ind = np.where((data_to_fill[:, 0] == data_with_nans[i, 0]) &
                               (data_to_fill[:, 1] == data_with_nans[i, 1]) &
                               (data_to_fill[:, 2] == data_with_nans[i, 2])
                               )[0]

        try:
            to_fill_ind.size()
            data_with_nans[to_fill_ind] = data_with_nans[i]
        except TypeError:
            continue

    pass


def trim_xyz_file(data, bounds):
    """
    sometime the xyz file is too large, trim it down to new dimensions
    :param header:
    :param data:
    :return:
    """
    # convert data into numpy array for easier working
    for i, bound in enumerate(bounds):
        too_small = np.where(data[:, i] < bound[0])[0]
        too_large = np.where(data[:, i] > bound[1])[0]
    
    to_remove = np.unique(np.concatenate((too_small, too_large), 0))
    
    # delete in place
    data = np.delete(data, to_remove, 0)
    
    # make new header
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

    return parsed_header, data


def write_new_xyz(header, data, fidout_template):
    """
    write out a new xyz file with proper header and data
    """
    fid = fidout_template.format(
            x=(header["end_x"] - header["orig_x"])*1E-3,
            y=(header["end_y"] - header["orig_y"])*1E-3,
            )

    with open(fid, "w") as f:
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

            
def write_xyz_from_np(fid_in, fid_out):
    """
    write out a new xyz file with proper header and data
    """
    data = np.load(fid_in)
    with open(fid_out, "w") as f:
        for line in data:
            x, y, z, vp, vs, rho, qp, qs = line.tolist()
            f.write(
             "{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
                   x, y, z, vp, vs, rho, qp, qs)
            )


if __name__ == "__main__":
    data_out = convert_interp(latlon_fid='nz_full_extrapolate_crust_latlon_coordinates.txt',
                              xyz_npy_fid='nz_full_extrapolate_eberhart2015_crust.xyz.npy',
                              spacing=2000., plot=True)
    # data_out = convert_interp('shallow_latlon.txt',
    #                           'nz_full_eberhart2015_shallow.xyz.npy',
    #                           spacing=2000.
    #                           )
    # data_out = convert_interp('crust_latlon.txt',
    #                           'nz_full_eberhart2015_crust.xyz.npy',
    #                           spacing=2000.
    #                           )
    # data_out = convert_interp('mantle_latlon.txt',
    #                           'nz_full_eberhart2015_mantle.xyz.npy',
    #                           spacing=8000.
    #                           )





