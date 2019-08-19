
"""
There are some values required by the internal mesher that control the sizes
of spectral elements in the horizontal and vertical directions. This collection
of functions should help determine those values in a standardized manner
"""
import numpy as np
from pyatoa.utils.operations.calculations import myround
from pyatoa.utils.operations.source_receiver import lonlat_utm


def provide_mesh_par_info():
    """
    Determine the number of vertical layers for a given number of doublings such
    that the top layer contains roughly cube shaped elements
    """
    # User Parameters
    lat_min = -43.
    lat_max = -37.
    lon_min = 173.
    lon_max = 179.
    utm_zone = -60
    shortest_period_s = 10.

    # Find the mesh dimensions in units of km
    x_min, y_min = lonlat_utm(lon_min, lat_min, utm_zone)
    x_max, y_max = lonlat_utm(lon_max, lat_max, utm_zone)
    x_length_km = (x_max - x_min) * 1E-3
    y_length_km = (y_max - y_min) * 1E-3

    ratio_x_to_y = x_length_km / y_length_km
    print(f"The ratio of x to y is {ratio_x_to_y}")

    # Tell the user the smallest value of NEX for shortest period
    # Equation taken from the Specfem3D manual
    nex_min = 288 * 2 / shortest_period_s
    print(f"Smallest number of elements for {shortest_period_s}s is {nex_min}")

    return x_length_km, y_length_km


def num_elements_horizontal(dx, dy, shortest_period_s, error_pct=0.95):
    """
    Determine the the number of elements along the two horizontal directions
    These are determined based on a few constraints

    Constraints:
    1: to give (relatively) square elements, the ratio of element sizes nx/ny
    should be equal to the ratio of mesh side lengths dx/dy, i.e.
        nx/ny = dx/dy
    2: the shortest period resolved by a mesh controls the minimum number of
    elements for a side length. This is given by Eq. 3.1 in the Specfem3D
    Cartesian manual, which is:
        shortest_period_s = 2 * (288/nx)
    3: nx and ny must be integer factors of 8, which is mandated by the division
    of elements to different slices which are then divvied up to the various
    processors

    :type dx: float
    :param dx: length of the x-dimension of the mesh
    :type dy: float
    :param dy: length of the y-dimension of the mesh
    :type shortest_period_s: int
    :param shortest_period_s: shortest resolvable target period in seconds
    :type error_pct: float
    :param error_pct: The allowable difference related to constraint 1, i.e.
        (nx/nx) / (dx/dy) >= error_pct
    """
    # Set the constraints
    target_ratio = dx/dy
    minimum_n = 2 * 288 / shortest_period_s

    # Loop over values of nx and ny until the error percentage criteria is met
    nx = minimum_n
    while True:
        # Get a good starting value that is a factor of 8
        nx = myround(nx, base=8, choice='up')
        # Find the nearest value of ny
        ny = myround(nx/target_ratio, base=8, choice='near')
        trial_ratio = nx/ny

        # If the trial ratio meets the goal, return nx, ny, else increase by 8
        if trial_ratio / target_ratio >= error_pct:
            return nx, ny
        else:
            nx += 8


def two_doubling_layers(alpha, mesh_depth):
    """
    Determine the number of layers for each doubling layer. I found it too
    difficult to generalize this problem in a clean way, so this only works
    if ndoublings = 2

    if this were to be generalized, each layer i would have a coefficient 2**i
    e.g. for i in range(0, ndoublings+1, 1):
             n = 2 ** i


    This is controlled by the desire to keep elements cube-shaped,
    i.e. at the top layer, the element size in the horizontal controls the
    vertical element length, and at each doubling layer, which coarsens the
    horizontal element size (by 2), the vertical layer should also coarsen.
    :type alpha: float
    :param alpha: dx/nx (horizontal element size), must have units of mesh_depth
    :type mesh_depth: float
    :param mesh_depth: desired depth extent of the mesh, same units as alpha
    """
    # The number of layers is controlled by the mesh depth and target resolution
    target_layers_n = myround(mesh_depth / alpha, base=1, choice='up')

    # Provide all possible combinations of layers, which is controlled by
    # the fact that with each new doubling layer, the vertical distance of each
    # element should be doubled
    # Restrict that the top layer can only provide 1/3 of the total layers
    print("TOP\tMID\tBOT")
    for top_layer_n in range(1, target_layers_n // 3, 1):
        for mid_layer_n in range(1, target_layers_n // 3, 1):
            bottom_layer_n = target_layers_n / 4 - 2 * mid_layer_n - top_layer_n
            if bottom_layer_n >= target_layers_n // 3:
                print(f"{top_layer_n}\t{mid_layer_n}\t{bot_layer_n}")


def cut_topography(data, x_or_lon_min, y_or_lat_min, x_or_lon_max,
                   y_or_lat_max):
    """
    Take a topography data array and cut it to the correct size and shape of
    the mesh based on mesh dimensions and element spacing.


    :type data: np.array
    :param data: data to be cut, should  be in the format [x, y, z]

    :return:
    """
    # Determine where the data falls outside the bounds
    x_too_small = np.where(data[:, 0] < x_or_lon_min)[0]
    x_too_large = np.where(data[:, 0] > x_or_lon_max)[0]
    y_too_small = np.where(data[:, 1] < y_or_lat_min)[0]
    y_too_large = np.where(data[:, 1] > y_or_lat_max)[0]

    # Get rid of duplicates
    to_remove = np.unique(np.concatenate(
        (x_too_small, x_too_large, y_too_small, y_too_large), 0)
    )

    # delete in place
    data = np.delete(data, to_remove, 0)

    return data


def save_topo(data, tag="topo"):
    """
    Saves the topo as an npy and an ascii, with a descriptive file name

    :param data:
    :return:
    """
    # Find the unique array elements since the data is in a grid format
    x_unique = np.unique(data[:, 0])
    y_unique = np.unique(data[:, 1])

    # Get some unique identifiers
    nx = len(x_unique)
    ny = len(y_unique)
    x0 = x_unique[0]
    y0 = y_unique[0]
    dx = x_unique[1] - x_unique[0]
    dy = y_unique[1] - y_unique[0]

    fid = f"{tag}_{x0:.0f}_{y0:.0f}_nx{nx}_ny{ny}_dx{dx:.3e}_dy{dy:.3e}"
    np.save(fid + '.npy', data)

    # Assuming the data structure is already in the correct format, that is
    # x is iterated over before y
    np.savetxt(fid + '.txt', data[:, 2], fmt='%d')


def interfaces():
    """
    Just determine the layer number that corresponds to a depth
    """




