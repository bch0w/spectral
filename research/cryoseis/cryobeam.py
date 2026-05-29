"""
Beamformer tuned for 2024 Gulkanaseis deployment
Beamforming source code inspired by: https://github.com/schipp/fast_beamforming
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import torch
from itertools import product
from obspy import UTCDateTime, read, Stream
from pyproj import Transformer


def lonlat_to_utm(lon, lat):
    """Convert lon/lat to UTM coordinates"""
    zone = int((lon + 180) / 6) + 1
    hemisphere = "N" if lat >= 0 else "S"

    crs_utm = f"EPSG:{32600 + zone if hemisphere == 'N' else 32700 + zone}"
    transformer = Transformer.from_crs("EPSG:4326", crs_utm, always_xy=True)
    easting, northing = transformer.transform(lon, lat)

    return easting, northing


def array_geometry(plot=False):
    """
    Define array geometry in UTM coordinates. Set the origin to be the center
    of the array
    """
    # Gulkana Seismic Deployment 2024
    stations = {
        "100": (63.265612, -145.415025),
        "101": (63.265362, -145.41407),
        "102": (63.265448, -145.41381),
        "103": (63.265203, -145.415355),
        "104": (63.265062, -145.414323),
        "105": (63.26511,	-145.413887),
        "106": (63.26495,	-145.415472),
        "107": (63.264683, -145.414657),
        "108": (63.264878, -145.414107),
        "109": (63.26392,	-145.41601),
        "110": (63.26466,	-145.41557),
        "111": (63.264278, -145.415858),
    }    


    # Convert to UTM coordinates so we are in a Cartesian reference frame
    stations_utm  =  {}
    xmax = 0
    xmin = np.inf
    ymax = 0
    ymin = np.inf
    for name, coords in stations.items():
        lat, lon = coords
        x, y = lonlat_to_utm(lon, lat)[:2]
        stations_utm[name] = (x, y)
        # Used to find array bounds for centering
        if x > xmax:
            xmax = x
        if x < xmin:
            xmin = x
        if y > ymax:
            ymax =  y
        if y < ymin:
            ymin = y

    # find array center
    xmiddle = (xmax + xmin) / 2
    ymiddle = (ymax + ymin) / 2

    for name, coords in stations_utm.items():
        x, y = coords
        stations_utm[name] = (x - xmiddle, y - ymiddle)

    if plot:
        for name, coords in stations.items():
            x, y = coords
            plt.scatter(x, y, c="k")
            plt.text(x, y, name)
        
        # Create empty cells on the mesh for clearer visualization
        # empty_cells = np.zeros_like(xx) * np.nan
        # plt.pcolormesh(xx, yy, empty_cells, ec="k", lw=.01)

        plt.xlabel("Easting [m]")
        plt.ylabel("Northing [m]")
        plt.title("Gulkanaseis 2024")
        plt.show()

    return stations_utm

def make_grid(xmin=-1E3, xmax=1E3, dx=100, ymin=-1E3, ymax=1E3, dy=100):
    """
    Create the grid search grid where the origin is defined by stations.
    This is for Matched Field Processing
    """
    grid_x = np.arange(xmin, xmax, dx)
    grid_y = np.arange(ymin, ymax, dy)
    xx, yy = np.meshgrid(grid_x, grid_y, indexing="xy")
    gridpoints = np.array(list(product(xx, yy)))
    return xx, yy


def get_data(station, component="Z", path="./data", fmin=10, fmax=250):
    """Cut waveform data from all the stations to use for beamforming"""
    # Time to cut
    starttime = UTCDateTime("2024-09-08T21:04:44")
    endtime = UTCDateTime("2024-09-08T21:04:50")

    fid = f"GS.{station}.DH{component}.{starttime.year}.{starttime.julday}"
    st = read(os.path.join(path, fid))

    st.trim(starttime, endtime)
    st.taper(.05)
    st.filter("bandpass", freqmin=fmin, freqmax=fmax)

    return st

def get_slowness_space(dbaz=4, smin=0, smax=0.5, ds=.02, plot=False):
    """
    Define slowness grid to search over

    Parameters
    ----------
    dbaz (float): discretization of backazimuth
    smin (float): minimum slowness value
    smax (float): maximum slowness value
    ds (float): discretization of slowness
    """
    azs = torch.arange(0, 2 * torch.pi, dbaz * 2 * torch.pi / 360)
    slows = torch.arange(smin, smax, ds)
    slowness_space = torch.tensor(
        [(torch.sin(az) * s, torch.cos(az) * s) for az, s in product(azs, slows)]
    )

    # visualise slowness space
    if plot:
        fig, axs = plt.subplots(1, 2)
        fig.subplots_adjust(wspace=0.3)

        # corresponding azs and slowns in slowness_space
        azs_calc = torch.atan2(slowness_space[:, 0], slowness_space[:, 1])
        slows_calc = torch.norm(slowness_space, dim=1)

        ax = axs[0]
        sct = ax.scatter(*slowness_space.T, c=azs_calc, s=10, lw=0, cmap="twilight")
        ax.set(aspect="equal", xlabel="x [s/km]", ylabel="y [s/km]")
        x0, y0, w, h = ax.get_position().bounds
        cbar_ax = fig.add_axes([x0 + w + 0.01, y0, 0.01, h])
        cbar = fig.colorbar(sct, cax=cbar_ax, label="Backazimuth [rad]")

        ax = axs[1]
        sct = ax.scatter(*slowness_space.T, c=slows_calc, s=10, lw=0, cmap="plasma_r")
        ax.set(aspect="equal", yticklabels=[], xlabel="x [s/km]")
        x0, y0, w, h = ax.get_position().bounds
        cbar_ax = fig.add_axes([x0 + w + 0.01, y0, 0.01, h])
        cbar = fig.colorbar(sct, cax=cbar_ax, label="Slowness [s/km]")
        plt.show()

    return slowness_space



def main():
    fmin = 1
    fmax = 10

    stations = array_geometry(plot=False)
    slowness_space = get_slowness_space(plot=False)
    # plot(stations)

    # Get filtered data
    st = Stream()
    coords = []
    for station, coord in stations.items():
        st += get_data(station, fmin=fmin, fmax=fmax)
        coords.append(coord)

    # Reorganize data so it is in the correct order
    data = torch.tensor(np.array([tr.data for tr in st]))
    coords = torch.tensor(coords)

    sampling_rate = st[0].stats.sampling_rate
    window_length = st[0].stats.npts / sampling_rate

    # Begin beamforming
    
    # reference point for relative plane-wave traveltimes
    # arbitrary but must be reasonably nearby
    reference_point = torch.mean(coords, axis=0)
    sampling_rate = st[0].stats.sampling_rate

    # Compute expected frequencies based on data length and sampling rate
    f = torch.fft.fftfreq(int(window_length * sampling_rate), 
                          d=1 / sampling_rate)
    omega = 2 * torch.pi * f
    omega_limited = omega[torch.where((f > fmin) & (f < fmax))]
    
    waveforms_prep = torch.zeros(
        (1, data.shape[0], int(window_length * sampling_rate))
    )
    waveforms_prep[0] = data[:, :]

    #  Compute spectra of waveforms
    data_spectra_all = torch.fft.fft(waveforms_prep, dim=-1)

    data_spectra_all_limited = \
        data_spectra_all[:, :, torch.where((f > fmin) & (f < fmax))[0]]

    traveltimes_space = torch.einsum("nx, sx -> sn", 
                                     reference_point - coords, 
                                     slowness_space)
    synth_spectra_space = torch.exp(
        -1j * torch.einsum("sn, w -> snw", traveltimes_space, omega_limited)
    )

    S = torch.einsum(
        "siw, sjw -> sijw", synth_spectra_space, torch.conj(synth_spectra_space)
    )

    # compute the cross-spectra between all stations for all time windows.
    # notice the additional dimension for time windows t.
    K = torch.einsum(
        "tiw, tjw -> tijw", data_spectra_all_limited, torch.conj(data_spectra_all_limited)
    )

    diag_idxs = torch.arange(K.shape[1])
    zero_spectra = torch.zeros(omega_limited.shape, dtype=torch.cfloat)

    K[:, diag_idxs, diag_idxs, :] = zero_spectra

    beampowers_for_time = torch.einsum("sjiw, tijw -> ts", S, K).real
    beampowers_for_time.shape

    # plot full beampower distribution for the first time window for quality check
    bp_to_plot = beampowers_for_time[0]

    # replace all nans with 0s
    bp_to_plot[torch.isnan(bp_to_plot)] = 0
    # convert to dB (negative values become nan)
    bp_to_plot = 10 * torch.log10(bp_to_plot / torch.max(bp_to_plot))

    fig, ax = plt.subplots()
    sct = ax.scatter(
        *slowness_space.T,
        c=bp_to_plot,
        s=50,
        lw=0,
        cmap="Purples",
        vmin=-5,
        vmax=0,
    )

    ax.set(
        xlim=(-0.5, 0.5), ylim=(-0.5, 0.5), aspect="equal", title="Beampower distribution"
    )

    # mark some slownesses for orientation
    for radius in [0.1, 0.3, 0.5]:
        circle = plt.Circle((0, 0), radius, color="#555", fill=False, ls="--", lw=0.5)
        ax.add_artist(circle)
        # add slowness labels
        ax.text(
            radius * 0.7,
            radius * 0.7,
            f"{radius:.1f} s/km",
        )

    x0, y0, w, h = ax.get_position().bounds
    cbar_ax = fig.add_axes([x0 + w + 0.01, y0, 0.01, h])
    cbar = fig.colorbar(sct, cax=cbar_ax, label="Beampower [dB]")
    plt.show()


if __name__ == "__main__":
    main()