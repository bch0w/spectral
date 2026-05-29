"""
https://github.com/schipp/fast_beamforming/blob/main/beamforming_planewave_data.ipynb
"""
# download a few days of gräfenberg array data
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

import matplotlib.pyplot as plt
import torch
from obspy.signal.util import util_geo_km
from itertools import product


client = Client("BGR")
t1 = UTCDateTime("2020-12-20T00:00:00.000")
t2 = UTCDateTime("2020-12-24T00:00:00.000")
st_orig = client.get_waveforms("GR", "GR*", "*", "LHZ", t1, t2)

# download metadata for gräfenberg
inv = client.get_stations(
    network="GR",
    station="GR*",
    location="*",
    channel="LHZ",
    starttime=t1,
    endtime=t2,
    level="response",
)

# remove instrument response
st_orig.remove_response(inventory=inv)

# pre-process seismograms for reasonable stable results
st = st_orig.copy()
st.detrend("demean")
st.detrend("linear")
st.filter("highpass", freq=0.01, zerophase=True)
st.taper(max_percentage=0.01, type="hann")

# force all traces to have same length
# 0s are not an issue for beamforming
st.merge(fill_value=0)
st.trim(t1, t2, pad=True, fill_value=0)

# prepare waveforms tensor
waveforms = torch.tensor([tr.data for tr in st])

# stations that data is present for
stations = [tr.stats.station for tr in st]

# extract locations in same order from inv
coordinates = []
station_order = []
for station in stations:
    for net in inv:
        for sta in net:
            if sta.code == station:
                coordinates.append((sta.longitude, sta.latitude))
                station_order.append(sta.code)
coordinates = torch.tensor(coordinates)

# convert coordinates to cartesian, relative to geographical array center
coordinates_cartesian = []
for sta in coordinates:
    coordinates_cartesian.append(
        util_geo_km(
            torch.mean(coordinates[:, 0]), torch.mean(coordinates[:, 1]), sta[0], sta[1]
        )
    )
coordinates_cartesian = torch.tensor(coordinates_cartesian)

fig, ax = plt.subplots()
ax.scatter(*coordinates_cartesian.T, marker="v")
ax.set_aspect("equal")
for i, txt in enumerate(station_order):
    ax.annotate(txt, (coordinates_cartesian[i, 0], coordinates_cartesian[i, 1]))
ax.set(title="Gräfenberg array", xlabel="x [km]", ylabel="y [km]")


# define slowness space to explore
# defined in backazimuth [º] and abs(slowness) [s/km]
backazimuth_spacing = 4
slowness_min = 0
slowness_max = 0.5
slowness_spacing = 0.02

azs = torch.arange(0, 2 * torch.pi, backazimuth_spacing * 2 * torch.pi / 360)
slows = torch.arange(slowness_min, slowness_max, slowness_spacing)
slowness_space = torch.tensor(
    [(torch.sin(az) * s, torch.cos(az) * s) for az, s in product(azs, slows)]
)

# visualise slowness space
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

# Beamforming setup

# Primary microseism band
fmin = 0.05
fmax = 0.1

# Define time windows to beamform in
# allows for overlapping windows

times = torch.tensor(st[0].times())
window_length = 2 * 3600  # sec
window_limits = 0, times[-1] - window_length
window_overlap = 0.5

# prepare beamforming in time windows moving over signals
# and restricting data to frequency band of interest

# windows to beamform in
windows = torch.arange(
    window_limits[0],
    window_limits[1] - window_length,
    window_length * (1 - window_overlap),
)

print(f"{windows.shape=}")

# reference point for relative plane-wave traveltimes
# arbitrary but must be reasonably nearby
reference_point = torch.mean(coordinates_cartesian, axis=0)

# define start and end indices of each window
sampling_rate = st[0].stats.sampling_rate
starttime_idxs = (windows * sampling_rate).int()
endtime_idxs = starttime_idxs + int(window_length * sampling_rate)

print(f"{starttime_idxs.shape=}")

# compute expected frequencies, and limit to range of interest
f = torch.fft.fftfreq(int(window_length * sampling_rate), d=1 / sampling_rate)

omega = 2 * torch.pi * f
omega_limited = omega[torch.where((f > fmin) & (f < fmax))]
print(f"{omega_limited.shape=}")

# cut timewindows from correlations
waveforms_prep = torch.zeros(
    (len(windows), waveforms.shape[0], int(window_length * sampling_rate))
)

print(f"{waveforms_prep.shape=}")
for window_idx, (starttime_idx, endtime_idx) in enumerate(
    zip(starttime_idxs, endtime_idxs)
):
    waveforms_prep[window_idx] = waveforms[:, starttime_idx:endtime_idx]
# compute data spectra and limit to frequency band of interest
data_spectra_all = torch.fft.fft(waveforms_prep, dim=-1)

print(f"{data_spectra_all.shape=}")
data_spectra_all_limited = data_spectra_all[
    :, :, torch.where((f > fmin) & (f < fmax))[0]
]
print(f"{data_spectra_all_limited.shape=}")


# BEAMFORMING
traveltimes_space = torch.einsum(
    "nx, sx -> sn", reference_point - coordinates_cartesian, slowness_space
)

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
breakpoint()

plt.show()