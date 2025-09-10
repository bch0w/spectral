# Plotters

Quick instructions for how to use these plotting routines for rapidly and 
flexibily visualizing data:

## prettyplot.py

Makes waveforms and (optionally) spectrograms on the same figure. Some options
for modest processing. Requires only ObsPy and it's dependencies.

### Usage

Call from the command line

```bash
python prettyplot.py -h  # calls help message
```

Use command line arguments to set options:
```bash
python prettyplot.py <path/to/fid> --resample 100 --fmin 1 --fmax 10
```

### Features

1. TauP: Add theoretical phase arrivals using a TauP model. All `tp-` 
    parameters are related
2. Spectrogram: Add a spectrogram for the time period you are plotting, along-
    side waveform figure, with the same x-axis. 
    All `sp-` parameters are related.
3. Tmarks: Put vertical lines at user-defined times to mark events/arrivals/etc.
4. Time Axis: Allows for relative time (sec/min/hr) or absolute time 
    (date + time). Absolute time can also be shifted to match local time.
5. (DEV) Stream Gage: For Bryant's Gulkanaseis project only, plots USGS stream 
    gage data with the waveform data. Might be useful example for plotting other
    geophysical data.

### Examples

#### Waveform 
Read file from MiniSeed (response already removed). Use absolute 
time converted to UTC-08 (AKDT), taper the waveform by 5%, trim the data
(in local time) to a 2 minute section, save to current dir. with auto-generated
filenames

```bash
python prettyplot.py GS.109.DHZ.2024.253 --fmin 0.5 --time 'a-08' --taper 0.05 \ 
    --xlim 2024-09-08T17:40:00 2024-09-08T17:42:00 --save 'auto' 
```

#### Multiwaveform
Plot multiple waveforms on the same axis, set the colors of each using the CN 
colors of matplotlib (C0, C1 ... CN), use relative scale in seconds on the time 
axis, resample to 25 Hz, bandpass 1--10Hz.

```bash
python prettyplot.py GS.109.DHZ.2024.253 GS.109.DHN.2024.253 GS.109.DHE.2024.253 \
    --colors 'C?' --resample 25  --taper 0.05 --fmin 1 --fmax 10 --time 'a-08' \
    --xlim 2024-09-08T17:40:45 2024-09-08T17:41:30 --save multiwaveform.png 
```

#### Waveform + Spectrogram
Plot waveform and add spectrogram on top. Set the amplitudes of the PSDs in 
decibel-scale, set the frequency axis of the spectrogram in log scale, 
use the `nipy_spectral` colormap for spectrogram

```bash
python prettyplot.py GS.109.DHZ.2024.253 \
    --fmin 0.5 \
    --time 'a-08' \
    --taper 0.05 \
    --xlim 2024-09-08T17:40:00 2024-09-08T17:42:00 \
    --spectrogram \
    --sp_dbscale \
    --sp_logscale \
    --sp_cmap 'nipy_spectral' \
    --save 'auto' 
```

#### Time Marks
Sets arbitrary vertical lines with distinct colors to mark user-defined portions
of the waveform. 

```bash
python prettyplot.py GS.109.DHZ.2024.253  \
    --time 'a-08' \
    --xlim 2024-09-08T17:40:00 2024-09-08T17:42:00 \
    --tmarks 2024-09-08T17:40:10 2024-09-08T17:40:20 2024-09-08T17:40:30 \
    --tmarks_c r b g
```
