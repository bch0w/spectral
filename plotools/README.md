# PLOT TOOLS

Quick instructions for how to use these plotting routines for rapidly and 
flexibily visualizing data. 

## Install

To install you can clone this repo and then make a Conda environment with ObsPy

```bash
git clone https://github.com/bch0w/spectral.git
conda create -n spectral
conda activate spectral
conda install obspy
python spectral/plottols/prettyplot.py -h
```

## PrettyPlot

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

![ex1](https://private-user-images.githubusercontent.com/23055374/488063348-214f1146-824a-48e1-af10-f0ebd477dcf5.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NTc1NDYxNjksIm5iZiI6MTc1NzU0NTg2OSwicGF0aCI6Ii8yMzA1NTM3NC80ODgwNjMzNDgtMjE0ZjExNDYtODI0YS00OGUxLWFmMTAtZjBlYmQ0NzdkY2Y1LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA5MTAlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwOTEwVDIzMTEwOVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTRmMmE1MTNlOGU4YTU5Y2Y0MTU0OTViOGIwYTJkYzc5YzMwMjhlNmVkOGZlOGU3M2M0NGI1NmNiMTFjNTRlZWImWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.agWBJnR2kYPZgpPdFlNJ4WkW4mXYUdLgI0dI15zgisU)

#### Multiwaveform
Plot multiple waveforms on the same axis, set the colors of each using the CN 
colors of matplotlib (C0, C1 ... CN), use relative scale in seconds on the time 
axis, resample to 25 Hz, bandpass 1--10Hz.

```bash
python prettyplot.py GS.109.DHZ.2024.253 GS.109.DHN.2024.253 GS.109.DHE.2024.253 \
    --colors 'C?' --resample 25  --taper 0.05 --fmin 1 --fmax 10 --time 'a-08' \
    --xlim 2024-09-08T17:40:45 2024-09-08T17:41:30 --save multiwaveform.png 
```

![ex2](https://private-user-images.githubusercontent.com/23055374/488063618-29cec8a0-4728-4c6f-8b61-5b53c3ddf67c.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NTc1NDYxNjksIm5iZiI6MTc1NzU0NTg2OSwicGF0aCI6Ii8yMzA1NTM3NC80ODgwNjM2MTgtMjljZWM4YTAtNDcyOC00YzZmLThiNjEtNWI1M2MzZGRmNjdjLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA5MTAlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwOTEwVDIzMTEwOVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTMxZWIzMDJmODNmN2ExZjI0NWFlNzc1YjJhMGYyYWY4ZTI2YmEwNDc1NTg3OTk0Mjk5Y2NkYzRmY2VkZTUyM2YmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.cYXujv7vuYYaACJTwV95sLX-Qi9svvI6-mkfHloHjp8)

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

![ex3](https://private-user-images.githubusercontent.com/23055374/488063805-88523a65-79cd-46c6-860a-2753fde4b359.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NTc1NDYxNjksIm5iZiI6MTc1NzU0NTg2OSwicGF0aCI6Ii8yMzA1NTM3NC80ODgwNjM4MDUtODg1MjNhNjUtNzljZC00NmM2LTg2MGEtMjc1M2ZkZTRiMzU5LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA5MTAlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwOTEwVDIzMTEwOVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWZjM2MzMmE5NWNkYzZkN2VmZmU3ZDllNDcyMDA1MWY4Y2U0MDA1YzRlOTIxNjdjNjFmY2QwZjJhNTkxOTM4MGImWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.DQJ5D-HlPZCPEBGOsAKvCNi_Ozbr-K7TyuYKItAG2H4)

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

![ex4](https://private-user-images.githubusercontent.com/23055374/488064007-8f3e7707-9850-43d2-8b9c-18f9dd8d26c5.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NTc1NDYxNjksIm5iZiI6MTc1NzU0NTg2OSwicGF0aCI6Ii8yMzA1NTM3NC80ODgwNjQwMDctOGYzZTc3MDctOTg1MC00M2QyLThiOWMtMThmOWRkOGQyNmM1LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA5MTAlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwOTEwVDIzMTEwOVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWJkMzJkNzQyOTg1MWVkZTUzMzQ5MDAwODYxMTI4YWNhYTNmNmM5Yjc0M2E0MGMwZTI4MDU0NDU4NzNkODA5MzMmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.t8FfwCn86rPIVsMKGS70Cs2PauLtT_4zBmR2sRXVD9c)
