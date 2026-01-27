# PLOT TOOLS

Quick instructions for how to use these plotting routines for rapidly and 
flexibily visualizing data. 

## Install

All of these instructions should be done in a terminal or Windows Powershell.

1. Install the Python package manager `Micromamba`, [instructions for different operating systems can be found here](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) 
    > Note: You may also download Miniconda, but I suggest Micromamba. If you go with Conda, you will have to replace all the `mamba` commands in these instrucftions with `conda`.

2. Create a new Conda environment by running the following and answering the prompts:
    ```bash
    mamba env create prettyplot obspy
    ````
3. Activate your Mamba environment (you will have to do this everytime you open a new terminal)
    ```bash
    mamba activate prettyplot
    ```
4. Save the plotting script from my GitHub repository wherever you want. You can swap out the name if you want access to other scripts. 

    ```bash
    curl https://raw.githubusercontent.com/bch0w/spectral/master/plotools/prettyplot.py -o prettyplot.py
    ```

5. Ensure that the script works by bringing up the help message
    ```bash
    python prettyplot.py -h
    ```

6. If you got to this stage and the help message worked then you can run the script. If you want to be able to run the script anywhere on your computer, you can get the general command with

    ```bash
    echo $(which python) $(pwd)/prettyplot.py
    ```
Whatever is returned, you can copy-paste that in your terminal anywhere and you should be able to plot things.


## PrettyPlot Usage

PrettyPlot makes waveforms and (optionally) spectrograms on the same figure. Some options for modest processing. Requires only ObsPy and it's dependencies.

Call from the command line

```bash
python prettyplot.py -h  # calls help message
```

Use command line arguments to set options:
```bash
python prettyplot.py <path/to/fid> --resample 100 --fmin 1 --fmax 10
```

> NOTE: `<path/to/fid>` should point to your seismic data, which can be in any [format that ObsPy recognizes](https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#obspy.core.stream.read). Your file should only be a single component, otherwise PrettyPlot will only plot the component in the first index.

### Features

1. TauP: Add theoretical phase arrivals using a TauP model. All `tp-` 
    parameters are related
2. Spectrogram: Add a spectrogram for the time period you are plotting, alongside waveform figure, with the same x-axis. 
    All `sp-` parameters are related.
3. Tmarks: Put vertical lines at user-defined times to mark events/arrivals/etc.
4. Time Axis: Allows for relative time (sec/min/hr) or absolute time 
    (date + time). Absolute time can also be shifted to match local time.
5. (DEV) Stream Gage: For Bryant's Gulkanaseis project only, plots USGS stream 
    gage data with the waveform data. Might be useful example for plotting other
    geophysical data.

---
### Examples

Images hosted [here](https://github.com/bch0w/spectral/issues/1#issuecomment-3276805610) incase links change.

#### Waveform 
Read file from MiniSeed (response already removed). 0.5Hz or 2s highpass, use absolute 
time converted to UTC-08 (AKDT), taper the waveform by 5%, trim the data
(in local time) to a 2 minute section, save to current dir. with auto-generated
filenames

```bash
python prettyplot.py GS.109.DHZ.2024.253 --fmin 0.5 --time 'a-08' --taper 0.05 \
    --xlim 2024-09-08T17:40:00 2024-09-08T17:42:00 --save 'auto' 
```
<img width="1600" height="800" alt="ex1" src="https://github.com/user-attachments/assets/5bade01b-721c-4160-a3f1-07cb21a28fc1" />


#### Multiwaveform
Plot multiple waveforms on the same axis, set the colors of each using the CN 
colors of matplotlib (C0, C1 ... CN), use relative scale in seconds on the time 
axis, resample to 25 Hz, bandpass 1--10Hz.

```bash
python prettyplot.py GS.109.DHZ.2024.253 GS.109.DHN.2024.253 GS.109.DHE.2024.253 \
    --colors 'C?' --resample 25  --taper 0.05 --fmin 1 --fmax 10 --time 'a-08' \
    --xlim 2024-09-08T17:40:45 2024-09-08T17:41:30 --save multiwaveform.png 
```
<img width="1600" height="800" alt="ex2" src="https://github.com/user-attachments/assets/ffffc492-5aac-4445-889f-d809b3f885e0" />


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
<img width="1600" height="1200" alt="ex3" src="https://github.com/user-attachments/assets/9a00b8a5-796f-4ab5-a3f2-24cbc54729e7" />


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
<img width="1600" height="800" alt="ex4" src="https://github.com/user-attachments/assets/465fa6a4-3075-4b58-bfd4-459cc517a38e" />

