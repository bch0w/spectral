# ISHTAR data mine paper: seismic noise and wind at co-located stations

Working notes for the analysis behind the data mine paper. It records what
each script does, how the paper figures were made, what to be careful
about, and where the work stopped. Written November 2026.

**Start with [Where this stands](#where-this-stands) when coming back to the
project.**

## Contents

- [Where this stands](#where-this-stands)
- [The experiment and the data](#the-experiment-and-the-data)
- [How the code is organized](#how-the-code-is-organized)
- [Paper figures and how to remake them](#paper-figures-and-how-to-remake-them)
- [Scripts: PPSD pipeline](#scripts-ppsd-pipeline)
- [Scripts: wind-state pipeline](#scripts-wind-state-pipeline)
- [Dependencies outside this repo](#dependencies-outside-this-repo)
- [Caveats](#caveats)
- [Results worth remembering](#results-worth-remembering)
- [Reproducibility record](#reproducibility-record)

## Where this stands

**Chosen for the paper** (copies in `windygrams/FIGURES/FOR_PAPER/` in the
working directory):

1. `k24k_wind_data.png`: K24K Covered, Uncovered, Windscreen and Permanent
   during rising wind. This is `wind_covers.png` from `plot_wind_states.py`,
   renamed.
2. `wind_compare_mitdu_mitdc_hvmitd_1hz.png`: MITDU, MITDC and HV.MITD
   during rising and high wind, with STA/LTA transients marked.

The PPSD figures in `ppsd/` (spectra relative to the permanent station,
composite PPSDs, wind summary, wind roses) are drafts or candidates;
[`ppsd/FIGURE_PLAN.md`](ppsd/FIGURE_PLAN.md) has their status and decisions.

**To do before submission**

- [ ] Check the text and captions once the real network code is assigned.
      The figures need no change: the temporary `XX` code is left out of
      every label, and permanent stations keep their real codes (e.g.
      `AK.K24K..BHZ`).
- [ ] Write the caption for the MITD comparison figure. Facts it needs are
      [below](#mitd-comparison-figure).
- [ ] Decide whether to mention the MITDC-only spike near 17:04 UTC in
      panel (b) of that figure.
- [ ] Say "much less affected" and not "unaffected" for HV.MITD. Its
      spectrogram still brightens above about 10 Hz with the wind.

**Open science questions**

- [ ] **Narrowband 20-45 Hz lines on the nodes.** Waiting on temperature
      data from colleagues. See
      [Narrowband lines](#narrowband-lines-unresolved).
- [ ] Next steps for those lines: night-time peak-frequency tracking per
      node from the cached PSDs (does wander go U > C > W, is the phase tied
      to sunrise, clear vs cloudy days), check the MITDW and DESDW (30 April)
      field notes, then correlate with temperature.
- [ ] Is the low wind at H23KA, MDMA, PAXA and RNDA real or an anemometer
      problem? They never exceeded 2.6-2.9 m/s in hourly means.
- [ ] Check the field notes for KAEDW: its HHE sits about 8 dB below KAEDC
      HHE with a frequency-dependent, time-varying offset that looks like
      the sensor component. Currently kept in the figure and noted in the
      text.

**Loose ends in the code**

- [ ] `pysep/pplot.py` changes are on PySEP branch
      `pplot-envelope-spectrogram`, not yet merged into `master`. The
      wind-state figures need them. See
      [Dependencies](#dependencies-outside-this-repo).
- [ ] The `psd_viz` fixes are on a local-only branch, also needed.
- [ ] Older figures in `windygrams/FIGURES/<condition>/` predate the
      antialiased spectrograms, the title placement and the `XX`-free
      labels. Re-run `windygrams.py` if any of them gets reused.

## The experiment and the data

Temporary nodal stations were installed next to permanent stations at eight
sites to compare noise and how wind affects it. Each site has:

- `<SITE>C`: Covered node
- `<SITE>U`: Uncovered node
- `<SITE>W`: node in a windscreen (always say "windscreen")
- `<SITE>A`: anemometer, wind speed `LWS` and direction `LWD` at 1 Hz
- the permanent station, with its own network and channel codes

The channel and common time window of each site:

| Site | Node | Permanent | Common window (2025) |
|------|------|-----------|----------------------|
| DESD | HH, 100 Hz | HV.DESD EH | 26 Mar to 30 Apr |
| KAED | HH, 100 Hz | HV.KAED EH (Z only, looks broken) | 26 Mar to 29 Apr |
| MITD | HH, 100 Hz | HV.MITD HH, 100 Hz | 26 Mar to 30 Apr |
| H23K | DH, 250 Hz | AK.H23K BH, 50 Hz | 9 Jul to 23 Jul |
| MDM | DH, 250 Hz | AK.MDM BH, 50 Hz | 10 Jul to 23 Jul |
| K24K | EH, 100 Hz | AK.K24K BH, 50 Hz | 2 Sep to 3 Oct |
| PAX | EH, 100 Hz | AK.PAX BH, 50 Hz | 2 Sep to 3 Oct |
| RND | EH, 100 Hz | AK.RND BH, 50 Hz | 5 Sep to 3 Oct |

The AK permanent stations have no data above 25 Hz, which is why the K24K
figure stops at 22.5 Hz. Some AK day files also contain `LHZ`, so scripts
select the channel after reading.

**Data are not in this repo.** On the original machine:

- `/Users/prof/Work/akvenus/MSEEDS`: day files named
  `NET.STA.LOC.CHA.YEAR.JDAY`, e.g. `XX.K24KU..EHZ.2025.261`
- `/Users/prof/Work/akvenus/STATIONXML`: `ISHTAR.xml` for the nodes and
  `NET.STA.xml` for each permanent station
- `ppsd/SPECTRA`: combined PPSD `.npz` files from `run_psd_viz.py`
- `windygrams/WIND`: 60 s mean wind files, made on demand

Julian days are 2025 days of the year: day 95 is 5 April, day 100 is
10 April, day 117 is 27 April, day 261 is 18 September.

## How the code is organized

```text
datamine_paper/
  README.md             this file
  setup_links.sh        recreates the working directory symlinks
  ppsd/                 long-term spectra figures
    FIGURE_PLAN.md      figure roadmap, decisions and status
  windygrams/           wind-state waveform and spectrogram figures
    configs/            YAML configs for windygrams.py
    checks/             one-off shell scripts
    results/            small derived tables that figures read
  external/             patches for the other repos this work depends on
```

**Scripts run from working directories, through symlinks.** This repo holds
code; `/Users/prof/Work/akvenus/ppsd` and `.../windygrams` hold links to the
data and all outputs. Each working directory has symlinks to the scripts in
this repo. Editing either path changes the same file.

Most windygrams scripts change into the directory of their own symlink when
they start, and all use relative paths such as `DATA` and `XML`. **So run
them through the symlink, not from the repo, or they will not find data.**
The `ppsd/` scripts use the same relative paths but do not change directory,
so run them from the ppsd working directory. On a new machine or after moving
files, run:

```bash
./setup_links.sh /path/to/akvenus
```

It expects `MSEEDS` and `STATIONXML` directly inside that directory, and
never overwrites a real file.

**Two Python environments are used:**

- `/Users/prof/micromamba/envs/work`: everything in `windygrams/`, and
  `plot_wind.py`. Python 3.14.4, ObsPy 1.5.0, NumPy 2.4.3, SciPy 1.17.1,
  Matplotlib 3.10.8, pandas 3.0.2, PyYAML 6.0.3, with PySEP installed
  from `/Users/prof/Repos/pysep`
- `/Users/prof/micromamba/envs/psd_viz`: everything in `ppsd/`. The plotting
  scripts import `watc_ppsd`, which only this environment has

The shell alias `pplot` points to an old prettyplot. Use
`/Users/prof/micromamba/envs/work/bin/pplot` for PySEP's version.

## Paper figures and how to remake them

### Wind-state figures

All are made by one command, from the windygrams working directory:

```bash
/Users/prof/micromamba/envs/work/bin/python plot_wind_states.py
```

Output goes to `FIGURES/paper/`. The script does not write to
`FIGURES/FOR_PAPER/`, which was filled by hand, so copy a remade figure
there if it changes.

| File | Content | Windows set by |
|------|---------|----------------|
| `wind_covers.png` | K24K C, U, W, P, rising wind. **In paper** as `k24k_wind_data.png` | `COVERS_K24K` |
| `wind_compare_mitdu_mitdc_hvmitd_1hz.png` | MITDU, MITDC, HV.MITD, rising and high wind, marks. **In paper** | `STATES_MITD_TRANSIENTS` |
| `wind_compare_mitdu_hvmitd_1hz.png` | Same without MITDC | `STATES_MITD_TRANSIENTS` |
| `wind_covers_mitd.png` | MITD C, U, W, P, rising wind | `COVERS_MITD` |
| `wind_states.png` | K24KU in four wind states | `STATES` |
| `wind_states_mitd.png` | MITDU in four wind states | `STATES_MITD` |
| `wind_states_stations.png` | Four wind states, a different site each | `STATES_STATIONS` |
| `wind_states_mitd_transients.png` | MITDU, windows rich in transients, 0.1 Hz | `STATES_MITD_TRANSIENTS` |
| `wind_states_mitd_transients_1hz.png` | Same at 1 Hz, with marks | `STATES_MITD_TRANSIENTS` |
| `wind_states_mitdc_*`, `_mitdw_*`, `_hvmitd_*` | Same for the other MITD sensors | `MARKED_FIGURES` |

The windows are constants at the top of `plot_wind_states.py`.

#### K24K figure

Window: 18 September 2025 (day 261), 14:10-18:10 UTC. Caption used:

> Seismic response to rising wind at co-located sensors at site K24K on 18
> September 2025 (14:10-18:10 UTC): (a) covered node, (b) uncovered node,
> (c) node with windscreen, and (d) permanent station AK.K24K. Each panel
> shows a spectrogram (top) above the vertical-component ground velocity
> (bottom, black), with the instrument response removed and a 0.1 Hz
> highpass filter applied. Spectrograms are truncated at 0.5 Hz to exclude
> the strongest microseism energy, and at 22.5 Hz, just below the Nyquist
> frequency of the permanent station. Red lines show wind speed from the
> co-located anemometer, recorded at 1 Hz and smoothed with a centered 60 s
> running mean. Axis limits and color scales are the same in all panels. As
> wind speed increases from about 1 to 8-10 m/s after 15:00 UTC, velocity
> amplitudes and power above about 1 Hz grow most strongly on the uncovered
> node, moderately on the covered node, and least on the windscreen node
> and permanent station.

#### MITD comparison figure

- Top row: rising wind, day 100 (10 April 2025), 13:50-17:50 UTC
- Bottom row: high wind, day 95 (5 April 2025), 07:20-11:20 UTC
- Columns: MITDU, MITDC, HV.MITD; vertical component
- Processing: response removed to velocity, 4-corner causal highpass at
  1 Hz applied before both the waveform and the spectrogram
- Spectrogram axis 1-45 Hz, -190 to -110 dB
- Waveform limit ±5 × 10⁻⁶ m/s: the strongest bursts in (a) clip, on purpose
- Wind axis 0-6 m/s, the figure's maximum 1-min mean rounded up
- Orange marks are ground transients from `score_transients.py`: triggers on
  HV.MITD and on MITDC or MITDW within 3 s. STA/LTA was run on 2-20 Hz data,
  independently of the 1 Hz display filter. They are drawn on the waveforms
  only. MITDU detected 7 of 13 in the rising window, 5 of 12 in the high
  wind window
- Page fit: 10.5 × 7 in, reduces to 7.5 in wide on US letter with 7-8 pt text

### PPSD figures

From the ppsd working directory, using the `psd_viz` environment:

```bash
P=/Users/prof/micromamba/envs/psd_viz/bin/python
$P run_psd_viz.py
$P plot_permanent_relative.py   # FIGURES/permanent_relative.png
$P plot_covered_relative.py     # FIGURES/covered_relative.png
$P plot_ppsd_composite.py       # FIGURES/ppsd_composite_Z.png
$P plot_wind_summary.py         # FIGURES/wind_summary.png
$P plot_wind.py                 # FIGURES/WIND/ALL_windrose.png etc.
```

`run_psd_viz.py` needs the `compare-mixed-grids` branch of `psd_viz`. The
saved `ppsd_composite_Z.png` was made with a wider frequency axis (about
0.005-40 Hz) than the script's 0.1-10 Hz default; the exact options were not
recorded, so set `--fmin` and `--fmax` to match if you want that version.

## Scripts: PPSD pipeline

All in `ppsd/`.

### `get_bb_data.py`

Downloads the permanent station data from EarthScope. Options are edited in
the file: `choice` selects HVO or the AEC phases. Not a command-line tool.

### `run_psd_viz.py`

Batch driver for `psd_tools.py` from the `psd_viz` repo. Scans `DATA` and
runs three steps: `compute` (one PPSD `.npz` per channel in `SPECTRA`),
`pdf` (a PPSD figure per channel) and `compare` (one figure per main station
and component). A JSON state file in `logs/` records each task, so re-runs
skip what succeeded and also what failed.

- Use `--retry-failed` or `--force` to redo work, `--dry-run` to preview
- Paths to the interpreter and `psd_tools.py` are hard-coded constants
  (`PYTHON`, `SCRIPT`) and can be overridden with `--python` and `--script`
- `old/run_psd_viz.sh` in the ppsd working directory is the superseded
  shell version and was not moved into this repo

### `plot_permanent_relative.py`

Noise of each node relative to its permanent station. Rows are sites,
columns are components. Each station's hourly PSD is smoothed to 1/12-octave
bins, hours recorded by both stations are paired, and the panel shows the
median difference in dB with an interquartile band. All stations at a site
are trimmed to their common window. KAED and RND are excluded by default.

Other scripts import its loading and binning code, so changes here affect
`plot_covered_relative.py`, `plot_ppsd_composite.py` and
`plot_wind_summary.py`.

### `plot_covered_relative.py`

Thin wrapper around the above with `--reference C --omit permanent --sites
KAED,RND`. Used because HV.KAED is unusable and AK.RND shows a shared offset.

### `plot_ppsd_composite.py`

Probabilistic PSD for one component at every station: rows are sites,
columns are Permanent, Covered, Uncovered and Windscreen, with the Peterson
noise models overlaid. `--component`, `--cmap`, `--pmax` control the look.

### `plot_wind.py`

Wind roses and wind-speed timelines from the anemometers, plus a CSV of
daily statistics. Averages 1 Hz samples in 10-min windows, with direction
averaged as a vector. **Screens anemometer faults**: speed below 0 or at
least 100 m/s, or direction outside 0-360, and drops the whole averaging
window around a fault, since neighboring "valid" samples are also corrupt.

### `plot_wind_summary.py`

Wind-induced noise: each station's median spectrum in windy hours (at least
4 m/s) minus its own median in calm hours (below 1 m/s). Only DESD, KAED, MITD
and K24K ever got windy enough. Below about 0.4 Hz the difference also
reflects when windy hours happened, e.g. microseism changes.

## Scripts: wind-state pipeline

All in `windygrams/`. The order of use was: find windows, score them, then
plot.

### `smooth_wind.py`

Writes a centered running mean of the 1 Hz wind speed to `WIND/`, for
overlaying on waveforms. Default 60 s. Edges are renormalized. Called
automatically by the plotting scripts when a file is missing. Files in
`WIND/` are reused, so delete them if the smoothing should change.

### `find_windows.py`

Ranks 4 h windows at each site from 10-min mean wind speed:

- high: largest minimum speed (wind that stays strong)
- low: smallest maximum speed (calm throughout)
- rising or falling: largest difference between the first and last hour

Windows must fit inside one UTC day and have full data coverage. 10-min bins
holding anemometer faults are dropped (same fault rule as `plot_wind.py`).
Its docstring refers to a `WINDOWS` dict, but windows now live in
`configs/windygrams.yaml` under `windows:`.

### `windygrams.py` and `configs/`

Exploratory figures for the windows in a YAML config: waveform, spectrogram
and wind for every chosen window, station and component, using PySEP's
`PrettyPlot`. Output goes to `FIGURES/<wind condition>/<window>/`.

```bash
python windygrams.py                     # configs/windygrams.yaml
python windygrams.py k24k_uz.yaml        # a smaller selection
python windygrams.py k24k_uz.yaml --dry_run  # print the pplot commands
```

Config sections are `paths`, `select`, `processing`, `pplot`, `sites` and
`windows`. A config can start with `base: other.yaml` and list only what
changes; nested sections merge key by key and lists are replaced. Quote
times in YAML (`"13:10"`), since unquoted ones are read as numbers.

- `windygrams.yaml`: all windows, stations and components
- `k24k_uz.yaml`: K24K, Uncovered and Permanent, Z, with envelope
- `soi.yaml`: >1 Hz, -190 to -140 dB view, used for the ~35 Hz DESD line
- `hv_desd_rescale.yaml`: DESD low-wind window at the same scale; linked into
  `CHECKS/desd_35hz/` because its `base` path is relative to that folder
- `checks/checkdesdu.sh`: one pplot command saved from `--dry_run` for DESDU

### `score_states.py`

Scores the top 6 candidates per site and wind state by how well the
5 Hz highpassed Uncovered amplitude tracks the wind (Spearman rho, amplitude
ratio, fraction of transient minutes). Used to choose `STATES`. Writes
`CHECKS/wind_states/scored_windows.csv`.

### `score_transients.py`

Ranks the top 15 candidates per site and state by their impulsive
transients, and measures detection rates.

1. Bandpass the vertical component 2-20 Hz at C, U, W and P
2. Recursive STA/LTA, STA 1 s, LTA 30 s, trigger on 5, off 1.5
3. Cluster triggers across sensors within 3 s
4. A ground transient is a cluster with P and C or W
5. For the detection rate of X, the condition is P plus another node besides
   X, so each sensor is judged against independent ones

Writes to `CHECKS/transients/`, whose two CSVs are symlinks into
`results/`:

- `ranked_windows.csv`: statistics of 477 windows; column meanings are in
  the script's docstring
- `event_times.csv`: transient onsets, read by `plot_wind_states.py`

**The figures with orange marks need `event_times.csv`.**

The four windows of `STATES_MITD_TRANSIENTS` were chosen by hand from this
ranking. Windows with the most transients were favored, except:

- Rising, day 100: chosen for stronger wind and rho 0.91 over day 87, which
  has more transients but less wind
- High, day 95: windows with at least 6 m/s have only 1-2 transients
- Low, day 117: day 86 had the most transients, but its 1-10 Hz noise steps
  down near 05:13 UTC at steady wind, so something else is affecting it

### `plot_wind_states.py`

Makes the paper figures. `PanelPlot` subclasses PySEP's `PrettyPlot` so each
panel draws into axes of a larger figure.

- `setup_plot`: uses the given axes
- `plot_waveforms`: optional extra highpass for the plotted waveform only
  (`wf_fmin`), so the spectrogram can keep a different band
- `plot_tmarks`: faint marks on the waveform axis only
- `finalize`: skipped; `make_figure` saves
- `make_figure`: lays out the grid. Its `opts` control columns, limits,
  ticks, spacing, the wind scale and dated time labels
- `style_panel`: shared axes, labels, and hides spectrogram grid lines,
  which could be mistaken for spectral lines
- `*_panels` functions: build the panel lists

Wind figures use a shared axis (0-14 m/s) unless `wind_auto` is set. Useful
options: `--component`, `--fmin`, `--wf_fmin`, `--sp_fmin`, `--sp_linear`,
`--ylim_states`, `--ylim_mitd`, `--ylim_marked`, `--ylim_covers`, `--sp_vmin`,
`--sp_vmax`, `--path_out`. The 1 Hz figures set their filter in `main`, so
`--fmin` does not change them.

### `longterm_psd.py`

Weeks-long spectrograms to see how common a narrowband line is. Welch PSDs
of 10-min segments, cached as `.npz`, converted with the scalar sensitivity.
Makes a frequency-time figure and a median-by-local-hour figure. Outputs and
caches are in `CHECKS/desd_35hz/` and `CHECKS/narrowband/`.

## Dependencies outside this repo

Both are saved as patches in `external/` so they survive if the other repos
change.

### PySEP `pplot.py`: branch `pplot-envelope-spectrogram`

Committed as `70f6f2d` on branch `pplot-envelope-spectrogram` of
adjtomo/pysep (pushed, not merged), on top of `5441bc7`. The changes add
`sp_vmin`, `sp_interp`, `tr_ylim`, `tr_front`, `wf_env`, `wf_env_smooth`,
`wf_env_color`, and move the title above the spectrogram. **The wind-state
figures fail without them.** Check out that branch in
`/Users/prof/Repos/pysep` (or `master` once it is merged). The same changes
are kept as a patch against `5441bc7`, as a fallback:

```bash
git -C /path/to/pysep apply $PWD/external/pysep_pplot.patch
```

### `psd_viz`: local branch `compare-mixed-grids`

`/Users/prof/repos/psd_viz`, not merged or pushed. Two commits let `compare`
use different frequency grids (AK BH against XX DH/EH/HH) and add `--domain`.
`run_psd_viz.py` comparisons fail for AK sites on `main`. The patch holds only
those commits. The working tree also has personal plot edits (marked BCBC) that
were deliberately kept out of commits.

## Caveats

**Processing**

- `smooth_wind.py` does **not** screen anemometer faults, unlike
  `plot_wind.py` and `find_windows.py`. A fault value such as 1000 m/s would
  be smeared across 60 s. None appear in the paper windows, but check the
  wind trace of any new window. PAX has mixed data types on some days, which
  makes `smooth_wind.py` fail to merge; the scoring scripts cast to float
- The highpass is 4-corner and causal (`zerophase=False`), so a 1 Hz
  corner only reduces power by 3 dB at 1 Hz. The band just above 1 Hz stays
  bright because the microseism extends to about 2 Hz
- ObsPy warns that the node response unit `MV` is unknown and uses the
  response as is. Node and permanent spectra agree in quiet conditions, so
  amplitudes look right, but this is worth knowing if amplitudes are
  questioned
- Waveform limits clip the largest signals on purpose in several figures
- Spectrograms use automatic windows, 81.92 s for 4 h at 100 Hz

**Interpreting the transient analysis**

- Detection rates depend on the STA/LTA settings. These are fixed choices
  and were not tuned or tested for sensitivity
- Marks require a trigger on the permanent station, so they look best there.
  Use the rates, not the plots, to compare sensors
- "High wind" is relative to each site: about 7 m/s for K24K candidates but
  under 2 m/s at PAX
- Fewer transients in high wind (median 4 against 9 per window at MITD) may
  partly reflect time of day, not only missed detections
- Windows were selected for having many transients, so they are examples,
  not an unbiased sample

**Choosing or using windows**

- At 0.1 Hz the microseism controls waveform amplitude at MITD and KAED; the
  wind shows mainly in the spectrograms there. K24K shows it more directly
- A window must fit in one UTC day, because plotting reads a single day file
- `CHECKS/`, `FIGURES/` and `WIND/` are working directory outputs, not in git

## Results worth remembering

**Wind and detection.** Share of ground transients detected over all sites:

| Wind state | Uncovered | Covered | Windscreen |
|------------|-----------|---------|------------|
| Low | 79% | 87% | 84% |
| Rising | 56% | 82% | 83% |
| Falling | 56% | 81% | 85% |
| High | 39% | 72% | 81% |

At MITD the Uncovered node goes from 88% (low) to 38% (high), and the
windscreen node from 89% to 53%, while Covered stays at 93-99%. Uncovered-only
triggers at MITD rise from about 0.4 per 4 h in low wind to about 120 in high.

**Noise in the windiest hour.** Median 5-20 Hz velocity (m/s):

| Sensor | K24K, day 261, 17:10-18:10 | MITD, day 96, 19:00-20:00 |
|--------|------|------|
| Uncovered | 2.97e-7 | 1.81e-7 |
| Covered | 8.86e-8 | 1.09e-7 |
| Windscreen | 6.18e-8 | 1.74e-7 |
| Permanent | 6.10e-8 | 8.39e-8 |

**The windscreen does not behave the same everywhere.** At K24K it matches
the permanent station. At MITD it is nearly as noisy as the uncovered node.

### Narrowband lines (unresolved)

Thin lines between about 20 and 45 Hz drift nightly and over weeks. Clear at
DESD, KAED, MITD and MDM; faint at K24K and PAX; unclear at RND and H23K. They
are absent on HV permanent stations, and AK stations cannot test them above
25 Hz. They were first noticed as a ~35 Hz line on DESDU.

The leading idea is a temperature-dependent resonance of the node-ground
coupling. Perforated windscreens may suppress it at DESD and KAED, but MITDW
is a counterexample. The nodes have no external electrical connection to
blame. Temperature data were to come from colleagues.

Separate from these, steady lines at about 15-20 Hz appear on HV.MITD, and an
11.5 Hz line appears on K24KC from about 15:15 UTC on day 261, with about
9 and 11.5 Hz lines on AK.K24K. None of these has been explained.

## Reproducibility record

Checked in November 2026 after reorganizing, through the new symlinks:

- All 12 figures from `plot_wind_states.py`: pixel-identical to the existing
  files, including both `FOR_PAPER` figures
- `event_times.csv` from `score_transients.py`: identical
- `score_states.py`: matches the earlier scores to within CSV rounding
- `permanent_relative.png`, `covered_relative.png`, `wind_summary.png`:
  pixel-identical
- `ALL_windrose.png`, `ALL_wind_speed.png` are **not** pixel-identical to the
  saved versions (2-4% of pixels differ). The script is unchanged by the
  reorganization; the cause was not investigated
- `ppsd_composite_Z.png`: see the note on its frequency range above
- `run_psd_viz.py`, `get_bb_data.py`, `longterm_psd.py` and full
  `windygrams.py` runs were not repeated. `--help` and `--dry_run` work
