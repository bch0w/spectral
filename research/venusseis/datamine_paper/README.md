# ISHTAR data mine paper: seismic noise and wind

Scripts for the data mine paper figures. Run them through the symlinks in
`/Users/prof/Work/akvenus/{ppsd,windygrams}` (recreate with
`./setup_links.sh`). Env `psd_viz` for `ppsd/`, `work` for `windygrams/`.

**In the paper:** `k24k_wind_data.png` (= `wind_covers.png`) and
`wind_compare_mitdu_mitdc_hvmitd_1hz.png`, both from `plot_wind_states.py`
and copied by hand to `windygrams/FIGURES/FOR_PAPER/`.

## ppsd/

- `get_bb_data.py`: downloads permanent-station data from EarthScope.
  Options are edited in the file (`choice`).
- `run_psd_viz.py`: batch driver for `psd_viz`'s `psd_tools.py`; computes
  PPSDs into `SPECTRA/` and makes per-channel and comparison figures.
  Needs the `compare-mixed-grids` branch of `psd_viz`.
- `plot_permanent_relative.py`: median noise of each node relative to its
  permanent station, per site and component → `permanent_relative.png`.
  Other ppsd plot scripts import its loading code.
- `plot_covered_relative.py`: same, relative to the covered node
  → `covered_relative.png`.
- `plot_ppsd_composite.py`: PPSD of every station for one component, with
  the Peterson noise models → `ppsd_composite_Z.png`.
- `plot_wind.py`: wind roses, wind-speed timelines and daily stats; screens
  anemometer faults → `WIND/ALL_windrose.png`, `ALL_wind_speed.png`.
- `plot_wind_summary.py`: median spectrum in windy (≥4 m/s) minus calm
  (<1 m/s) hours per station → `wind_summary.png`.
- `FIGURE_PLAN.md`: status and decisions for the PPSD figures.

## windygrams/

- `smooth_wind.py`: 60 s running mean of 1 Hz wind into `WIND/`; called by
  the plot scripts. Does not screen anemometer faults.
- `find_windows.py`: ranks 4 h windows per site as high, low, rising or
  falling wind from 10-min means.
- `windygrams.py`: exploratory waveform + spectrogram + wind figures for
  the windows in a YAML config (`configs/`) → `FIGURES/<state>/`.
  `--dry_run` prints the pplot commands.
- `score_states.py`: scores candidate windows by how well uncovered-node
  amplitude tracks the wind; used to pick the wind-state windows.
- `score_transients.py`: STA/LTA (2-20 Hz) on all sensors of a site,
  clustered within 3 s; ranks windows by transients and measures detection
  rates → `results/ranked_windows.csv`, `results/event_times.csv`.
- `plot_wind_states.py`: makes all wind-state figures, including both paper
  figures, → `FIGURES/paper/`. Windows are constants at the top; orange
  transient marks come from `event_times.csv`.
- `longterm_psd.py`: weeks-long spectrograms and median-by-hour plots for
  tracking narrowband lines → `CHECKS/narrowband/`, `CHECKS/desd_35hz/`.

## Code outside this repo

- PySEP `pplot.py` changes: branch `pplot-envelope-spectrogram` (pushed,
  not merged). Wind-state figures fail without it.
- `psd_viz`: local-only branch `compare-mixed-grids`, not pushed.
- Both are also saved as patches in `external/`.

## Unsolved

- **Wandering narrowband lines, ~20-45 Hz, on the nodes.** Drift nightly
  and over weeks; clear at DESD, KAED, MITD, MDM; absent on HV permanent
  stations. Best guess: temperature-dependent node-ground coupling
  resonance. Waiting on temperature data. Next: track peak frequency per
  node at night (does it wander more on U than C, W? tied to sunrise?),
  check MITDW and DESDW field notes, correlate with temperature.
- Other unexplained lines: steady ~15-20 Hz on HV.MITD; 11.5 Hz on K24KC
  after 15:15 UTC on day 261; ~9 and 11.5 Hz on AK.K24K.
- Low wind at H23KA, MDMA, PAXA, RNDA (never >2.9 m/s hourly): real or
  anemometer problem?
- KAEDW HHE sits ~8 dB below KAEDC HHE, time- and frequency-dependent;
  check field notes.
- MITDC-only spike near 17:04 UTC in panel (b) of the MITD figure: mention
  or not?
- `ALL_windrose.png` and `ALL_wind_speed.png` no longer reproduce exactly
  (2-4% of pixels differ); cause unknown.
