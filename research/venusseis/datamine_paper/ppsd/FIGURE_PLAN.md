# Figure plan: data mine paper

Goal: a few figures summarizing noise levels per site, highlighting how
co-located stations (e.g., MDM vs MDMC/MDMU/MDMW) differ.

## Main text

1. **[DRAFT DONE] (2) Spectra relative to the permanent station**
   - Script: `plot_permanent_relative.py` (one script per figure)
   - Rows = sites, columns = components (`--components Z,H`, where H is
     the power mean of N and E; `Z,N,E` for all three)
   - Curve = median of paired hourly differences (XX minus permanent, dB,
     on common 1/12-octave bins), 25th-75th percentile band; all stations
     trimmed to the site's common window
   - Z, N and E always shown as separate columns
   - Writes one figure, `FIGURES/permanent_relative.png`: rows DESD, MITD
     (HV), then H23K, K24K, MDM, PAX (AK); KAED and RND left out by
     default (`--exclude-sites`)
   - Legend: Permanent, Covered (C), Uncovered (U), Windscreen (W)
   - Keep the low-frequency end: the XX sensors are 2 s sensors, so the
     rise below ~0.2 Hz is out-of-band and worth showing
   - **KAED + RND figure:** `plot_covered_relative.py` ->
     `FIGURES/covered_relative.png`; same method, Uncovered and
     Windscreen relative to the Covered (C) station, permanent station
     left out (shared code: `plot_permanent_relative.py --reference C
     --omit permanent`)
   - **Checked 2026-11-03: XX.KAEDW HHE** sits ~8 dB below KAEDC HHE.
     Not a metadata error: all 9 KAED XX channels have identical
     responses in `XML/ISHTAR.xml` (same sensitivity, poles/zeros,
     stages, one epoch, E at 90 deg). The offset is in the raw data, is
     frequency dependent (about -10.6 dB at 0.2-1 Hz, -1.5 dB at 5-10 Hz,
     -6 to -9 dB at 20-40 Hz) and shrinks over the deployment (about -9 dB
     in week 1 to about -5 dB by the end). N and Z match KAEDC within
     0.5 dB. Below ~1 Hz the stations see the same ground motion, so this
     points at the KAEDW E sensor component itself (e.g. leveling, mass
     position or damping), not site noise. **Decision: keep it in the
     figure and note it in the caption/text.** Still worth checking the
     field notes for KAEDW
   - Sites flagged for a separate special-analysis figure:
     - **KAED:** reference EHZ looks broken; no N/E reference
     - **RND:** all three XX stations rise together above ~2 Hz, to about
       +30 dB at 10 Hz on every component, which points at the reference
       (AK.RND) rather than the XX installs
     - **H23K (minor):** narrow peaks near 2.7 Hz on all components and
       around 6 Hz on E

2. ~~(1) Band-averaged heatmap~~ **dropped 2026-11-03** (tried, didn't
   like the representation; script deleted)

3. **[NOW] Composite PPSD, vertical component, all stations**
   - Script: `plot_ppsd_composite.py` -> `FIGURES/ppsd_composite_Z.png`
     (`--component N`/`E` for the others)
   - Rows = all 8 sites incl. KAED and RND (same start-date order as
     figure 2); columns = Permanent, Covered, Uncovered, Windscreen
   - Each panel: probability of hourly acceleration PSDs (1/12-octave
     bins x 1 dB bins) over the site's common window, median line,
     Peterson NLNM/NHNM dashed; shared color scale (`--pmax`)
   - Shows HV.KAED.EHZ's broken response directly, and that RND's XX
     stations (all three) rise to about -105 dB at 10 Hz while AK.RND
     looks like the other AK stations

4. **[NOW] (4) Noise vs wind speed (Seppi et al. 2026 Fig. 9 style)**
   - Hourly PSDs grouped by mean wind speed at the site's anemometer
     (`<SITE>A`, read with `plot_wind.py`); each station's windy median
     spectrum minus its own calm median (decided: own calm spectrum, not
     the permanent station)
   - Per-site multi-bin version (`plot_wind_binned.py`) tried for MITD and
     K24K, then **dropped 2026-11-03** in favor of the summary below
   - **Chosen design:** `plot_wind_summary.py` ->
     `FIGURES/wind_summary.png`; rows = sites, columns = Z/N/E, one curve
     per installation = median spectrum in wind >= 4 m/s minus < 1 m/s
     (`--calm`, `--windy`)
   - Only DESD, KAED, MITD and K24K ever got windy: H23KA, MDMA, PAXA and
     RNDA never exceeded 2.6-2.9 m/s (hourly means), so those sites have
     no windy hours. Worth checking whether that's real (sheltered sites)
     or an anemometer problem
   - Caveat: below ~0.4 Hz the windy-minus-calm difference also reflects
     when windy hours happened (e.g. microseism changes), not wind itself

## Supplement

- **Full-page grid** of PPSD comparisons: every station and every component
  (rows = sites, columns = Z/N/E)

## Decisions and caveats

- **Reference:** the site's permanent station (NLNM tried first and
  dropped). KAED (broken HV.KAED.EHZ, no N/E) and RND (shared offset vs
  AK.RND) use the Covered station instead, in their own figure.
- **Same time window:** each site is trimmed automatically to the dates
  all its stations share (shown in each row of figure 2)
- **Frequency range:** 0.1-10 Hz (AK reference stations only go to 25 Hz)

## Figure folders

- `FIGURES/` top level: paper figures only (e.g. `permanent_relative.png`)
- `FIGURES/ppsd/`: per-station PPSD PDFs (`run_psd_viz.py` pdf step)
- `FIGURES/compare/`: per-site, per-component comparisons (compare step);
  `compare/legacy/` holds old shell-script versions
- `FIGURES/WIND/`, `FIGURES/MITD_collage/`: your own folders, unchanged

## Status

- `run_psd_viz.py`: compute, PDF and compare steps done; full run started
  2026-11-03
- psd_viz fixes on the `compare-mixed-grids` branch (not merged, not pushed)
