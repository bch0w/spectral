#!/bin/bash
# shell script, for station list, download station data then create spectrogram
for station in TBAS GWTS GISS GKBS TUDS WFSS KFHS MWDS ORCS TTHS TPPS WAIS TIRS RPCS ROTS TUHS NSPS HNPS WPWS TDHS GHHS NAAS NCDS NGHS NCHS HCDS KAFS OPCS WKHS MWFS WAKS THHS UTKS MNGS TBCS MMCS HBHS TKHS
do
python make_spectrogram.py "${station}" z
done
