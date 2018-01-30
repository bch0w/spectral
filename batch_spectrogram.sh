#!/bin/bash
# shell script, for station list, download station data then create spectrogram
#for station in TBAS GWTS GISS GKBS TUDS WFSS KFHS MWDS ORCS TTHS TPPS WAIS TIRS RPCS ROTS TUHS NSPS HNPS WPWS TDHS
#GHHS NAAS NCDS NGHS NCHS HCDS
for station in KAFS OPCS WKHS MWFS WAKS THHS UTKS MNGS TBCS MMCS HBHS TKHS
do
python getdata_fdsn.py --station "${station}" --channel BN* --start 2016-11-13 --end 2016-11-14 --response True
python make_spectrogram.py "${station}" z
done
