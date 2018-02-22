cd ..
#strong motion for kaikoura
#for station in TBAS TUDS TDHS ECLS GKBS GISS WICS WFSS KFHS HNPS THPS ORCS MWDS TTHS WAIS TPPS RPCS TIRS ROTS OPCS KAFS WKHS TKHS TUHS WCDS
#broadband for M5 events
#for station in MXZ HAZ PUZ URZ MWZ RTZ KNZ BKZ PXZ BFZ RATZ TSZ WAZ KHEZ VRZ HIZ TLZ TOZ OPRZ MRZ 
#do
#python eq_duration.py "${station}" 2014p240655
#done
for station in MXZ HAZ PUZ URZ MWZ RTZ KNZ BKZ PXZ BFZ RATZ TSZ WAZ KHEZ VRZ HIZ TLZ TOZ OPRZ MRZ
do                                                                               
python eq_duration.py "${station}" 2015p822263                                   
done  
for station in MXZ HAZ PUZ URZ MWZ RTZ KNZ BKZ PXZ BFZ RATZ TSZ WAZ KHEZ VRZ HIZ TLZ TOZ OPRZ MRZ
do                                                                               
python eq_duration.py "${station}" 2014p240655                                   
done  
for station in MXZ HAZ PUZ URZ MWZ RTZ KNZ BKZ PXZ BFZ RATZ TSZ WAZ KHEZ VRZ HIZ TLZ TOZ OPRZ MRZ
do                                                                               
python eq_duration.py "${station}" 2017p059122
done
