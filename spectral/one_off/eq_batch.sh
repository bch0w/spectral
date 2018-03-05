cd ..
#strong motion for kaikoura
#for station in TBAS TUDS TDHS ECLS GKBS GISS WICS WFSS KFHS HNPS THPS ORCS MWDS TTHS WAIS TPPS RPCS TIRS ROTS OPCS KAFS WKHS TKHS TUHS WCDS
#broadband for M5 events
for event in 2922302 
	do
		echo "${event}"
		for station in TRVZ NTVZ MXZ HAZ PUZ URZ MWZ RTZ KNZ BKZ PXZ BFZ RATZ TSZ WAZ KHEZ VRZ HIZ TLZ TOZ OPRZ MRZ 
		do
			echo "${station}"
				python eq_duration_FOR_SSW.py "${station}" "${event}"
		done
	done
