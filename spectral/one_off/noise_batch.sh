cd /seis/prj/fwi/bchow/spectral
#for station in TSZ MRZ BFZ WAZ KHEZ VRZ HIZ TLZ TOZ MKAZ KUZ GRZ WSRZ RATZ WCZ OUZ MWZ BKZ KNZ RTZ HAZ PUZ OPRZ URZ PXZ MXZ
for station in TSZ MRZ BFZ WAZ KHEZ VRZ HIZ TLZ TOZ MKAZ KUZ GRZ WSRZ RATZ WCZ OUZ MWZ BKZ KNZ RTZ HAZ OPRZ URZ PXZ MXZ
do
	python noise_analysis.py --station "${station}" --channel HHZ --start 2017-11-01 --end 2018-02-01 --dec 5
done
