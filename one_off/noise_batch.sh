cd /seis/prj/fwi/bchow/spectral
for station in TSZ MRZ BFZ WAZ KHEZ VRZ HIZ TLZ TOZ MKAZ KUZ GRZ WSRZ RATZ WCZ OUZ MWZ BKZ KNZ RTZ HAZ PUZ OPRZ URZ PXZ MXZ
do
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-01-01 --end 2015-02-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-02-01 --end 2015-03-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-03-01 --end 2015-04-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-04-01 --end 2015-05-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-05-01 --end 2015-06-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-06-01 --end 2015-07-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-07-01 --end 2015-08-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-08-01 --end 2015-09-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-09-01 --end 2015-10-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-10-01 --end 2015-11-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-11-01 --end 2015-12-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHE --start 2015-12-01 --end 2016-01-01 --dec 5
done

for station in TSZ MRZ BFZ WAZ KHEZ VRZ HIZ TLZ TOZ MKAZ KUZ GRZ WSRZ RATZ WCZ OUZ MWZ BKZ KNZ RTZ HAZ PUZ OPRZ URZ PXZ MXZ
do
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-01-01 --end 2015-02-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-02-01 --end 2015-03-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-03-01 --end 2015-04-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-04-01 --end 2015-05-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-05-01 --end 2015-06-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-06-01 --end 2015-07-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-07-01 --end 2015-08-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-08-01 --end 2015-09-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-09-01 --end 2015-10-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-10-01 --end 2015-11-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-11-01 --end 2015-12-01 --dec 5
	python noise_analysis.py --station "${station}" --channel HHN --start 2015-12-01 --end 2016-01-01 --dec 5
done
