# QUERY FITS API TO RETRIEVE GEONET cGPS DATA

for station in MNHR DNVK PORA TURI WPUK OROA WPAW PNUI KAHU RAKW KERE KAWK MCNL CKID
do	
	for comp in u #e n 
	do
		curl -H "Accept: application/vnd.geo+json;version=1" "http://fits.geonet.org.nz/observation?typeID=${comp}&siteID=${station}" -o ${station}_${comp}.txt
	done
done

