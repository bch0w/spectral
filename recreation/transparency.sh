for FIG in *".png" 
do
	convert ${FIG} -fuzz 10% -transparent white ${FIG}
done
