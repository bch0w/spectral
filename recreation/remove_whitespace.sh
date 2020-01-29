for FIG in *".png" 
do
	convert ${FIG} -fuzz 1% -trim +repage ${FIG}
done
