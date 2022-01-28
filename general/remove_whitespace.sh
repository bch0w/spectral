for FIG in "$@"
do
	convert ${FIG} -fuzz 1% -trim +repage ${FIG}
done
