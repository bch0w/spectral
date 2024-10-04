for FIG in "$@"
do
	convert ${FIG} -fuzz 1% -trim +repage trim_${FIG}
done
