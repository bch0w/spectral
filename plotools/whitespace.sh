for FIG in "$@"
do
	magick ${FIG} -fuzz 1% -trim +repage trim_${FIG}
done
