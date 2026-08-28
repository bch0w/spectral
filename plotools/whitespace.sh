for FIG in "$@"
do
	magick convert ${FIG} -fuzz 1% -trim +repage trim_${FIG}
done
