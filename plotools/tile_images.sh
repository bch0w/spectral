# sh tile_images.sh <ncol> <nrow> <files>
FID=${3}
echo "montage -mode concatenate -tile ${1}x${2} ${FID} tiled.png"
