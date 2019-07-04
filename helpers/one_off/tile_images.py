import glob
import subprocess


compare_a = "default"
compare_b = "hikurangi_strict"

path_wildcard = "./figures/*_{}*"
pngs = glob.glob(path_wildcard.format(compare_a))
for png_a in pngs:
    png_b = png_a.replace(compare_a, compare_b)
    png_c = png_a.replace(compare_a, "compare")
    subprocess.run(["montage", png_a, png_b, "-tile",
        "2x1", "-geometry", "+0+0", png_c])


