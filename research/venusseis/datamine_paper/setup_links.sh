#!/usr/bin/env bash
#
# Link the repo's scripts and configs into the local working directories.
#
# The scripts run from a working directory that holds the (untracked) data
# links and outputs: each script changes into the directory of its symlink,
# so run them through these links and not from the repo. Re-run this after
# moving or adding files. It only replaces symlinks; a real file at a link
# path is left alone and reported.
#
# Usage: ./setup_links.sh [path to the akvenus directory]

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AKVENUS="${1:-/Users/prof/Work/akvenus}"
PPSD="$AKVENUS/ppsd"
WINDY="$AKVENUS/windygrams"

# Make a symlink, refusing to overwrite a real file or directory
link() {
    local target="$1" path="$2"
    if [ -e "$path" ] && [ ! -L "$path" ]; then
        echo "SKIPPED (real file exists): $path"
        return
    fi
    ln -sfn "$target" "$path"
    echo "$path -> $target"
}

mkdir -p "$PPSD" "$WINDY/CHECKS/transients"

for dir in "$PPSD" "$WINDY"; do
    link ../MSEEDS "$dir/DATA"
    link ../STATIONXML "$dir/XML"
done

for fid in "$REPO"/ppsd/*.py "$REPO"/ppsd/FIGURE_PLAN.md; do
    link "$fid" "$PPSD/$(basename "$fid")"
done

for fid in "$REPO"/windygrams/*.py "$REPO"/windygrams/configs/*.yaml \
        "$REPO"/windygrams/checks/*.sh; do
    [ "$(basename "$fid")" = hv_desd_rescale.yaml ] && continue
    link "$fid" "$WINDY/$(basename "$fid")"
done

# This config names its base relative to the CHECKS folder it is run from
mkdir -p "$WINDY/CHECKS/desd_35hz"
link "$REPO/windygrams/configs/hv_desd_rescale.yaml" \
    "$WINDY/CHECKS/desd_35hz/hv_desd_rescale.yaml"

for fid in "$REPO"/windygrams/results/*.csv; do
    link "$fid" "$WINDY/CHECKS/transients/$(basename "$fid")"
done
