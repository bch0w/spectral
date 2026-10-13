"""
Standalone renaming script to change node naming conventions. Reads waveform
files whose names follow `input_convention`, renames network and station
codes in both the filename and the trace headers, and writes the results to
an output directory named by `output_convention`. Files are processed in
parallel with `concurrent.futures`.
"""
import os
import re
import sys
from string import Formatter
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob
from obspy import read

# USER PARAMETERS
network_dict = {
        "XX": "XX"
        }
station_dict = {
        "101": "HC01",
        "102": "HC02",
        "103": "HC03",
        "104": "HC04",
        "105": "HC05",
        "106": "HC06",
        "107": "HC07",
        "108": "HC08",
        "109": "HC09",
        "110": "HC10",
        "111": "HC11",
        "112": "HC12",
        "113": "HC13",
        "115": "HC15",
        "116": "HC16",
        "201": "JB01",
        "202": "JB02",
        "203": "JB03",
        "204": "JB04",
        "205": "JB05",
        "206": "JB06",
        "207": "JB07",
        "208": "JB08",
        "209": "JB09",
        "210": "JB10",
        "211": "JB11",
        "212": "JB12",
        "213": "JB13",
        }
# Filename templates. Wrap each field name in curly braces; anything outside
# braces is matched literally. `sta` is required in the input; `net` is
# optional. Both are renamed using the dicts above. If `net`, `loc` or `cha`
# are present, they are also written into the trace headers. Any other field
# names (e.g., `typ`, `year`, `jday`) are free-form and are just passed
# through to the output.
# Every field in `output_convention` must also appear in `input_convention`.
# If two fields have no separator between them, give at least one of them a
# fixed character width after a colon, e.g., "{year:4}{jday:3}".
# Examples:
#   SeisComP: "{net}.{sta}.{loc}.{cha}.{typ}.{year}.{jday}"
#   Other:    "{sta}_{year:4}{jday:3}_{cha}.mseed"
input_convention = "{net}.{sta}.{loc}.{cha}.{typ}.{year}.{jday}"
output_convention = "{net}.{sta}.{loc}.{cha}.{typ}.{year}.{jday}"
indir = "./"
outdir = "./"
nproc = os.cpu_count()  # number of parallel worker processes
overwrite = False  # if False, skip files whose output already exists


def template_fields(template):
    """
    List the field names used in a filename template.

    :type template: str
    :param template: filename template, e.g., "{net}.{sta}.{cha}"
    :rtype: list of str
    :return: field names in the order they appear
    """
    return [name for _, name, _, _ in Formatter().parse(template) if name]


def template_to_regex(template):
    """
    Convert a filename template into a compiled regular expression where
    each field becomes a named group and literal text is escaped. Fields
    with a width (e.g., "{year:4}") match exactly that many characters,
    otherwise fields match lazily and may be empty (e.g., a blank location
    code).

    :type template: str
    :param template: filename template, e.g., "{net}.{sta}.{cha}"
    :rtype: re.Pattern
    :return: pattern that matches a full filename against the template
    """
    pattern = ""
    for literal, name, spec, _ in Formatter().parse(template):
        pattern += re.escape(literal)
        if name:
            width = f"{{{spec}}}" if spec else "*?"
            pattern += f"(?P<{name}>.{width})"
    return re.compile(pattern)


def template_to_glob(template):
    """
    Convert a filename template into a glob wildcard by replacing each field
    with `*`.

    :type template: str
    :param template: filename template, e.g., "{net}.{sta}.{cha}"
    :rtype: str
    :return: glob wildcard, e.g., "*.*.*"
    """
    return "".join(glob_escape(literal) + ("*" if name else "")
                   for literal, name, _, _ in Formatter().parse(template))


def glob_escape(text):
    """
    Escape glob special characters so they are matched literally.

    :type text: str
    :param text: literal text from a filename template
    :rtype: str
    :return: escaped text
    """
    return re.sub(r"([*?\[])", r"[\1]", text)


def check_conventions():
    """
    Validate the input and output filename templates before processing.

    :raises ValueError: if a template is malformed, `sta` is missing from
        the input, a field is repeated in the input, a field width is not
        an integer, two adjacent input fields both lack a width, or the
        output uses a field not present in the input
    """
    parsed = list(Formatter().parse(input_convention))
    for i, (literal, name, spec, _) in enumerate(parsed):
        if spec and not spec.isdigit():
            raise ValueError(f"width of '{{{name}:{spec}}}' must be an "
                             f"integer")
        # Two fields with no separator can't be split unless one is fixed
        if i and name and not literal:
            prev_name, prev_spec = parsed[i - 1][1:3]
            if prev_name and not (spec or prev_spec):
                raise ValueError(f"'{{{prev_name}}}{{{name}}}' in "
                                 f"input_convention is ambiguous, give one "
                                 f"a width, e.g., '{{{prev_name}:4}}'")

    in_fields = template_fields(input_convention)
    out_fields = template_fields(output_convention)
    if "sta" not in in_fields:
        raise ValueError("input_convention must contain '{sta}'")
    dupes = {f for f in in_fields if in_fields.count(f) > 1}
    if dupes:
        raise ValueError(f"input_convention repeats field(s): {dupes}")
    missing = set(out_fields) - set(in_fields)
    if missing:
        raise ValueError(f"output_convention field(s) {missing} not in "
                         f"input_convention")


def parse_fid(fid):
    """
    Split a filename into its components according to `input_convention`
    and map network and station codes to their new values.

    :type fid: str
    :param fid: full path to the input waveform file
    :rtype: dict
    :return: filename components with `net` and `sta` already renamed
    :raises ValueError: if the filename does not match `input_convention`
    :raises KeyError: if the network or station code is not in the lookup
        dictionaries
    """
    match = template_to_regex(input_convention).fullmatch(
        os.path.basename(fid))
    if match is None:
        raise ValueError(f"filename does not match input_convention: {fid}")
    codes = match.groupdict()

    if "net" in codes:
        if codes["net"] not in network_dict:
            raise KeyError(f"network '{codes['net']}' not in network_dict "
                           f"({fid})")
        codes["net"] = network_dict[codes["net"]]
    if codes["sta"] not in station_dict:
        raise KeyError(f"station '{codes['sta']}' not in station_dict "
                       f"({fid})")
    codes["sta"] = station_dict[codes["sta"]]
    return codes


def rename(fid):
    """
    Rename a single waveform file, updating trace headers and writing the
    result to `outdir`. Skips the file if the output exists and `overwrite`
    is False.

    :type fid: str
    :param fid: full path to the input waveform file
    :rtype: str
    :return: status message describing what was done
    """
    codes = parse_fid(fid)
    fid_out = os.path.join(outdir, output_convention.format(**codes))

    if os.path.abspath(fid_out) == os.path.abspath(fid):
        return f"SKIP (output identical to input): {fid}"
    if os.path.exists(fid_out) and not overwrite:
        return f"SKIP (exists): {fid_out}"

    st = read(fid)
    for tr in st:  # 1 tr per st expected
        if "net" in codes:
            tr.stats.network = codes["net"]
        tr.stats.station = codes["sta"]
        if "loc" in codes:
            tr.stats.location = codes["loc"]
        if "cha" in codes:
            tr.stats.channel = codes["cha"]

    # Write to a temp file then move, so an interrupted run never leaves a
    # partial file that a later run would mistake as already complete
    fid_tmp = f"{fid_out}.tmp"
    st.write(fid_tmp, format="MSEED")
    os.replace(fid_tmp, fid_out)

    return f"{fid} -> {fid_out}"


def main():
    """
    Collect input files, validate all filenames up front, then rename them
    in parallel.
    """
    try:
        check_conventions()
    except ValueError as e:
        sys.exit(str(e))

    # Only pick up regular files that look like waveform files, which avoids
    # this script, temp files and subdirectories (e.g., `outdir`)
    wildcard = os.path.join(indir, template_to_glob(input_convention))
    fids = sorted(f for f in glob(wildcard)
                  if os.path.isfile(f) and not f.endswith(".tmp"))
    if not fids:
        sys.exit(f"no matching files found in {indir}")

    # Validate before doing any work so a bad dict entry fails fast rather
    # than leaving a partially renamed dataset
    errors = []
    for fid in fids:
        try:
            parse_fid(fid)
        except (KeyError, ValueError) as e:
            errors.append(e.args[0])
    if errors:
        print("\n".join(errors))
        sys.exit(f"{len(errors)} file(s) could not be parsed, exiting")

    os.makedirs(outdir, exist_ok=True)

    nfail = 0
    with ProcessPoolExecutor(max_workers=nproc) as executor:
        futures = {executor.submit(rename, fid): fid for fid in fids}
        for future in as_completed(futures):
            try:
                print(future.result())
            except Exception as e:
                nfail += 1
                print(f"FAIL: {futures[future]}: {e}")

    print(f"finished {len(fids)} files with {nfail} failures")
    if nfail:
        sys.exit(1)


if __name__ == "__main__":
    main()
