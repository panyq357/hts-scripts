import argparse
import datetime
import re
import hashlib
import pathlib

parser = argparse.ArgumentParser(description='Find all md5 files and check them all.')

parser.add_argument('--glob', type=str, dest="glob", default="**/md5.txt", help='md5 file glob')
parser.add_argument('-o', type=str, dest="out", default=f"check_md5.{datetime.date.today()}.tsv", help='path to output TSV')

args = parser.parse_args()

PATTERN = re.compile("(\\S+)\\s+(\\S+)")

open(args.out, "wt").close()

for md5_file in pathlib.Path(".").glob(args.glob):

    with open(md5_file, "rt") as f:
        lines = f.readlines()

    for line in lines:

        file_md5, file_name = PATTERN.match(line).groups()
        full_name = (md5_file.parent / file_name).absolute()

        with open(full_name, "rb") as f:
            if file_md5 == hashlib.md5(f.read()).hexdigest():
                status = "OK"
            else:
                status = "FAILED"

        with open(args.out, "at") as f:
            f.write(f"{file_name}\t{file_md5}\t{md5_file}\t{status}\t{full_name}\n")

