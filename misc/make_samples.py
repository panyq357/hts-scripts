# Generate samples yaml suitable for reads2bwa.smk.
#
# Author: panyq357
# Date: 2024-01-02

from collections import defaultdict
from pathlib import Path
import re

import yaml

config = {
    # Glob pattern of raw FASTQ files.
    "glob": "**/*.gz",

    # Regex for extract group and read orientation info.
    "pattern": r".*/([^_]+)_.*\.(R1|R2).fastq.gz",

    # Output filename.
    "out": "samples.yaml"
}


def main():

    PATTERN = re.compile(config["pattern"])

    samples = defaultdict(lambda: defaultdict(list))
    for file in Path(".").glob(config["glob"]):
        sample, r1r2 = PATTERN.match(str(file)).groups()
        samples[sample][r1r2].append(file.absolute())

    samples = dict(samples)
    for sample in samples:
        samples[sample] = dict(samples[sample])
        if len(samples[sample]["R1"]) != len(samples[sample]["R2"]):
            raise Exception(f"len(R1) != len(R2): {sample}")
        elif len(samples[sample]["R1"]) == 1:
            samples[sample]["R1"] = str(samples[sample]["R1"][0])
            samples[sample]["R2"] = str(samples[sample]["R2"][0])
        else:
            samples[sample]["R1"] = list(map(str, sorted(samples[sample]["R1"])))
            samples[sample]["R2"] = list(map(str, sorted(samples[sample]["R2"])))

    with open(config["out"], "wt") as f:
        f.write(yaml.dump({"samples": samples}, sort_keys=True))

if __name__ == "__main__":
    main()
