# ref: <https://gist.github.com/dpryan79/39c70b4429dd4559d88fb079b8669721>

from pathlib import Path

import pyBigWig

d = {
    "chr1": "chr01",
    "chr2": "chr02",
    "chr3": "chr03",
    "chr4": "chr04",
    "chr5": "chr05",
    "chr6": "chr06",
    "chr7": "chr07",
    "chr8": "chr08",
    "chr9": "chr09",
    "chr10": "chr10",
    "chr11": "chr11",
    "chr12": "chr12"
}

def main():

    for bw_path in Path("rawdata").iterdir():

        bw = pyBigWig.open(str(bw_path))
        hdr = [(d[chrom], length) for chrom, length in bw.chroms().items() if chrom in d]
        bwOutput = pyBigWig.open("resources/" + bw_path.stem + ".bw", "w")
        bwOutput.addHeader(hdr)
        for chrom, length in bw.chroms().items():
            if chrom in d.keys():
                ints = bw.intervals(chrom, 0, length)
                if len(ints):
                    bwOutput.addEntries([d[chrom]] * len(ints), [x[0] for x in ints], ends=[x[1] for x in ints], values=[x[2] for x in ints])
        bw.close()
        bwOutput.close()

if __name__ == "__main__":
    main()
