from pathlib import Path
from os import chdir
from argparse import ArgumentParser

from Bio import SeqIO

parser = ArgumentParser(description='Convert Genbank file to GFF and FASTA for visualization of features in IGV')

parser.add_argument('--genbank', required=True, type=str, dest="genbank", help="path to Genbank file")
parser.add_argument('--prefix', type=str, dest="prefix", default=None, help='output file prefix')
parser.add_argument('--outdir', type=str, dest="outdir", default=None, help='output directory')

args = parser.parse_args()

records = list(SeqIO.parse(args.genbank, "gb"))

if len(records) > 1:
    print("Multiple records found in Genbank file, only first will be processed.")

record = records[0]

if args.prefix is not None:
    seqname=args.prefix
else:
    try:
        seqname=record.name
    except AttributeError:
        print("No record name found in Genbank file, please provide a prefix.")
        exit(1)

fasta_lines = [f">{seqname}\n", f"{record.seq}\n"]

gff_lines = []
for f in record.features:

    start = f.location.start + 1
    end = f.location.end
    feature = f.type

    if f.location.strand == 1:
        strand = "+"
    elif f.location.strand == -1:
        strand = "-"
    else:
        strand = "."

    if "label" not in f.qualifiers:
        continue
    else:
        label = f.qualifiers["label"][0].replace(" ", "-")  # IGV parse spaces incorrectly.

    gff_lines.append(f"{seqname}\t.\tmisc_feature\t{start}\t{end}\t.\t{strand}\t.\tName=\"{label}\";GenbankType=\"{feature}\"\n")

if args.outdir is not None:
    outdir = Path(args.outdir)
    if not outdir.exists():
        outdir.mkdir(parents=True)
    chdir(args.outdir)

with open(f"{seqname}.fa", "wt") as f:
    f.writelines(fasta_lines)
with open(f"{seqname}.gff", "wt") as f:
    f.writelines(gff_lines)

