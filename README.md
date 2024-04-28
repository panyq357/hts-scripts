This directory contains some stand-alone scripts which might be useful.

- `check_md5.py`: find all md5 files in current directory recursively, check them and aggregate results.
- `genbank_to_fasta_gff.py`: convert a Genbank file to FASTA and GFF files.
- `make_consensus.bash`: make a new reference genome using ref FASTA and VCF.
- `make_samples.py`: find all FASTQ files in current directory recursively, make a samples.yaml suitable for reads2bwa.smk pipeline.
- `minimap2_homo_fetch.bash`: fetch sequence from ref FASTA, use minimap2 to find homologous sequence in target FASTA.
- `gsea_io.R`: some R functions to convert DESeq2 counts and coldata to GCT and CLS files that may suitable for GSEA (not tested).
- `rename_bigwig_chrom.py`: rename chromesome name in bigwig files.

