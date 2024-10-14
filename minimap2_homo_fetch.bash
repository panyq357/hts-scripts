#!/usr/bin/env bash

# Given a position range in ref genome, use minimap2 to fetch homologous sequence in target genome.

# Usage: ./minimap2_homo_fetch.bash input.fa target.genome.fa out.fa

minimap2 -a $2 $1 | samtools view -b | bedtools bamtobed -i stdin | bedtools getfasta -fi $2 -bed stdin > $3

