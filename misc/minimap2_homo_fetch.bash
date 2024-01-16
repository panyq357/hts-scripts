# Given a position range in ref genome, use minimap2 to fetch homologous sequence in target genome.
#
# author: panyq357
# date: 2023-10-23

out_name="out.fa"      # Prefix of output file.
chr="chr01"            # Chromosome in ref genome.
pos_start=$((100000))  # Position in ref genome.
pos_end=$((101000))

ref_genome="path/to/ref.genome.fa" 
target_genome="path/to/target.genome.fa"

samtools faidx ${ref_genome} ${chr}:${pos_start}-${pos_end} \
| minimap2 -a ${target_genome} - \
| samtools view -b \
| bedtools bamtobed -i stdin \
| bedtools getfasta -fi ${target_genome} -bed stdin > ${out_name}

