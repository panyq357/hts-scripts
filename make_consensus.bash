
# ref: http://samtools.github.io/bcftools/howtos/consensus-sequence.html

ref_genome="/home/panyq/Tools/index-scripts/os/rap-db/2023-11-03/results/genome/os.rap-db.genome.fa"
vcf="results/concat_and_norm/ZH.norm.vcf.gz"
out_fa="results/ZH.consensus.fa"
out_chain="results/ZH.chain.txt"

if [[ ! -e $(dirname out_fa) ]] ; then
    mkdir -p $(dirname out_fa)
fi

bcftools norm --threads 10 --fasta-ref ${ref_genome} ${vcf} \
| bcftools filter --IndelGap 5 \
| bgzip > norm.filter.vcf.gz

tabix norm.filter.vcf.gz
    
bcftools consensus --fasta-ref ${ref_genome} --chain ${out_chain} norm.filter.vcf.gz  > ${out_fa}

rm norm.filter.vcf.gz norm.filter.vcf.gz.tbi

