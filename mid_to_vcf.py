import pandas as pd

config = {
    "mid_in_range": snakemake.input["mid_in_range"],
    "vcf_in_range": snakemake.output["vcf_in_range"]
}


def main():

    mid = pd.read_table(config["mid_in_range"])

    vcf_row_list = []
    for i, row in mid.iterrows():
        ref = row["RefBase"]
        alt = row["SnpBase"]
        vcf_row = row[8:].map(lambda x: base_to_zero_one(x, ref, alt))
        vcf_row = pd.concat([
            pd.Series({"#CHROM": row["Chromosome"], "POS": row["Position"], "ID": f"{row['Chromosome']}-{row['Position']}", "REF": ref, "ALT": alt, "QUAL": ".", "FILTER": ".", "INFO": ".", "FORMAT": "GQ"}),
            vcf_row
        ])
        vcf_row_list.append(vcf_row)

    vcf = pd.concat(vcf_row_list, axis=1).T

    vcf.to_csv(config["vcf_in_range"], index=False, sep="\t")


def base_to_zero_one(base, ref, alt):
    """
    e.g. Convert "-", "G", "A" to "./.", "0/0", "1/1".
    """

    if base == ref:
        return "0/0"
    elif base == alt:
        return "1/1"
    else:
        return "./."

if __name__ == "__main__":
    main()
