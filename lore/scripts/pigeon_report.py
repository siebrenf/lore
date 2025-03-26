import os
import pandas as pd

header = True
with open(snakemake.output[0], "w") as f1:
    for fname in snakemake.input:
        sample = os.path.basename(fname).rsplit("_saturation.txt", 1)[0]
        df1 = pd.read_table(fname, index_col=None)
        if header:
            row = [""] + df1["reads"].astype(str).to_list()
            f1.write("\t".join(row) + "\n")
            header = False
        for col in [
            "unique_genes",
            "unique_isoforms",
            "unique_genes_known",
            "unique_isoforms_known",
        ]:
            row = [f"{sample} {col}"] + df1[col].astype(str).to_list()
            f1.write("\t".join(row) + "\n")
# df1["m_reads"] = df1["yield_count"] / 1_000_000
# df1["%_reads"] = df1["yield_fraction"] * 100
# df = pd.concat([df, df1[["m_reads", "%_reads"]]])

# df.to_csv(snakemake.output.stats[0], sep="\t", index=True)
