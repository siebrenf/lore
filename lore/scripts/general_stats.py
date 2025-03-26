import pandas as pd

df = pd.DataFrame()

# isoseq_correct: yield_count & yield_fraction
df1 = pd.read_table(snakemake.input.isoseq_correct_stats[0], index_col=0)
df1["m_reads"] = df1["yield_count"] / 1_000_000
df1["%_reads"] = df1["yield_fraction"] * 100
df = pd.concat([df, df1[["m_reads", "%_reads"]]])

df.to_csv(snakemake.output.stats[0], sep="\t", index=True)
