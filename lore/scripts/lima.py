import os
import pandas as pd


# counts
files = [
    "/bank/experiments/2025-03-lore/lima/segmented.lima.counts",
]
counts_out = "/mnt/cn/vol/mbdata/siebrenf/lore/multiqc_report_data/multiqc_lima_counts2.txt"
header = ["Sample", "M Reads (lima)", "Quality Score (lima)"]
with open(counts_out, "w") as f:
    f.write("\t".join(header)+"\n")
    for file_name in files:
        sample = os.path.basename(file_name).rsplit(".", 2)[0]
        df = pd.read_table(file_name, index_col=None)
        data = [
            sample,
            str(df.at[0, "Counts"]/1_000_000),
            str(df.at[0, "MeanScore"]),
        ]
    f.write("\t".join(data)+"\n")


# summaries
files = [
    "/bank/experiments/2025-03-lore/lima/segmented.lima.summary",
]
counts_out = "/mnt/cn/vol/mbdata/siebrenf/lore/multiqc_report_data/multiqc_lima_summaries2.txt"
header = [
    "Sample",
    "% Failed",
    "% Too short",
    "% Low score",
    "% Low end score",
    "% Low # passes",
    "% Low score lead",
    "% Low ref span",
    "% 5p--5p pairs",
    "% 3p--3p pairs",
]
with open(counts_out, "w") as f1:
    f1.write("\t".join(header)+"\n")
    for file_name in files:
        sample = os.path.basename(file_name).rsplit(".", 2)[0]
        keys = []
        data = [sample]
        total_reads = 0
        with open(file_name) as f2:
            for line in f2:
                line = line.strip()
                if len(line) == 0:
                    continue
                k, v = line.split(":")
                k = k.strip()
                v = v.strip()
                if v == "":
                    pass
                elif k == "Reads input":
                    total_reads = int(v)
                elif k in [
                    "Reads above all thresholds (A)",
                    "Without SMRTbell adapter",
                    "Wrong different pair",
                    "With same pair",
                    "With different pair",
                ]:
                    pass
                else:
                    keys.append(k)
                    data.append(
                        round(100 * (int(v) / total_reads), 2)
                    )
        f1.write("\t".join(data)+"\n")
