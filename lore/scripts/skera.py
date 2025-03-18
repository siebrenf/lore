import os
import json


# summary
files = [
    "/bank/experiments/2025-03-lore/skera/segmented.summary.json",
]
# summary_out = "/mnt/cn/vol/mbdata/siebrenf/lore/multiqc_report_data/multiqc_skera_summary2.txt"
summary_out = "/bank/experiments/2025-03-lore/qc/multiqc_skera_summary.txt"
header = [
    "sample",
    "m_hifi_reads",     # for the summary table
    "m_s-reads",        # for the summary table
    "avg_len_s-reads",  # for the skera section
    "pct_full_array",     # for the skera section
    "avg_array_size",   # for the skera section
]
with open(summary_out, "w") as f1:
    f1.write("\t".join(header) + "\n")
    for file_name in files:
        sample = os.path.basename(file_name).rsplit(".", 2)[0]
        data = [sample]
        with open(file_name) as f2:
            summary = json.load(f2)
        for md in summary["attributes"]:
            if md["id"] in ["reads", "s_reads"]:
                data.append(str(round(int(md["value"])/1_000_000, 2)))
            elif md["id"] in ["mean_len_s_reads"]:
                data.append(str(md["value"]))
            else:
                data.append(str(round(float(md["value"]), 2)))

        f1.write("\t".join(data) + "\n")
