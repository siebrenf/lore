import os
import json

all_files = [
    "/bank/experiments/2025-03-lore/isoseq_refine/segmented.filter_summary.report.json",
    "/bank/experiments/2025-03-lore/isoseq_refine/segmented.report.csv",  # TODO: different section?

    "/bank/experiments/2025-03-lore/isoseq_correct/segmented.report.json",

    "/bank/experiments/2025-03-lore/isoseq_collapse/segmented.unsorted.report.json",  # TODO: different section?
]


files = [
    "/bank/experiments/2025-03-lore/isoseq_refine/segmented.filter_summary.report.json",
    "/bank/experiments/2025-03-lore/isoseq_correct/segmented.report.json",
    "/bank/experiments/2025-03-lore/isoseq_collapse/segmented.unsorted.report.json",
]
summary_out = "/bank/experiments/2025-03-lore/qc/multiqc_isoseq.tsv"
header = ["sample"]
data = ["segmented"]
with open(summary_out, "w") as f1:
    for file_name in files:
        with open(file_name) as f2:
            contents = json.load(f2)
        for md in contents["attributes"]:
            k = str(md["id"])
            if k in [
                "sample_name",  # isoseq_refine + isoseq_collapse
                "processed_chunks",  # isoseq_correct (technical stat)
                "processed_reads",  # isoseq_correct (technical stat)
                "processed_bases",  # isoseq_correct (technical stat)
                "filtered_reads",  # isoseq_correct (we dont filter)

                "edit_read_count",  # isoseq_correct (breakdown stat)  # TODO: these could be used for a barplot section
                "unchanged_read_count",  # isoseq_correct (breakdown stat
                "missing_read_count",  # isoseq_correct (breakdown stat
                "found_but_failing_read_count",  # isoseq_correct (breakdown stat
                "failing_read_count",  # isoseq_correct (breakdown stat
            ]:
                continue
            header.append(k)
            data.append(str(md["value"]))
    f1.write("\t".join(header) + "\n")
    f1.write("\t".join(data) + "\n")
