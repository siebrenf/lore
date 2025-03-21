import os
import json


header = True
with open(snakemake.output.barcodes[0], "w") as f1, open(
    snakemake.output.stats[0], "w"
) as f2:
    for n in range(len(snakemake.input)):
        sample = os.path.basename(snakemake.input[n]).rsplit(".", 2)[0]
        data1 = {
            "sample": sample,
            "edit_read_count": None,
            "unchanged_read_count": None,
        }
        data2 = {"sample": sample}

        with open(
            snakemake.input[n]
        ) as f3:  # f'{config["isoseq_correct_dir"]}/{sample}.report.json'
            contents = json.load(f3)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = f'{dct["value"]:.2f}'
            if key in data1:
                data1[key] = value
            if key not in [
                "processed_chunks",
                # "processed_bases",
            ]:
                data2[key] = value

        if header:
            f1.write("\t".join(data1.keys()) + "\n")
            f2.write("\t".join(data2.keys()) + "\n")
            header = False
        f1.write("\t".join(data1.values()) + "\n")
        f2.write("\t".join(data2.values()) + "\n")
