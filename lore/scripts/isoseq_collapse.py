import os
import json


header = True
with open(snakemake.output[0], "w") as f1:
    for n in range(len(snakemake.input)):
        sample = os.path.basename(snakemake.input[n]).rsplit(".", 3)[0]
        data = {"sample": sample}

        with open(
            snakemake.input[n]
        ) as f2:  # f'{config["isoseq_correct_dir"]}/{sample}.unsorted.report.json'
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = str(dct["value"])
            if key != "sample_name":
                data[key] = value

        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False
        f1.write("\t".join(data.values()) + "\n")
