import os
import json


header = True
with open(snakemake.output[0], "w") as f1:
    for n in range(len(snakemake.input)):
        sample = os.path.basename(snakemake.input[0]).rsplit(".", 2)[0]
        data = {
            "sample": sample,
            "edit_read_count": None,
            "unchanged_read_count": None,
        }
        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False

        with open(snakemake.input[n]) as f2:  # f'{config["isoseq_correct_dir"]}/{sample}.report.json'
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = str(dct["value"])
            if key in data:
                data[key] = value
        
        f1.write("\t".join(data.values()) + "\n")
