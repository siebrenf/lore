import os
import json


header = True
with open(snakemake.output[0], "w") as f1:
    for n in range(len(snakemake.input)):
        sample = os.path.basename(snakemake.input[n]).rsplit(".", 2)[0]
        data = {"sample": sample}

        with open(
            snakemake.input[n]
        ) as f2:  # f'{config["isoseq_correct_dir"]}/{sample}.bcstats.json'
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = f'{dct["value"]:.2f}'
            data[key] = value

        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False
        f1.write("\t".join(data.values()) + "\n")
