import json
import os


header = True
with open(snakemake.output[0], "w") as f1:
    for n in range(len(snakemake.input.raw)):
        sample = os.path.basename(snakemake.input.raw[n]).rsplit(".", 2)[0]
        data = {"sample": sample}

        with open(
            snakemake.input[n]
        ) as f2:  # f'{config["pigeon_classify_dir"]}/{sample}.report.json'
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = f'{dct["value"]:.2f}'
            data[key] = value

        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False
        f1.write("\t".join(data.values()) + "\n")

        # same for the filtered report
        data = {"sample": f"{sample} (filtered)"}
        fname = os.path.join(
            os.path.dirname(snakemake.input[n]),
            f"{sample}_classification.filtered.report.json"
        )
        with open(fname) as f2:
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = f'{dct["value"]:.2f}'
            data[key] = value

        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False
        f1.write("\t".join(data.values()) + "\n")
