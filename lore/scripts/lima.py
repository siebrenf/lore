import os


header = True
with open(snakemake.output[0], "w") as f1:
    for n in range(len(snakemake.input)):
        sample = os.path.basename(snakemake.input[n]).rsplit(".", 2)[0]
        data = {"sample": sample}

        prefix = ""
        with open(
            snakemake.input[n]
        ) as f2:  # f'{config["skera_dir"]}/{sample}.lima.summary'
            for line in f2:
                if line == "\n":
                    continue  # empty line
                if line.endswith(":\n"):
                    prefix = line[-6:-2] + " "
                    continue  # header line

                key, value = line.split(": ")
                key = prefix + key.strip()
                value = value.strip()
                data[key] = value

        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False
        f1.write("\t".join(data.values()) + "\n")
