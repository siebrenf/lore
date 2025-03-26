import json
import os


header = True
with open(snakemake.output[0], "w") as f1:  # f'{config["qc_dir"]}/reads_per_step.tsv'
    for n in range(len(snakemake.input.skera)):
        sample = os.path.basename(snakemake.input.skera[n]).rsplit(".", 2)[0]
        data = {
            "": sample,
            # skera (deconcatenate reads)
            "ccs_reads": None,
            "s_reads": None,
            # lima (primers)
            "primed_reads": None,  # matching primers + passing minimal thresholds
            # refine (reads)
            "fl_reads": None,
            "flnc_reads": None,
            "polya_reads": None,
            # correct (barcodes)
            "non-missing_reads": None,
            "yield_reads": None,
        }
        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False

        f_skera = snakemake.input.skera[
            n
        ]  # f'{config["skera_dir"]}/{sample}.summary.json'
        with open(f_skera) as f2:
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = dct["id"]
            value = str(dct["value"])
            if key == "reads":
                data["ccs_reads"] = value
            elif key == "s_reads":
                data["s_reads"] = value

        f_lima = snakemake.input.lima[
            n
        ]  # f'{config["skera_dir"]}/{sample}.lima.summary'
        dct = {}
        with open(f_lima) as f2:
            for line in f2:
                if line.startswith("Reads above all thresholds (A) : "):
                    data["primed_reads"] = line.split(": ")[1].strip()
                    break

        f_refine = snakemake.input.refine[
            n
        ]  # f'{config["isoseq_refine_dir"]}/{sample}.filter_summary.report.json'
        with open(f_refine) as f2:
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = str(dct["value"])
            if key == "num_reads_fl":
                data["fl_reads"] = value
            elif key == "num_reads_flnc":
                data["flnc_reads"] = value
            elif key == "num_reads_flnc_polya":
                data["polya_reads"] = value

        f_correct = snakemake.input.correct[
            n
        ]  # f'{config["isoseq_correct_dir"]}/{sample}.report.json'
        remaining = None
        with open(f_correct) as f2:
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = str(dct["id"])
            value = str(dct["value"])
            if key == "processed_reads":
                remaining = dct["value"]
            elif key == "missing_read_count":
                remaining -= dct["value"]
                data["non-missing_reads"] = str(remaining)
            elif key == "yield_count":
                data["yield_reads"] = value

        f1.write("\t".join(data.values()) + "\n")
