import os
import json


header = True
with open(snakemake.output[0], "w") as f1:  # f'{config["qc_dir"]}/reads_per_step.tsv'
    for n in range(len(snakemake.input.skera)):
        sample = os.path.basename(snakemake.input.skera[0]).rsplit(".", 2)[0]
        data = {
            "sample": sample,
            # skera (deconcatenate reads)
            "ccs_reads": None,
            "s_reads": None,
            # lima (primers)
            # TODO: low scores
            # TODO: incorrect primer pairs
            # refine (reads)
            "fl_reads": None,
            "flnc_reads": None,
            "poly_reads": None,
            # correct (barcodes)
            "non-missing_reads": None,
            "yield_reads": None,
        }
        if header:
            f1.write("\t".join(data.keys()) + "\n")
            header = False
        
        f_skera = snakemake.input.skera[n]  # f'{config["skera_dir"]}/{sample}.summary.json'
        with open(f_skera) as f2:
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = dct["id"]
            value = str(dct["value"])
            if key == "reads":
                data["ccs_reads"] = value
            elif key == "s_reads":
                data["s_reads"] = value

        f_refine = snakemake.input.refine[n]  # f'{config["isoseq_refine_dir"]}/{sample}.filter_summary.report.json'
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
                data["poly_reads"] = value

        f_correct = snakemake.input.correct[n]  # f'{config["isoseq_correct_dir"]}/{sample}.report.json'
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
