import json
import os


data = {
    "Classifications": {
        "Sample": [],
        "Input": [],
        "Passed": [],
        "Removed": [],
        "Unique genes": [],
        "Unique transcripts": [],
    },
    "Classifications, by isoform": {
        "Sample": [],
        "Full splice match": [],
        "Incomplete splice match": [],
        "Novel in catalog": [],
        "Novel not in catalog": [],
        "Antisense": [],
        "Genic intron": [],
        "Genic genomic": [],
        "Intergenic": [],
        "Other": [],
    },
    "Classifications, by read": {
        "Sample": [],
        "Full splice match": [],
        "Incomplete splice match": [],
        "Novel in catalog": [],
        "Novel not in catalog": [],
        "Antisense": [],
        "Genic intron": [],
        "Genic genomic": [],
        "Intergenic": [],
        "Other": [],
    },
    "Junctions": {
        "Sample": [],
        "Known canonical": [],
        "Known non-canonical": [],
        "Novel canonical": [],
        "Novel non-canonical": [],
    },
    "Genes": {
        "Sample": [],
        "Known": [],
        "Novel": [],
    },
    "RT Switching": {
        "Sample": [],
        "All transcripts": [],
        "Unique transcripts": [],
        "All junctions": [],
        "Unique junctions": [],
    },
    "Filter reasons": {
        "Sample": [],
        "Intrapriming": [],
        "Monoexonic": [],
        "RT switching": [],
        "Low coverage/non-canonical": [],
    },
    "Classifications, by cell": {
        "Sample": [],
        "median_genes_per_cell": [],
        "median_transcripts_per_cell": [],
        "median_genes_per_cell_known": [],
        "median_transcripts_per_cell_known": [],
    },
    "Classifications, by mapping": {
        "Sample": [],
        "flnc_mapped_genome": [],
        "flnc_mapped_transcriptome": [],
        "flnc_mapped_transcriptome_excluding_ribomito": [],
    },
    "Classifications, by transcript": {
        "Sample": [],
        "transcripts_fsm": [],
        "transcripts_ism": [],
        "transcripts_nic": [],
        "transcripts_nnc": [],
    },
}
for n in range(len(snakemake.input.raw_summary)):
    for fname in [
        snakemake.input.raw_summary[n],  # f'{config["pigeon_classify_dir"]}/{sample}.summary.txt'
        snakemake.input.raw_summary[n].replace(".summary.txt", "_classification.filtered.summary.txt")
    ]:
        # Collect data from the summary.txt files.
        # The raw and filtered files are slightly different,
        #  which we correct for first
        if fname.endswith("_classification.filtered.summary.txt"):
            sample = os.path.basename(fname).rsplit("_classification.filtered.summary.txt", 1)[0]
            sample += " (filtered)"
            for key in data["RT Switching"]:
                if key == "Sample":
                    data["RT Switching"][key].append(sample)
                else:
                    data["RT Switching"][key].append(str(None))
            data["Classifications"]["Passed"].append(str(None))
            data["Classifications"]["Removed"].append(str(None))
        else:
            sample = os.path.basename(fname).rsplit(".", 2)[0]
            for key in data["Filter reasons"]:
                if key == "Sample":
                    data["Filter reasons"][key].append(sample)
                else:
                    data["Filter reasons"][key].append(str(None))

        subsection = "Classifications"
        data[subsection]["Sample"].append(sample)
        with open(fname) as f2:  
            for line in f2:
                line = line.strip()
                # skip empty lines
                if len(line) == 0 or line.startswith("---"):
                    pass

                # check subsection
                elif line == "Classifications":
                    pass  # edge case
                elif line in [
                    "Classifications",
                    "Classifications, by isoform",
                    "Classifications, by read",
                    "Junctions",
                    "Genes",
                    "Filter reasons",
                    "RT Switching",
                ]:
                    subsection = line
                    data[subsection]["Sample"].append(sample)

                # collect data
                else:
                    key, value = line.split(":")
                    key = key.strip()
                    if key == "Input transcripts":
                        key = "Input"
                    value = value.strip().split(" ")[0]  # strip suffix: " (x.xx%)"
                    data[subsection][key].append(value)

        # Collect data from the report.json files
        # The raw and filtered files are the same (yay)
        fname = fname.replace(".summary.txt", ".report.json")
        with open(fname) as f2:
            contents = json.load(f2)
        for dct in contents["attributes"]:
            key = dct["id"]
            if key in [
                "median_genes_per_cell",
                "median_transcripts_per_cell",
                "median_genes_per_cell_known",
                "median_transcripts_per_cell_known",
            ]:
                subsection = "Classifications, by cell"
                value = f'{dct["value"]:.2f}'
                data[subsection][key].append(value)
                if key == "median_genes_per_cell":
                    data[subsection]["Sample"].append(sample)
            elif key in [
                "flnc_mapped_genome",
                "flnc_mapped_transcriptome",
                "flnc_mapped_transcriptome_excluding_ribomito",
            ]:
                subsection = "Classifications, by mapping"
                value = f'{dct["value"]:.2f}'
                data[subsection][key].append(value)
                if key == "flnc_mapped_genome":
                    data[subsection]["Sample"].append(sample)
            elif key in [
                "transcripts_fsm",
                "transcripts_ism",
                "transcripts_nic",
                "transcripts_nnc",
            ]:
                subsection = "Classifications, by transcript"
                value = f'{dct["value"]:.2f}'
                data[subsection][key].append(value)
                if key == "transcripts_fsm":
                    data[subsection]["Sample"].append(sample)


# print(json.dumps(data, indent=4))
for subsection in data:
    fname = "pigeon_classify_" + subsection.lower().replace(",", "").replace(" ", "_") + ".tsv"
    fname = os.path.join(
        os.path.dirname(snakemake.output[0]),
        fname
    )
    with open(fname, "w") as f1:
        header = []
        for key in data[subsection]:
            header.append(key)
        f1.write("\t".join(header) + "\n")

        for n in range(len(data[subsection]["Sample"])):
            row = []
            for key in data[subsection]:
                # print(subsection, key, n)
                # print(data[subsection][key])
                row.append(data[subsection][key][n])
            f1.write("\t".join(row) + "\n")
