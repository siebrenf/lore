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
}
for n in range(len(snakemake.input.raw)):
    for fname in [
        snakemake.input.raw[n],  # f'{config["pigeon_classify_dir"]}/{sample}.summary.txt'
        snakemake.input.raw[n].replace(".summary.txt", "_classification.filtered.summary.txt")
    ]:
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
