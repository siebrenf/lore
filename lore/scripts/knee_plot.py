import logging
import os
import sys
import traceback
from contextlib import redirect_stderr, redirect_stdout

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# snakemake variables
logfile = snakemake.log[0]
tsv = snakemake.input[0]
estimate_percentile = snakemake.params.estimate_percentile
output = snakemake.output[0]

# capture exceptions & traceback
logging.basicConfig(
    filename=logfile,
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


sys.excepthook = handle_exception

# capture stdout & stderr
with open(logfile, "a") as log, redirect_stdout(log), redirect_stderr(log):
    df = pd.read_table(tsv)

    df = df[df["BarcodeType (Cell/Group or Molecular/UMI)"] == "Group/CellBarcode"]
    counts = df["NumberOfMolecules"].to_list()
    counts.sort(reverse=True)
    print(f"{len(counts)=}")

    df = df[df["RealCell"] == "cell"]
    ncells = df["FrequencyRank"].max()
    print(f"{ncells=}")

    plt.figure(figsize=(5, 3.5))
    c = plt.axes([0.125, 0.125, 0.8, 0.8])

    pink = (223 / 255, 25 / 255, 149 / 255)

    c.plot(range(len(counts)), counts, lw=1, color="grey", zorder=10)
    c.plot(range(ncells), counts[:ncells], lw=1.5, color=pink, zorder=11)

    if estimate_percentile is not None:
        lb = int(estimate_percentile)
        if lb < 80:
            print(
                "WARNING: lower bound for the percentile has a min of 80. Truncating to 80.",
                file=sys.stderr,
            )
            lb = 80
        colors = [
            (249 / 255, 157 / 255, 65 / 255),
            (255 / 255, 102 / 255, 204 / 255),
            (141 / 255, 109 / 255, 176 / 255),
            (68 / 255, 168 / 255, 223 / 255),
            (0 / 255, 184 / 255, 196 / 255),
            (106 / 255, 191 / 255, 105 / 255),
            (155 / 255, 167 / 255, 186 / 255),
            (225 / 255, 106 / 255, 44 / 255),
            (223 / 255, 25 / 255, 149 / 255),
            (95 / 255, 36 / 255, 159 / 255),
            (19 / 255, 131 / 255, 198 / 255),
            (0 / 255, 156 / 255, 162 / 255),
            (0 / 255, 157 / 255, 78 / 255),
            (103 / 255, 111 / 255, 127 / 255),
            (195 / 255, 74 / 255, 33 / 255),
            (158 / 255, 13 / 255, 59 / 255),
            (67 / 255, 31 / 255, 103 / 255),
            (0 / 255, 84 / 255, 150 / 255),
            (13 / 255, 77 / 255, 101 / 255),
            (0 / 255, 91 / 255, 66 / 255),
        ]
        percs, cutoffs = [], [0] * (99 - lb + 1)
        for percentile in range(99, lb - 1, -1):
            percs.append(np.percentile(counts, percentile) * 10)
        print(f"{percs=}")
        for count in counts:
            for j, p in enumerate(percs):
                if count > p:
                    cutoffs[j] += 1
        for i, cutoff in enumerate(cutoffs):
            c.plot(
                cutoff,
                counts[cutoff],
                marker="o",
                ms=4,
                c=colors[i],
                label=str(list(range(99, lb - 1, -1))[i]),
                zorder=15,
            )
        c.legend(loc=3, prop={"size": 6})

    c.set_xlabel(r"Cell # (log$_{10}$)")
    c.set_xscale("log")
    c.set_yscale("log")
    c.set_ylabel(r"log$_{10}$(# of UMIs)")
    sample_name = os.path.basename(tsv).rsplit(".", 2)[0]
    c.set_title(f"{sample_name}")
    plt.savefig(output, dpi=600)
