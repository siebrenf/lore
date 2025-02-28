#!/usr/bin/env python

import subprocess as sp
from os.path import dirname, join
import sys
from snakemake import __version__ as version_sm
from pblrp import __version__ as version_pblrp


def cli():
    """
    wrap snakemake with the Snakefile preset
    """
    script_dir = dirname(__file__)
    pkg_dir = dirname(script_dir)
    snakefile = join(pkg_dir, "Snakefile")
    cmd = " ".join([f"snakemake --snakefile {snakefile} --use-conda"] + sys.argv[1:])

    if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
        print(f"Pac. Bio. Long Read Pipeline v{version_pblrp}")
        print(f"snakemake v{version_sm}")
        return

    retcode = sp.call(cmd, shell=True)
    print("")  # newline
    if retcode != 0:
        sys.exit(retcode)


if __name__ == "__main__":
    cli()
