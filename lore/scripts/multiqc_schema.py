import logging
import os
import re
import sys
import traceback
from contextlib import redirect_stdout, redirect_stderr


# snakemake variables
logfile = snakemake.log[0]
f_in = snakemake.params.schema
f_out = snakemake.output.schema[0]
config = snakemake.params.config

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
    # shutil.copy(f_in, f_out)

    with open(f_in) as cookie_schema:
        cookie = cookie_schema.read()

    while len(list(re.finditer("{{.*}}", cookie))):
        match_obj = next(re.finditer("{{.*}}", cookie))
        span = match_obj.span()
        match = eval(match_obj.group(0)[2:-2])  # without the curly braces
        cookie = cookie[: span[0]] + match + cookie[span[1] :]

    with open(f_out, "w") as out_file:
        out_file.write(cookie)
