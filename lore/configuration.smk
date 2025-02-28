import os

from snakemake.logging import logger


# directories
assert "results" in config
assert os.path.exists(config['results'])

for directory in ["ccs_dir", "lima_dir"]:
    if config.get(directory) is None:
        suffix = directory.split("_")[0]
        config[directory] = os.path.join(config["results_dir"], suffix)

for k, v in config.items():
    if isinstance(k, str) and k.endswith("_dir"):
        if not os.path.exists(v):
            os.makedirs(v)

# print config
logger.info("Config")
for key in config:
    if config[key] not in ["", False, 0, "None"]:
        logger.info(f"{key: <23}: {config[key]}")
logger.info("")
