import os

from snakemake.logging import logger


# directories
if "results_dir" not in config:
    config["results_dir"] = os.path.join(os.getcwd(), "results")
if not os.path.exists(config["results_dir"]):
    os.makedirs(config["results_dir"])

for directory in [
    "benchmark_dir",
    "ccs_dir",  # 0-CCS/
    "sreads_dir",  # 1-Sreads/
    "lima_dir",
    "isoseq_tag_dir",
    "isoseq_refine_dir",
    "isoseq_correct_dir",
    "dedup_dir",  # 2-DeduplicatedReads/
    "isoseq_collapse_dir",
    "pbmm2_dir",
    "pigeon_dir",  # 3-CollapsedTranscripts/
    "seurat_dir",  # 4-SeuratMatrix/
]:
    if config.get(directory) is None:
        suffix = directory.rsplit("_", 1)[0]
        config[directory] = os.path.join(config["results_dir"], suffix)

for k, v in config.items():
    if isinstance(k, str) and k.endswith("_dir"):
        if not os.path.exists(v):
            os.makedirs(v)

# reference genome
config["genome"] = os.path.basename(config["genome_fasta"].rsplit(".fa", 1)[0])

# print config
logger.info("===================LORE Config===================")
for key in config:
    if config[key] not in ["", False, 0, "None"]:
        logger.info(f"{key: <23}: {config[key]}")
logger.info("=================================================")
logger.info("")
