import os

from snakemake.logging import logger


# reference genome
config["genome"] = os.path.basename(config["genome_fasta"].rsplit(".fa", 1)[0])
config["genome_dir"] = os.path.dirname(config["genome_fasta"])

# directories
if "results_dir" not in config:
    config["results_dir"] = os.path.join(os.getcwd(), "results")
if not os.path.exists(config["results_dir"]):
    os.makedirs(config["results_dir"])

for directory in [
    "benchmark_dir",
    "scripts_dir",
    "ccs_dir",  # 0-CCS/
    "skera_dir",  # 1-Sreads/
    "lima_dir",
    "isoseq_tag_dir",
    "isoseq_refine_dir",
    "isoseq_correct_dir",
    "isoseq_dedup_dir",  # 2-DeduplicatedReads/
    "isoseq_collapse_dir",
    "pbmm2_dir",
    "pigeon_sort_dir",
    "pigeon_classify_dir",  # 3-CollapsedTranscripts/
    # "pigeon_filter_dir",
    "pigeon_report_dir",
    "seurat_dir",  # 4-SeuratMatrix/
    "qc_dir",
]:
    if config.get(directory) is None:
        suffix = directory.rsplit("_", 1)[0]
        config[directory] = os.path.join(config["results_dir"], suffix)

for k, v in config.items():
    if isinstance(k, str) and k.endswith("_dir"):
        if not os.path.exists(v):
            os.makedirs(v)

# print config
logger.info("===================LORE Config===================")
for key in config:
    if config[key] not in ["", False, 0, "None"]:
        logger.info(f"{key: <23}: {config[key]}")
logger.info("=================================================")
logger.info("")
