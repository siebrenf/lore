"""
source: https://github.com/bvaldebenitom/SoloTE/blob/main/SoloTE_RepeatMasker_to_BED.py
commit: https://github.com/bvaldebenitom/SoloTE/commit/b454b27b051dba076c1d02f93ba0018982dce1ef
"""
import pandas
from colorama import Fore, init
from tqdm import tqdm


rmsk_filename = snakemake.output.fasta
rmsk_bed_filename = snakemake.output.bed

rmsk_out = pandas.read_csv(
    rmsk_filename,
    compression="gzip",
    skiprows=3,
    header=None,
    sep=" ",
    skipinitialspace=True,
)
rmsk_out.columns = [
    "SW_score",
    "percDiv",
    "percDel",
    "percIns",
    "querySeq",
    "queryStart",
    "queryEnd",
    "queryLeft",
    "strand",
    "matchingRepeat",
    "repeatClass_Family",
    "repeatBegin",
    "repeatStart",
    "repeatEnd",
    "ID",
]
repeatinfo = rmsk_out["repeatClass_Family"].str.split("/", expand=True)
repeatinfo = repeatinfo.fillna(value="-")
rmsk_out["repeatID"] = (
    rmsk_out["matchingRepeat"] + ":" + repeatinfo[1] + ":" + repeatinfo[0]
)
rmsk_out["percDiv"] = rmsk_out["percDiv"].astype(str)
rmsk_out["queryStart"] = rmsk_out["queryStart"].astype(str)
rmsk_out["queryEnd"] = rmsk_out["queryEnd"].astype(str)
rmsk_out["fixedStrand"] = rmsk_out["strand"].str.replace("C", "-")
rmsk_out["te_name"] = (
    rmsk_out["querySeq"]
    + "|"
    + rmsk_out["queryStart"]
    + "|"
    + rmsk_out["queryEnd"]
    + "|"
    + rmsk_out["repeatID"]
    + "|"
    + rmsk_out["percDiv"]
    + "|"
    + rmsk_out["fixedStrand"]
)
rmsk_out_bed = rmsk_out[
    ["querySeq", "queryStart", "queryEnd", "te_name", "percDiv", "fixedStrand"]
]
rmsk_out_bed = rmsk_out_bed[
    rmsk_out["repeatClass_Family"].str.contains("LINE|SINE|LTR|DNA|RC")
]
rmsk_out_bed = rmsk_out_bed[
    ~rmsk_out_bed["querySeq"].str.contains("chrna|_fix|_random|_alt|chrUn")
]
rmsk_out_bed.to_csv(rmsk_bed_filename, sep="\t", header=None, index=False)
