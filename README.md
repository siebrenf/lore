# Long-read Repeat Element pipeline for PacBio single-cell MAS-seq data


## LoRE workflow

<details>
<summary>Expand minimal workflow</summary>

![broken image](imgs/rulegraph.png)
</details>

<details>
<summary>Expand maximal workflow</summary>

![broken image](imgs/rulegraph_full.png)
</details>


## Underlying MAS-Seq workflow overview

<details>
<summary>Expand</summary>

![broken image](imgs/workflow.png)
</details>


## How install LORE:

Clone the repository:
```[bash]
git clone https://github.com/siebrenf/lore.git
```

Create the conda environment:
```[bash]
conda env create -n lore -f lore/requirements.yaml
conda activate lore
pip install -e .
```


## How to run LORE:

Change directory into the LORE folder.

Activate the conda environment:
```[bash]
conda activate lore
```

Update the `config.yaml`. 
- Adapters, primers and barcodes can be downloaded by lore, or can be placed inside the results directory. 
  The default results directory is `./results`.
- The genome and gene annotation need to be obtained manually.
  You will need to specify their locations in the config, as well as the symbol for the mitochondria.
- URLs to additional documentation for most rules (steps) in the workflow can be found in the code.
- Optional outputs (currently) include bigwigs (for track visualization) and a QC report.
  Both are recommended, but adds (some) computational load.

Test your config:
```[bash]
snakemake --snakefile lore/Snakefile --configfile config.yaml --dry-run
```

Run your config:
```[bash]
nice snakemake --use-conda --snakefile lore/Snakefile --configfile config.yaml --resources parallel_downloads=1 mem_mb=100_000 -j 60 > log.txt 2>&1
```


## Further reading:
  - [PacBio docs](https://isoseq.how/getting-started.html#recommended-single-cell-iso-seq-workflow)
  - [PacBio repos](https://github.com/PacificBiosciences/pbbioconda)
  - [PacBio datasets](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/)
  - [PacBio glossary](https://www.pacb.com/wp-content/uploads/2015/09/Pacific-Biosciences-Glossary-of-Terms.pdf)


## WIP:
  - TE pipelines:
    - https://doi.org/10.1016/j.isci.2023.108214
      - https://github.com/javiercguard/teNanoporePipeline  # usable
      - RepeatMaster
      - Dfam
    - https://doi.org/10.1093/nar/gkac794
      - https://github.com/bergmanlab/TELR  # looks good!
    - https://doi.org/10.1186/s13059-023-02911-2
      - https://github.com/DrosophilaGenomeEvolution/TrEMOLO  # looks good + snakemake!
    - https://doi.org/10.1186/s13100-017-0088-x
      - LoRTE  # python 2.7, dead link
    - de novo Repeat library construction:
      - could be useful for non-model organisms/strains
      - https://doi.org/10.1186/s12864-021-08117-9
      - https://github.com/kacst-bioinfo-lab/TE_ideintification_pipeline  # nice workflow figure
  
  - what to do with pbmm2 unmapped reads:
    - second pass to the alignment stage to specifically map repeat elements?
    - omit?
  - integrate genomepy to:
    - get a genome
    - create the minimap2 index
    - get an annotation
