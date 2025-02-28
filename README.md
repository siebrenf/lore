# LOng-read Repeat Element pipeline for PacBio single-cell MAS-seq data

## MAS-Seq bioinformatics workflow
![](imgs/workflow.png)

## Repeat elements
WIP

## How to run LORE:

Test your config:
```[bash]
snakemake --snakefile lore/Snakefile --configfile config.yaml --dry-run
```

Run your config
```[bash]
snakemake --snakefile lore/Snakefile --configfile config.yaml --resources mem_mb=12_000 -j 14 > log.txt 2>&1
```

## Further reading:
  - [PacBio docs](https://isoseq.how/getting-started.html#recommended-single-cell-iso-seq-workflow)
  - [PacBio repos](https://github.com/PacificBiosciences/pbbioconda)
  - [PacBio datasets](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/)

## WIP:
  - RE mapping: add a second pass to the alignment stage to specifically map repeat elements. 
  - Optionally merge SMRT cells:
    https://isoseq.how/umi/cli-workflow.html#step-4b---merge-smrt-cells