To create the rulegraph:
```[bash]
snakemake --snakefile lore/Snakefile --configfile config.yaml --forceall --rulegraph | dot -Tpng  -Gsize=2,4\! -Gdpi=500 > imgs/rulegraph.png
```
