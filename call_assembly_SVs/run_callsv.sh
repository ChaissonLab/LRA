#!/bin/bash
# change the "path", "ref", "centromere" in callassemblysv.json 
# the input sorted sam would be: lra.maternal.sam, lra.paternal.sam

snakemake -j 2 -ps callassemblysv.snakefile 
