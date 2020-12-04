## Users' Guide
LRA, the long read aligner for sequences and assembly contigs
LRA is a sequence alignment program that aligns long reads from single-molecule sequencing (SMS) instruments, or megabase-scale contigs from SMS assemblies. LRA implements seed chaining sparse dynamic programming with a convex gap function to read and assembly alignment, which is also extended to allow for inversion cases. Through the Truvari analysis of LRA, Minimap2 and NGM-LR alignments. LRA achieves higher f1 score over HG002 HiFi, CLR and ONT datasets. 

## Installation

Install by bioconda: `conda install -c bioconda lra`

Install from github: `git clone --recursive https://github.com/ChaissonLab/LRA.git -b master && cd LRA && make`

## General Usage
1. LRA needs to first build a two-tiered minimizer indexes (global and local) for the reference before mapping. Both can be built at once using commands:
```
lra index -CCS/CLR/ONT/CONTIG ref.fa
```
LRA has different parameters setting for the index when aligning reads from different sequencing instruments (CCS/CLR/ONT/CONTIG). You can also custimize the parameters. Details see `lra index --help`. LRA takes a few minutes to index the human reference genome

Alternatively the global and local indexes may be built separately: 
```
lra global -CCS/CLR/ONT/CONTIG ref.fa
lra local -CCS/CLR/ONT/CONTIG ref.fa
```

2. LRA takes reads fasta, fastq or bam format in the mapping step. The output format can be SAM, PAF, BED, and pairwise alignment. Details see `lra align --help`. The usage of multiple threads can be specified by `-t`. LRA uses the same base algorithm for mapping all datatypes with different parameters settings. It is recommended to choose among `CCS/CLR/ONT/CONTIG` based on the accuracy and average length of the input reads. 

```
lra align -CCS/CLR/ONT/CONTIG ref.fa read.fa -t 16 -p s > output.sam  
lra align -CCS/CLR/ONT/CONTIG ref.fa read.fa -t 16 -p p > output.paf  
lra align -CCS/CLR/ONT/CONTIG ref.fa read.fa -t 16 -p b > output.bed
```

If you have read.fa.gz, you may pip the read.fa to lra.

```
zcat read.fa.gz | lra align -CCS ref.fa /dev/stdin -t -p s > output.sam
```

## Output Format

LRA uses a set of customized tags in SAM and PAF output.

|Tag|Type  |Description                                  						|
|:--|:----:|:-------------------------------------------------------------------|
|NM |i |Number of matches in the alignment                                      |
|NX |i |Number of mismatches in the alignment                                   |
|ND |i |Number of bases of deletions in the alignment                           |
|TD |i |Number of deletions in the alignment                                    |
|NI |i |Number of bases of insertions in the alignment                          |
|TI |i |Number of insertions in the alignment                                   | 
|NV |f |The alignment score                                                     |
|TP |A |Type of aln, P/primary, S/secondary, I/inversion                        |
|RT |i |runtime   									                            |
|CG |z |CIGAR string        						                            |	 
|AO |i |This number shows the order of the aligned segment when a read is split |



















