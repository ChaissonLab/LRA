LRA, the long read aligner.

Much of the functionality is still rolling out for lra, but your feedback is welcome.

Installation:

`git clone --recursive https://github.com/ChaissonLab/LRA.git -b master`

`make`


Running:

1. Index the genome. 
LRA uses a two-tiered word-indexing scheme (global and local).  
Both can be built at once using index:

`lra index hg38.fa`

`lra index hg38.fa -W 15 -K 17 -w 10 -k 5`

Alternatively they may be built separately:
`lra global hg38.fa`

`lra global hg38.fa -W 10`

`lra local hg38.fa`

`lra local hg38.fa -w 10`

Right now minimzers are used, but one of the purposes of lra is to test
out alternative indexing schemes.


2. Map.

The reads may be in fasta, sam, or bam format.

`lra align hg38.fa reads.fasta`

`lra align hg38.fa reads.bam`

`lra align hg38.fa reads.sam`

run `lra align -h` for a list of additional options.










