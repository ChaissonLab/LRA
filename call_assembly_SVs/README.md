To call SVs from lra alignment of assemblys, run `snakemake -ps callassemblysv.snakefile -j 4`. 
The input `bam` file must in the format of l`ra.${hap}.bam`. You can define the `hap` in callassemblysv.json as `[h1, h2]` or `[`maternal`, `paternal`]`.There would be five outputs under `data_path`: lra.INS.hom.vcf, lra.DEL.hom.vcf, lra.INS.het.vcf, lra.DEL.het.vcf, lra.vcf
You must modify the `ref`(reference path), `data_path`(where you store the alignment `bam`), `"samtools"` (path to `samtools`) in callassemblysv.json. 
