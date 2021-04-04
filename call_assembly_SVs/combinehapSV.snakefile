import os
import tempfile
import subprocess
import os.path
import sys

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
# Config
configfile: "combinehapSV.json"
aligner = "mm2"
SVTYPE = ["INS", "DEL"]
hap = config["hap"]
path = config["data_path"]
samtools = config["samtools"]

rule all:
  input:
    SVTYPE_hom = expand(path + "/vcf/" + aligner + ".{SVTYPE}.hom.vcf",  SVTYPE = SVTYPE),
    SVTYPE_het_hap = expand(path + "/vcf/" + aligner + ".{hap}.{SVTYPE}.het.vcf", hap = hap, SVTYPE = SVTYPE),
    vcf = path + "/vcf/" + aligner + ".vcf"

rule vcfhead:
  input:
    originalvcf = path + "/vcf/" + aligner + ".{hap}.original.vcf"
  output:
    head = path + "/vcf/" + aligner + ".{hap}.head.vcf"
  shell:"""
    awk '$1~/^#/' {input.originalvcf} > {output.head}
  """

rule FilterBySizeAndSVTYPE:
  input:
    originalvcf = path + "/vcf/" + aligner + ".{hap}.original.vcf"
  output:
    INSfldvcf_original = path + "/vcf/" + aligner + ".{hap}.fld.INS.original.vcf",
    DELfldvcf_original = path + "/vcf/" + aligner + ".{hap}.fld.DEL.original.vcf",
  shell:"""
    awk '$1!~/^#/' {input.originalvcf} | awk '{{ if (length($4) > length($5) && length($4) - length($5) >= 40) {{print$0;}} }}' |\
    awk 'BEGIN{{OFS="\\t";}} {{print$1,$2,$2+length($4)-length($5),$3, $4,$5,$6,$7,$8,$9,$10;}}' > {output.DELfldvcf_original}
    awk '$1!~/^#/' {input.originalvcf} | awk '{{ if (length($5) > length($4) && length($5) - length($4) >= 40) {{print$0;}} }}' |\
    awk 'BEGIN{{OFS="\\t";}} {{print$1,$2,$2+length($5)-length($4),$3, $4,$5,$6,$7,$8,$9,$10;}}' > {output.INSfldvcf_original}
  """

rule FilterByCentromereAndAltChromAndAddID:
  input:
    fldvcf_original = path + "/vcf/" + aligner + ".{hap}.fld.{SVTYPE}.original.vcf",
  output:
    fldvcf = temp(path + "/vcf/" + aligner + ".{hap}.fld.{SVTYPE}.vcf")
  shell:"""
    bedtools intersect -v -a {input.fldvcf_original} -b hg19_centromere.wide.bed |\
    awk '{{if ($1!~/^GL/ || $1!~/^hs37d5/ || $1!~/^NC/) {{print$0;}}}}' |\
    awk -v aln=aligner -v type={wildcards.SVTYPE} -v h={wildcards.hap} \
    'BEGIN{{OFS="\\t";}}{{print$1,$2,$3, aln"."type"."h"."NR,$5,$6,$7,$8,$9,$10,$11;}}' > {output.fldvcf}
  """  

rule mergeCalls:
  input:
    fldvcf = path + "/vcf/" + aligner + ".{hap}.fld.{SVTYPE}.vcf"
  output:
    fldvcf_merge = temp(path + "/vcf/" + aligner + ".{hap}.fld.{SVTYPE}.merge.vcf")
  shell:"""
    bedtools intersect -wao -a {input.fldvcf} -b {input.fldvcf} |\
    awk 'BEGIN{{OFS="\\t"}}{{if ($23/($3-$2) >= 0.1 && $23/($14-$13) >= 0.1) {{print$0;}}}}' |\
    bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11 -c 15 -o distinct |\
    python3 mergeSV.py > {output.fldvcf_merge}
  """


rule Homcalls:
  input:
    fldvcf_m = path + "/vcf/" + aligner + ".maternal.fld.{SVTYPE}.merge.vcf",
    fldvcf_p = path + "/vcf/" + aligner + ".paternal.fld.{SVTYPE}.merge.vcf"
  output:
    SVTYPE_hom = path + "/vcf/" + aligner + ".{SVTYPE}.hom.vcf"
  shell:"""
    bedtools intersect -u -a {input.fldvcf_m} -b {input.fldvcf_p} -f 0.1 -r  |\
    awk 'BEGIN{{OFS="\\t";}} {{print$1,$2,$4,$5,$6,$7,"PASS",$9,$10,$11;}}' > {output.SVTYPE_hom}
  """

rule Hetcalls:
  input:
    fldvcf_m = path + "/vcf/" + aligner + ".maternal.fld.{SVTYPE}.merge.vcf",
    fldvcf_p = path + "/vcf/" + aligner + ".paternal.fld.{SVTYPE}.merge.vcf"
  output:
    het_m = path + "/vcf/" + aligner + ".maternal.{SVTYPE}.het.vcf",
    het_p = path + "/vcf/" + aligner + ".paternal.{SVTYPE}.het.vcf"
  shell:"""
    bedtools intersect -v -a {input.fldvcf_m} -b {input.fldvcf_p} -f 0.1 -r |\
    awk 'BEGIN{{OFS="\t";}} {{print$1,$2,$4,$5,$6,$7,"PASS",$9, $10,$11;}}' > {output.het_m}
    bedtools intersect -v -a {input.fldvcf_p} -b {input.fldvcf_m} -f 0.1 -r |\
    awk 'BEGIN{{OFS="\t";}} {{print$1,$2,$4,$5,$6,$7,"PASS",$9,$10,$11;}}' > {output.het_p}
  """

rule combineHet:
  input:
    het_hap = expand(path + "/vcf/" + aligner + ".{hap}.{SVTYPE}.het.vcf", aligner = aligner, hap = hap, SVTYPE = SVTYPE),
    SVTYPE_hom = expand(path + "/vcf/" + aligner + ".{SVTYPE}.hom.vcf", SVTYPE = SVTYPE),
    head = path + "/vcf/" + aligner + ".maternal.head.vcf"
  output:
    srtvcf = path + "/vcf/" + aligner + ".srt.vcf",
    vcf = path + "/vcf/" + aligner + ".vcf",
    vcfgz = path + "/vcf/" + aligner + ".vcf.gz"
  params: 
    fai = config["fai"]
  shell:"""
    # rm {output.srtvcf}
    cat {input.head} {input.het_hap} {input.SVTYPE_hom} >> {output.srtvcf}
    bedtools sort -header -g {params.fai} -i {output.srtvcf}  > {output.vcf}
    cat {output.vcf}| bgzip -c >| {output.vcfgz}
    tabix -f {output.vcfgz}
  """













  


