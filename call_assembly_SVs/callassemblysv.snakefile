import os
import tempfile
import subprocess
import os.path
import sys

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
# Config
configfile: "callassemblysv.json"
aligner = "lra"
hap = ["maternal", "paternal"]
SVTYPE = ["INS", "DEL"]
path = config["path"]
centromere = config["centromere"]


rule all:
  input:
    SVTYPE_hom = expand(path + aligner + ".{SVTYPE}.hom.vcf",  SVTYPE = SVTYPE),
    SVTYPE_het_hap = expand(path + aligner + ".{hap}.{SVTYPE}.het.vcf", hap = hap, SVTYPE = SVTYPE),
    head = expand(path + aligner + ".{hap}.head.vcf", hap = hap),
    vcf = path + aligner + ".vcf"


rule removehead: 
  input:
    sam = path + aligner + ".{hap}.sam" 
  output: 
    nohead_sam = temp(path + aligner + ".{hap}.noheader.sam"),
    head = temp(path + aligner + ".{hap}.head.sam")
  shell:"""
    awk '$1!~/\@/' {input.sam} > {output.nohead_sam}
    awk '$1~/\@/' {input.sam} > {output.head}
  """

rule parseAlignment:
  input:
    nohead_sam = path + aligner + ".{hap}.noheader.sam"
  output:
    parsedsam_nohead = temp(path + aligner + ".{hap}.parsed.nohead.sam")
  shell:""" 
    python3 ParseAlignment.py --inputsam {input.nohead_sam} --parsedsam {output.parsedsam_nohead}
  """

rule addhead:
  input:
    parsedsam_nohead = path + aligner + ".{hap}.parsed.nohead.sam",
    head = path + aligner + ".{hap}.head.sam"
  output:
    parsedsam = temp(path + aligner + ".{hap}.parsed.sam")
  shell:"""
      cat {input.head} {input.parsedsam_nohead} >> {output.parsedsam}
  """

rule CallSV:
  input:
    parsedsam = path + aligner + ".{hap}.parsed.sam"
  output:
    originalvcf = temp(path + aligner + ".{hap}.original.vcf"),
    head = temp(path + aligner + ".{hap}.head.vcf")
  params:
    ref = config["ref"]
  shell:"""
    python3 SamToVCF.py --ref {params.ref} --sam {input.parsedsam}  > {output.originalvcf}
    awk '$1~/\#/' {output.originalvcf} > {output.head}
  """

rule FilterBySizeAndSVTYPE:
  input:
    originalvcf = path + aligner + ".{hap}.original.vcf"
  output:
    INSfldvcf_original = temp(path + aligner + ".{hap}.fld.INS.original.vcf"),
    DELfldvcf_original = temp(path + aligner + ".{hap}.fld.DEL.original.vcf")
  shell:"""
    awk '$1!~/\#/' {input.originalvcf} | awk '{{ if (length($4) > length($5) && length($4) - length($5) >= 40) {{print$0;}} }}' |\
    awk 'BEGIN{{OFS="\\t";}} {{print$1,$2,$2+length($4)-length($5),$3, $4,$5,$6,$7,$8,$9,$10;}}' > {output.DELfldvcf_original}
    awk '$1!~/\#/' {input.originalvcf} | awk '{{ if (length($5) > length($4) && length($5) - length($4) >= 40) {{print$0;}} }}' |\
    awk 'BEGIN{{OFS="\\t";}} {{print$1,$2,$2+length($5)-length($4),$3, $4,$5,$6,$7,$8,$9,$10;}}' > {output.INSfldvcf_original}
  """

rule FilterByCentromereAndAltChromAndAddID:
  input:
    fldvcf_original = path + aligner + ".{hap}.fld.{SVTYPE}.original.vcf",
  output:
    fldvcf = temp(path + aligner + ".{hap}.fld.{SVTYPE}.vcf")
  params:
    centromere = config["centromere"]
  shell:"""
    bedtools intersect -v -a {input.fldvcf_original} -b {params.centromere} |\
    awk '{{if ($1!~/^GL/ || $1!~/^hs37d5/ || $1!~/^NC/) {{print$0;}}}}' |\
    awk -v aln=aligner -v type={wildcards.SVTYPE} -v h={wildcards.hap} 'BEGIN{{OFS="\\t";}}{{print$1,$2,$3, aln"."type"."h"."NR,$5,$6,$7,$8,$9,$10,$11;}}' > {output.fldvcf}
  """  

rule mergeCalls:
  input:
    fldvcf = path + aligner + ".{hap}.fld.{SVTYPE}.vcf"
  output:
    fldvcf_merge = temp(path + aligner + ".{hap}.fld.{SVTYPE}.merge.vcf")
  shell:"""
    bedtools intersect -wao -a {input.fldvcf} -b {input.fldvcf} |\
    awk 'BEGIN{{OFS="\\t"}}{{if ($23/($3-$2) >= 0.1 && $23/($14-$13) > 0.1) {{print$0;}}}}' |\
    bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11 -c 15 -o distinct |\
    python3 mergeSV.py > {output.fldvcf_merge}
  """


rule Homcalls:
  input:
    fldvcf_m = path + aligner + ".maternal.fld.{SVTYPE}.merge.vcf",
    fldvcf_p = path + aligner + ".paternal.fld.{SVTYPE}.merge.vcf"
  output:
    SVTYPE_hom = path + aligner + ".{SVTYPE}.hom.vcf"
  shell:"""
    bedtools intersect -u -a {input.fldvcf_m} -b {input.fldvcf_p} -f 0.3 -r  |\
    awk 'BEGIN{{OFS="\\t";}} {{print$1,$2,$4,$5,$6,$7,$9,$9,$10,$11;}}' > {output.SVTYPE_hom}
  """

rule Hetcalls:
  input:
    fldvcf_m = path + aligner + ".maternal.fld.{SVTYPE}.merge.vcf",
    fldvcf_p = path + aligner + ".paternal.fld.{SVTYPE}.merge.vcf"
  output:
    het_m = path + aligner + ".maternal.{SVTYPE}.het.vcf",
    het_p = path + aligner + ".paternal.{SVTYPE}.het.vcf"
  shell:"""
    bedtools intersect -v -a {input.fldvcf_m} -b {input.fldvcf_p} -f 0.3 -r |\
    awk 'BEGIN{{OFS="\t";}} {{print$1,$2,$4,$5,$6,$7,$9,$9,$10,$11;}}' > {output.het_m}
    bedtools intersect -v -a {input.fldvcf_p} -b {input.fldvcf_m} -f 0.3 -r |\
    awk 'BEGIN{{OFS="\t";}} {{print$1,$2,$4,$5,$6,$7,$9,$9,$10,$11;}}' > {output.het_p}
  """

rule combineHet:
  input:
    het_hap = expand(path + aligner + ".{hap}.{SVTYPE}.het.vcf", aligner = aligner, hap = hap, SVTYPE = SVTYPE),
    SVTYPE_hom = expand(path + aligner + ".{SVTYPE}.hom.vcf", SVTYPE = SVTYPE),
    head = path + aligner + ".maternal.head.sam"
  output:
    vcf = path + aligner + ".vcf"
  shell:"""
    cat {input.head} {input.het_hap} {input.SVTYPE_hom} >> {output.vcf}
  """













  


