#!/usr/bin/env python
import sys
import pysam
import argparse


ap=argparse.ArgumentParser(description="Parse a sam into a vcf")
ap.add_argument("--sam", help="Input sam file", default=None)
ap.add_argument("--bam", help="Input bam file", default=None)
ap.add_argument("--ref", help="Input reference file", required=True)
ap.add_argument("--sample", help="Sample", default="unknown")
ap.add_argument("--minLength", help="Minimum length of variant", type=int, default=1)
ap.add_argument("--chrom", help="Process this region.", default=None)
args=ap.parse_args()

sam=pysam.AlignmentFile(args.sam)

ref=pysam.FastaFile(args.ref)

fai=open(args.ref+".fai")
contigHeader=""
for line in fai:
    vals=line.split()
    contigHeader+="##contig=<ID="+vals[0]+",length="+vals[1]+">\n"

sys.stdout.write("""##fileformat=VCFv4.2
##fileDate=20190102
##source=SamToVCF
""")
sys.stdout.write("##reference="+args.ref+"\n")
sys.stdout.write(contigHeader)
sys.stdout.write("""##INFO=<ID=QNAME,Number=1,Type=String,Description="Name of query sequence">
##INFO=<ID=QSTART,Number=1,Type=Integer,Description="Position of query sequence">
##INFO=<ID=QSTRAND,Number=1,Type=String,Description="Contig strand">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}
""".format(args.sample))

for aln in sam.fetch(contig=args.chrom):
    refStart = aln.reference_start
    refEnd   = aln.reference_end
    refSeq   = ref.fetch(aln.reference_name, refStart, refEnd).upper()
    b=aln.get_aligned_pairs()
    i=0
    strand="+"
    if aln.is_reverse:
        strand = "-"
    lastQuery=0
    lastRef=0
    while i < len(b) and b[i][0] == None or b[i][1] == None:
        i+=1
    if i >= len(b):
        continue

    while i < len(b):

        varPos=None
        if b[i][0] != None and b[i][1] != None:
            while i < len(b) and b[i][0] != None and b[i][1] != None:
                lastQuery=b[i][0]
                lastRef=b[i][1]

                if aln.seq[b[i][0]] != refSeq[b[i][1]-refStart]:
                    refVarSeq=refSeq[b[i][1]-refStart]                    
                    readVarSeq=aln.seq[b[i][0]]
                    varPos=b[i][1]+1
                    queryPos=b[i][0]+1
                    i+=1
                    break
                else:
                    i+=1
            continue
        elif b[i][0] != None and b[i][1] == None:
            #
            # Starting an insertion.
            #
            startIndex=i

            while i < len(b) and b[i][0] != None and b[i][1] == None:
                i+=1
            if i >= len(b):
                break
            queryEnd=b[i][0]
            readVarSeq=aln.seq[lastQuery:queryEnd]
            refVarSeq=refSeq[lastRef-refStart]
            varPos=lastRef
            queryPos=lastQuery+1
        elif b[i][0] == None and b[i][1] != None:
            # starting a deletion
            startIndex=i
            delStart=lastRef
            while i < len(b) and b[i][0] == None and b[i][1] != None:
                lastRef=b[i][1]
                i+=1
            if i >= len(b):
                break

            readVarSeq=aln.seq[lastQuery]
            refVarSeq=refSeq[delStart-refStart:(lastRef+1)-refStart]
            varPos=delStart
            queryPos=lastQuery+1


        if varPos is not None and abs(len(readVarSeq)-len(refVarSeq)) >= args.minLength:
            svLen = len(readVarSeq) - len(refVarSeq)
            if svLen < 0:
                svType="DEL"
            else:
                svType="INS"
            var=[aln.reference_name, varPos, ".", refVarSeq, readVarSeq, "60", "PASS", "SVTYPE={};SVLEN={};QNAME={};QSTART={};QSTRAND={}".format(svType, svLen, aln.query_name, queryPos, strand), "GT", "1/1"]
            vars=[str(s) for s in var]
            sys.stdout.write("\t".join(vars)+"\n")
        
    
    


