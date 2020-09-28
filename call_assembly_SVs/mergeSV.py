#!/usr/bin/env python
import sys
import re

SV = []
for line in sys.stdin:
    Chrom,Start,End,ID,Ref,Alt,Qual,Filter,Info,Format,Sample,intersect = line.rstrip("\n").split("\t")
    SV.append([Chrom, Start, End, ID, Ref, Alt, Qual, Filter, Info, Format, Sample, intersect])

seen = set()
remove = [0] * len(SV)
for id, i in enumerate(SV):
    if i[3] not in seen:
        for j in re.split(",",i[11]):
            seen.add(j)
    else:
        remove[id] = 1
new = []
for id, i in enumerate(remove):
    if i == 0:
        new.append(SV[id][0:11])

for x in new:
  sys.stdout.write("\t".join(x)+"\n")


