#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem=100GB
#SBATCH --partition=cmb
#SBATCH --nodes=1
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jingwenr@usc.edu

path=/home/cmb-16/mjc/jingwenr/LRA_revision/SV_calling/canu_assembly
samtools view -@16 -F 260 ${path}/lra.paternal.bam > ${path}/lra.paternal.sam
python3 ParseAlignment.py  ${path}/lra.paternal.sam > ${path}/lra.paternal.parsed.sam

