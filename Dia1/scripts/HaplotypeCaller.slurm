#!/usr/bin/env bash
#SBATCH --job-name=HaplotypeCaller
#SBATCH --output=HaplotypeCaller.%j
#SBATCH --ntasks=36

module load tools/gatk/gatk-4.1.2.0

INPUTDIR=/home/marcopretti/ligandoma/LBB/results
OUTDIR=$INPUTDIR
sampleID='amostra-lbb'
BAM=$INPUTDIR/$sampleID'_hg38.recal.bam'

# databases and references
reference=/home/marcopretti/ligandoma/LBB/data/grch38.chr22.fasta

# Create dictionnary for reference
gatk CreateSequenceDictionary -R $reference

# HaploptypeCaller
echo -e "HaplotypeCaller\n"
gatk --java-options "-Xmx4g" HaplotypeCaller  \
    -R $reference \
    -I $BAM \
    -O $OUTDIR/$sampleID.HaplotypeCaller.vcf \
    -L chr22 

gatk VariantFiltration \
    -V $OUTDIR/$sampleID.HaplotypeCaller.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O $OUTDIR/$sampleID.HaplotypeCaller.filtered.vcf


echo -e "\nEnd of workflow\n"
date
module unload tools/gatk/gatk-4.1.2.0

exit
