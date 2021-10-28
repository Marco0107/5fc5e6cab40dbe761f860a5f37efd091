#!/usr/bin/env bash

## Pre-processing workflow based on GATK4 IARC pipeline
## https://github.com/IARCbioinfo/gatk4-DataPreProcessing-nf

module load tools/gatk/gatk-4.1.2.0
module load tools/picard/2.18.5
module load tools/samtools/samtools-1.8
module load tools/bwa/bwa-0.7.17
qualimap=/home/marcopretti/bin/qualimap_v2.2.1/qualimap

INPUTDIR=/home/marcopretti/LBB/data
OUTDIR=/home/marcopretti/LBB/results
BINDIR=/home/marcopretti/bin
sampleID='amostra-lbb'
hgversion=hg38

input1=$INPUTDIR/$sampleID'_R1.fq.gz'
input2=$INPUTDIR/$sampleID'_R2.fq.gz'
TMPDIR=$OUTDIR/tmpdir.${SLURM_JOB_ID}

## databases and references
database=/Broad_bundle/hg38bundle
reference=$INPUTDIR/grch38.chr22.fasta
ref_dbsnp=$database/dbsnp_146.hg38.vcf.gz
ref_1kg=$database/1000G_phase1.snps.high_confidence.hg38.vcf.gz
ref_mills=$database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

echo -e "Using $SLURM_NTASKS threads"
mkdir $OUTDIR
mkdir $TMPDIR

## FastqToSam
echo -e "Fastq to Sam\n"
picard FastqToSam \
    F1=$input1 \
    F2=$input2 \
    O=$TMPDIR/$sampleID.ubam \
    SM=$sampleID \
    RG=rg0013 TMP_DIR=$TMPDIR

echo -e "Sorting with samtools\n"
samtools sort -o $TMPDIR/$sampleID.sorted.ubam -@ $SLURM_NTASKS $TMPDIR/$sampleID.ubam

## PICARD_FixMate
echo -e "FixMate\n"
picard FixMateInformation \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    I=$TMPDIR/$sampleID.sorted.ubam \
    O=$TMPDIR/$sampleID.fix.bam \
    TMP_DIR=$TMPDIR


## PICARD_SamToFastq
echo -e "Sam to Fastq\n"
picard \
    SamToFastq \
    VALIDATION_STRINGENCY=LENIENT \
    I=$TMPDIR/$sampleID.fix.bam \
    F=$TMPDIR/${sampleID}_1.fastq \
    F2=$TMPDIR/${sampleID}_2.fastq TMP_DIR=$TMPDIR


## BWA_mapping
echo -e "\nBWA Mapping\n"
cd $TMPDIR

echo -e "Using bwa with hg38 parameters\n"
bwa mem -t 24 -R "@RG\\tID:${sampleID}\\tSM:${sampleID}\\tPL:Illumina" $reference $TMPDIR/${sampleID}_1.fastq $TMPDIR/${sampleID}_2.fastq > $TMPDIR/$sampleID.sam
samtools view -b $TMPDIR/$sampleID.sam | samtools sort > $TMPDIR/$sampleID.bam


## Mark Duplicates
echo -e "Mark Duplicates\n"
picard MarkDuplicates \
            INPUT=$TMPDIR/${sampleID}.bam \
            OUTPUT=$TMPDIR/${sampleID}_dedup.bam \
            METRICS_FILE=$OUTDIR/${sampleID}_${hgversion}_duplicate_metrics.txt \
            TMP_DIR=$TMPDIR

## GATK_BaseRecalibrator_ApplyBQSR
echo -e "Base Recalibrator\n"    
samtools index $TMPDIR/${sampleID}_dedup.bam

gatk BaseRecalibrator \
    -R ${reference}.fasta \
    -I $TMPDIR/${sampleID}_dedup.bam \
    --use-original-qualities \
    -O $TMPDIR/${sampleID}.recal_data.table \
    -known-sites ${ref_dbsnp} \
    -known-sites ${ref_1kg} \
    -known-sites ${ref_mills}

gatk ApplyBQSR \
    -R ${reference}.fasta \
    -I $TMPDIR/${sampleID}_dedup.bam \
    --use-original-qualities \
    -O $OUTDIR/${sampleID}_${hgversion}.recal.bam \
    --bqsr-recal-file $TMPDIR/${sampleID}.recal_data.table


## QUALIMAP_BamQC
echo -e "Bam Quality Check\n"
samtools flagstat $OUTDIR/${sampleID}_${hgversion}.recal.bam > $OUTDIR/${sampleID}_${hgversion}.stats.txt
grep -v "^@" ${interList} | cut -f-3,5 | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,0,"."}' > $TMPDIR/tmp.bed
unset DISPLAY # requirement to avoid java DISPLAY  error according to qualimap

$qualimap bamqc \
    --java-mem-size=24G \
    -c -gd HUMAN \
    -nt $SLURM_NTASKS \
    -bam $OUTDIR/${sampleID}_${hgversion}.recal.bam \
    --feature-file $TMPDIR/tmp.bed

echo -e "\nEnd of workflow\n"

date
module unload tools/gatk/gatk-4.1.2.0
module unload tools/picard/2.18.5
module unload tools/samtools/samtools-1.8
module unload tools/bwa/bwa-0.7.17

rm -rf $TMPDIR
exit
