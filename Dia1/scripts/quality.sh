#!/usr/bin/env bash
#SBATCH --job-name=Compass
#SBATCH --ntasks=12
#SBATCH --time=06:00:00
#SBATCH --output=fastqc.%j.out

module load tools/fastqc/fastqc_v0.11.7

fastq_dir=/home/marcopretti/ligandoma/LBB/data
nt=12
outdir=$fastq_dir

mkdir $fastq_dir/quality

date
echo "Quality check of reads" 
for i in $(ls $fastq_dir/*R?.fq.gz); do
	fastqc -t $nt $i -o $outdir/quality ;
done;

module unload tools/fastqc/fastqc_v0.11.7

exit
