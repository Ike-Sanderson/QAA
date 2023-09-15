#!/bin/bash
#SBATCH --account=bgmp   ### change this to your actual account for charging
#SBATCH --partition=compute       ### queue to submit to
#SBATCH --job-name=star-align    ### job name
#SBATCH --output=hostname_%j.out   ### file in which to store job stdout
#SBATCH --error=hostname_%j.err    ### file in which to store job stderr
#SBATCH --time=2:00:00                ### wall-clock time limit, in minutes
#SBATCH --mem=32G              ### memory limit per node, in MB
#SBATCH --nodes=1               ### number of nodes to use
#SBATCH --cpus-per-task=8       ### number of cores for each task

/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn Bi623/QAA/trimcontrol/1_2A_forward_paired.fq.gz /projects/bgmp/ikes/bioinfo/Bi623/QAA/trimcontrol/1_2A_reverse_paired.fq.gz \
--genomeDir /projects/bgmp/ikes/bioinfo/Bi623/QAA/mouse_GRCm39_dna_ens110_STAR_2.7.10b \
--outFileNamePrefix ikes