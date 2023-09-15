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

/usr/bin/time -v STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /projects/bgmp/ikes/bioinfo/Bi623/QAA/mouse_GRCm39_dna_ens110_STAR_2.7.10b \
--genomeFastaFiles /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf


#  -v STAR --runThreadN 8 --runMode alignReads \
# --outFilterMultimapNmax 3 \
# --outSAMunmapped Within KeepPairs \
# --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
# --readFilesCommand zcat \
# --readFilesIn /projects/bgmp/shared/Bi621/dre_WT_ovar12_R1.qtrim.fq.gz /projects/bgmp/shared/Bi621/dre_WT_ovar12_R2.qtrim.fq.gz \
# --genomeDir /projects/bgmp/ikes/bioinformatics/Bi621/PS/PS8/Danio_rerio.GRCz11.dna.ens109.STAR_2.7.10b \
# --outFileNamePrefix ikes