#!/bin/bash
#SBATCH --account=bgmp              ### This is my actual account for charging
#SBATCH --partition=bgmp            ### queue to submit to (backup option = compute)
#SBATCH --job-name=htseq-count      ### job name
#SBATCH --output=htseq_%j.out       ### file in which to store job stdout
#SBATCH --error=htseq_%j.err        ### file in which to store job stderr
#SBATCH --time=4:00:00              ### wall-clock time limit, in minutes
#SBATCH --mem=32G                   ### memory limit per node, in MB
#SBATCH --nodes=1                   ### number of nodes to use
#SBATCH --cpus-per-task=8           ### number of cores for each task

conda activate bgmp_qaa



#stranded yes
/usr/bin/time -v htseq-count --stranded=yes /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_1_2A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 1_2A_htseqout_AlignedStr.tsv

/usr/bin/time -v htseq-count --stranded=yes /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_24_4A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 24_4A_htseqout_AlignedStr.tsv

#stranded reverse
/usr/bin/time -v htseq-count --stranded=reverse /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_1_2A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 1_2A_htseqout_AlignedUnstr.tsv

/usr/bin/time -v htseq-count --stranded=reverse /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/ikes_24_4A_Aligned.out.sam /projects/bgmp/ikes/bioinfo/Bi623/QAA/Aligns/Mus_musculus.GRCm39.110.gtf > 24_4A_htseqout_AlignedUnstr.tsv