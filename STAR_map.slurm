#!/bin/bash -l
#SBATCH --account=tpadvits
#SBATCH --cpus-per-task=16
#SBATCH --error=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/snRNAseq_kif3a/STAR_kif3a.error
#SBATCH --output=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/snRNAseq_kif3a/STAR_kif3a.out
#SBATCH --workdir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/snRNAseq_kif3a
#SBATCH --job-name=STAR_kif3a
#SBATCH --partition=all
#SBATCH --ntasks=1


#### 
# module load STAR-2.6.0c


#CMD="/data/public/tpadvits/PROJECTS/PodocytePJ/Chip_seq/Pipeline/bwaANDqc_raw.sh"
bash /data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/STARsc_mapp.sh
