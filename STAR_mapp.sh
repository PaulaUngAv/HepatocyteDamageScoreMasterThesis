######## run STAR with WT1 RNAseq
# #STAR version 2.5

# # indeces
# STAR --runThreadN 7 \
#     --runMode genomeGenerate \
#     --genomeDir /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/STAR_indeces \
#     --genomeFastaFiles /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.dna.chromosome.1.fa \
#     --sjdbGTFfile /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.88.gtf \
#     --sjdbOverhang 61




# mapping
data_dir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/wtPodocytes/FASTA/

cd $data_dir

 ### single-end mode  
samples=$(find * -type f -not -name "*_*")
for fq in $samples;
do
  fq2=${fq/_R1/_R2}
  dirname=${fq##*/}
  dirname=${dirname/_R1_001.fastq.gz/}
  echo $dirname
  STAR --runThreadN 8 \
    --genomeDir /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/STAR_indeces \
#     --readFilesCommand zcat \
    --readFilesIn $fq $fq2 \
    --outFileNamePrefix /data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/NPHS2/STAR_mapped12w/$dirname \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMmapqUnique 42 \
    --quantMode  GeneCounts\
    --outReadsUnmapped Fastx
done

 ### paired-end mode  
samples=$(find * -type f -name "*_1.*")
for fq in $samples;
do
  fq2=${fq/_R1/_R2}
  dirname=${fq##*/}
  dirname=${dirname/_R1_001.fastq.gz/}
  echo $dirname
  STAR --runThreadN 8 \
    --genomeDir /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/STAR_indeces \
#     --readFilesCommand zcat \
    --readFilesIn $fq $fq2 \
    --outFileNamePrefix /data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/NPHS2/STAR_mapped12w/$dirname \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMmapqUnique 42 \
    --quantMode  GeneCounts\
    --outReadsUnmapped Fastx
done

# indexing and generating bigWig files
data_dir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/NPHS2/STAR_mapped

samples=$(ls $data_dir/*.bam)
for bam in $samples;
do
    # samtools index -b $bam
    file=$(basename $bam)
    bamCoverage -b $bam -o ${file/out.bam/.bw} --numberOfProcessors 8
done

# ## Count reads in genes by HTseq-count using reads mapped to exones
# #version 0.7.2
# data_dir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/NPHS2/STAR_mapped
# cd $data_dir
# samples=$(ls $data_dir/*.bam)
# dirname=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/WT1/glomeruli_Ribo_minus/HTseq_count
# cd $dirname
# for bam in $samples;
# do
#     file=$(basename $bam)
#     file=${file/Aligned.sortedByCoord.out.bam/_genecount_rexones.txt}
#     # counting in the reverse mode - the first read has to be mapped to the opposite strand as the feature and the second read to the same strand
#     htseq-count -f bam -r pos --stranded reverse $bam /data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.88.gtf > $file
#     echo $file
# done
# 
# ## Count reads in genes by HTseq-count using reads mapped to intrones 
# data_dir=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/NPHS2/STAR_mapped
# gtf=/data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/Mus_musculus.GRCm38.88.gtf
# cd $data_dir
# samples=$(ls $data_dir/*.bam)
# dirname=/data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/WT1/glomeruli_Ribo_minus/HTseq_count
# cd $dirname
# for bam in $samples;
# do
#     file=$(basename $bam)
#     file=${file/Aligned.sortedByCoord.out.bam/_genecount_rGene.txt}
#     # counting in the reverse mode - the first read has to be mapped to the opposite strand as the feature and the second read to the same strand
#     htseq-count -f bam -r pos -t gene --stranded reverse $bam $gtf > $file
#     echo $file
# done

# ## RSEQC quality check and genome coverage and various statistics
# # version 2.6.4
# bed=/data/public/tpadvits/global_data/Genome_annot/mus_musculus/GRCm38/mm10_GENCODE_VM11_basic_mod.bed
# for bam in $samples;
# do
#     file=$(basename $bam)
#     file=${file/Aligned.sortedByCoord.out.bam/_}
#     infer_experiment.py -r $bed -i $bam 
#     #geneBody_coverage.py -r $bed -i $bam  -o /data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/NPHS2/STAR_mapped/$file
#     #junction_annotation.py -i $bam -o $file -r $bed
#     #read_distribution.py -i $bam -r $bed > /data/public/tpadvits/PROJECTS/PodocytePJ/RNA_seq/NPHS2/STAR_mapped/$file"readDist.txt"
# done






    
