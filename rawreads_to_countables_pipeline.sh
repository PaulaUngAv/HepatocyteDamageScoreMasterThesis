## Extraction with SRA toolkit

# Download SRA run in SRA format
#./prefetch SRAnumber

# or on whole list of samples that belong to one accession
# Download SRA run in SRA format
#./prefetch --option-file SraAccList.txt

# check for integrity: should report 'ok' and 'consistent' for all parameters
#./vdb-validate SRAnumber

# convert from sra to fastq
#./fasterq-dump SRAnumber

#questions: 1. what does the output on the console mean: exp.
#[pungerav@beyer-ws01 bin]$ ./fasterq-dump SRR5572772
#spots read      : 17,993,089
#reads read      : 35,986,178
#reads written   : 35,986,178
#
# any options I should consider adding?

## Quality Control with FastQC

#./fastqc path_containing_output_from_fasterq_dump.fastq

# output in html form stored in directory containing fastq file
# cd to directory all the fastqc output files
# summarize results from all samples with multiqc

#multiqc .

# optional trimming and filtering step

## alignment with STAR

#STAR command line has the following format:
#STAR --option1-name option1-value(s)--option2-name option2-value(s) ...
#If an option can accept multiple values, they are separated by spaces,
#and in a few cases - by commas

#Two step process:
#1)generate genome index files: supply reference genome fasta files
# and gtf annotations file
# For GSE99010 using GRCm39 reference genome:
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/
# download corresponding annotations (beware of chromsome names)
# ftp://ftp.ncbi.nlm.nih.gov/genomes//all/annotation_releases/10090/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
#STAR --runMode genomeGenerate --genomeFastaFiles thefastfile1.fasta the fastafileN.fasta \
 #--sjdbGTFfile annotations_file.gtf --genomeDir GenomeDirectory/

#Example:
STAR --runMode genomeGenerate --genomeFastaFiles \
/data/public/pungerav/global_data/genomes/mus_musculus/GRCm39/ncbi-genomes-2021-12-01/GCF_000001635.27_GRCm39_genomic.fasta \
--sjdbGTFfile /data/public/pungerav/global_data/genomes/mus_musculus/GRCm39/ncbi-genomes-2021-12-01/GCF_000001635.27_GRCm39_genomic.gtf \
--genomeDir /data/public/pungerav/global_data/genomes/mus_musculus/GRCm39/ncbi-genomes-2021-12-01/index/ \
--sjdbOverhang 100 \
--runThreadN 20


# Error: lack of memory?
#Dec 02 11:40:23 ..... started STAR run
#Dec 02 11:40:23 ... starting to generate Genome files
#Dec 02 11:41:20 ... starting to sort Suffix Array. This may take a long time...
#Dec 02 11:41:34 ... sorting Suffix Array chunks and saving them to disk...
#terminate called after throwing an instance of 'std::bad_alloc'
#  what():  std::bad_alloc
#Aborted (Speicherabzug geschrieben)
#


#2) Mapping
#STAR --genomeDir path/to/generated_indices --readFilesIn path/to/files/tobemapped \
#--quantMode GeneCounts --outTmpKeep None --runThreadN 20 --outSAMtype None
#when running on server:
#option --runThreadN defines number of threads to be used for genome generation
# depends on number of available cores
