#RNAseq DATA analysis Pipe line
SRA Toolkit
      ↓
FastQC
      ↓
fastp
      ↓
HISAT2
      ↓
SAMtools
      ↓
featureCounts
      ↓
DESeq2


# download sra_list.txt and metadata.csv from NCBI-SRA database 
BioProject number : BioProject PRJNA667177

#Command : prefetch --progress --option-file SRR_Acc_List_mbcd.txt


#convert sra into fastq files by using 
: for srr in $(cat SRR_Acc_List_mbcd.txt)
do
    echo "Converting $srr"

    fasterq-dump $srr \
    --split-files \
    --threads 8 \
    --progress \
    -O .
done


#quality check after compression by 
: fastqc -t 4 *.fastq.gz 


#trimming  by using fastp
: for srr in $(cat SRR_Acc_List_mbcd.txt)
do
    fastp \
    -i ${srr}_1.fastq.gz \
    -I ${srr}_2.fastq.gz \
    -o trimmed/${srr}_trimmed_1.fastq.gz \
    -O trimmed/${srr}_trimmed_2.fastq.gz \
    -h trimmed/${srr}.html \
    -j trimmed/${srr}.json \
    -w 4
done






#Download reference genome and its annotation

: wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz

: wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz


#HISAT indexing for reference
:wget -c https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz


#Alignment by HISAT
: for srr in $(cat SRR_Acc_List_mbcd.txt)
do
    echo "Aligning $srr"

    hisat2 -p 4 \
    -x genome/grch38/genome \
    -1 trimmed/${srr}_trimmed_1.fastq.gz \
    -2 trimmed/${srr}_trimmed_2.fastq.gz \
    2> alignment/${srr}.log | \
    samtools sort -@ 4 -m 1G -o alignment/${srr}.sorted.bam

    samtools index alignment/${srr}.sorted.bam
done



# will give aligned .bam file



#Feature count 
; featureCounts -T 4 -p --countReadPairs \
-a reference/annotation.gtf \
-o counts/all_samples_counts.txt \
alignment/*.sorted.bam


#Output will 
:counts/all_samples_counts.txt
counts/all_samples_counts.txt.summary




#Next step DESeq2 in R



