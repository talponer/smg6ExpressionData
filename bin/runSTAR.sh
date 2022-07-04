#!/bin/bash

TXT=$1

## Make the genome
# STAR \
#     --runThreadN 30 \
#     --runMode genomeGenerate \
#     --genomeDir /mnt/data/databases/mouse/star/Mmusculus.GRCm38.91 \
#     --genomeFastaFiles /mnt/data/databases/mouse/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa \
#     --sjdbGTFfile /mnt/data/databases/mouse/gtf/Mmusculus.GRCm38.91.gtf \
#     --sjdbOverhang 100

for SAMPLE in $(cut -f2 $TXT | sort -u);
do
    FQS1=raw_data/${SAMPLE}r3_*_R1.fastq.gz
    FQFILE1=`echo ${FQS1} | sed 's/ /,/g'`

    FQS2=raw_data/${SAMPLE}r3_*_R2.fastq.gz
    FQFILE2=`echo ${FQS2} | sed 's/ /,/g'`
    
    OUTFILE="STAR/${SAMPLE}."
    STDERRFILE="logs/${SAMPLE}.STAR.stderr"
    if [[ ! -e "${OUTFILE}.Aligned.out.bam" ]]
    then
	echo -n "aligning $FQFILE1 $FQFILE2..."
	STAR \
	    --runMode alignReads \
	    --runThreadN 20 \
	    --genomeDir $GENOME \
	    --readFilesCommand gzip -dc \
	    --readFilesIn $FQFILE1 $FQFILE2 \
	    --outFilterType BySJout \
	    --outFilterMultimapNmax 20 \
	    --outMultimapperOrder Random \
	    --alignSJoverhangMin 8 \
	    --alignSJDBoverhangMin 1 \
	    --outFilterMismatchNmax 999 \
	    --alignIntronMin 20 \
	    --alignIntronMax 1000000 \
	    --alignMatesGapMax 1000000 \
	    --outFileNamePrefix $OUTFILE \
	    --outSAMmultNmax 1 \
##	    --outSAMtype BAM SortedByCoordinate \
	    --outSAMtype BAM Unsorted \
	    --outSAMunmapped None \
	    --outSAMstrandField intronMotif \
	    --outBAMsortingThreadN 10 \
	    --outFilterIntronMotifs RemoveNoncanonical \
	    --limitBAMsortRAM 40000000000 2> $STDERRFILE
	echo "done"
    else
	echo "skipping existing $OUTFILE"
    fi
done
