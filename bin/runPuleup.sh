#!/bin/bash

ANNO=$1
COORD=$2

echo "### Sorting BAM files"
for SAMPLE in $(cut -f2 $ANNO | sort -u);
do
    echo $SAMPLE
    INFILE="STAR/${SAMPLE}.Aligned.out.bam"
    OUTFILE="STAR/${SAMPLE}.Sorted.out.bam"
    if [ ! -f $OUTFILE ];
    then
	samtools sort -@ 20 -o $OUTFILE $INFILE
	samtools index $OUTFILE
    else
	echo '  already done, skipping'
    fi
done

echo "### making pileups"
for SAMPLE in $(cut -f2 $ANNO | sort -u);
do
    echo $SAMPLE
    for GENE in $(awk '{print $1}' $COORD);
    do
	REGION=$(awk -v "gene=$GENE" '$1 == gene{print $2}' 2021-10-29_geneCoord)
	NAME=$(awk -v "gene=$GENE" '$1 == gene{print $3}' 2021-10-29_geneCoord)
	STRAND=$(awk -v "gene=$GENE" '$1 == gene{print $4}' 2021-10-29_geneCoord)
	INFILE=tmp/$SAMPLE.$NAME.mergedSorted.bam
	OUTFILE="pileup/timeCourse/${SAMPLE}.${GENE}.${NAME}.mpileup"
	## 1. Get the properly mapped reads:
	## if gene on the minus strand
	if [ $STRAND == "-" ]; then
	    samtools view -b -f 147 STAR/$SAMPLE.Sorted.out.bam $REGION > tmp/$SAMPLE.1.bam
	    samtools view -b -f 99 STAR/$SAMPLE.Sorted.out.bam $REGION > tmp/$SAMPLE.2.bam
	    ## if gene on the plus strand
	else
	    samtools view -b -f 83 STAR/$SAMPLE.Sorted.out.bam $REGION > tmp/$SAMPLE.1.bam
	    samtools view -b -f 163 STAR/$SAMPLE.Sorted.out.bam $REGION > tmp/$SAMPLE.2.bam
	fi
	## 2. marge and sort them
	samtools merge -f tmp/$SAMPLE.merged.bam tmp/$SAMPLE.1.bam tmp/$SAMPLE.2.bam
	## 3. sort and index it
	samtools sort -@ 20 -o $INFILE tmp/$SAMPLE.merged.bam
	samtools index $INFILE
	## 4. get the depth
	samtools depth -a -r $REGION  $INFILE > $OUTFILE
    done
done
