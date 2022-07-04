ANNO=$1

[ ! -d bigWig ] && ( mkdir bigWig )
for GENO in $(awk '$1~"TS"{print $6}' $ANNO | sort -u);
do
    echo $GENO
    awk -v "geno=$GENO" '$6 == geno' $ANNO > tmp.txt
    for TIME in $(cut -f5 tmp.txt | sort -u);
    do
	echo "  $TIME"
	BAMS=$(awk -v "time=$TIME" 'BEGIN{ORS=" "}$5 == time{print "STAR/" $1 ".Sorted.out.bam"}' tmp.txt)
	OUT=$(echo "${GENO}.${TIME}")
	if [ ! -f STAR/$OUT.bam ]; then
	echo "    Merging replicas - $OUT"
	    samtools merge --threads 30 STAR/$OUT.bam $BAMS
	    samtools index STAR/$OUT.bam
	fi
	echo "    Making BW file"
	samtools view -b -f 128 -F 16 STAR/$OUT.bam > STAR/$OUT.fwd1.bam
	samtools view -b -f 80 STAR/$OUT.bam > STAR/$OUT.fwd2.bam
	samtools merge --threads 30 -f STAR/$OUT.fwd.bam STAR/$OUT.fwd1.bam STAR/$OUT.fwd2.bam
	samtools index STAR/$OUT.fwd.bam
	samtools view -b -f 144 STAR/$OUT.bam > STAR/$OUT.rev1.bam
	samtools view -b -f 64 -F 16 STAR/$OUT.bam > STAR/$OUT.rev2.bam
	samtools merge --threads 30 -f STAR/$OUT.rev.bam STAR/$OUT.rev1.bam STAR/$OUT.rev2.bam
	samtools index STAR/$OUT.rev.bam
	
	bamCoverage -p "max" \
		    --binSize 1 \
		    --skipNonCoveredRegions \
		    -b STAR/$OUT.rev.bam \
		    -o bigWig/$OUT.rev.bw 2> bigWig/$OUT.rev.log

	bamCoverage -p "max" \
		    --binSize 1 \
		    --skipNonCoveredRegions \
		    -b STAR/$OUT.fwd.bam \
		    -o bigWig/$OUT.fwd.bw 2> bigWig/$OUT.fwd.log

	rm STAR/$OUT.fwd*.bam
	rm STAR/$OUT.rev*.bam
    done
done


## Change chromosome name and resize to minimum coverage (across
## samples)

[ -f bigWig/sampleCoverage.txt ] && ( rm bigWig/sampleCoverage.txt )
## Count properly mapped reads:
for GENO in $(awk '$1~"TS"{print $6}' $ANNO | sort -u);
do
    echo $GENO
    awk -v "geno=$GENO" '$6 == geno' $ANNO > tmp.txt
    for TIME in $(cut -f5 tmp.txt | sort -u);
    do
	SAMPLE=${GENO}.${TIME}
	COUNT=$(samtools view -c -F 260 STAR/$SAMPLE.bam)
	echo -e "\t$TIME\t$COUNT"
	echo -e "$SAMPLE\t$COUNT" >> bigWig/sampleCoverage.txt
    done
done

awk 'BEGIN{FS=OFS="\t"}{if ( min>$2 || min==""){min=$2;}; c[$1]=$2}END{for(i in c){print i, c[i], min/c[i]}}' bigWig/sampleCoverage.txt | sort -k1,1 > bigWig/tmp
mv bigWig/tmp bigWig/sampleCoverage.txt

cd bigWig
[ ! -f mm10.chrSize ] && ( fetchChromSizes mm10 > mm10.chrSize )
for GENO in $(awk '$1~"TS"{print $6}' ../$ANNO | sort -u);
do
    echo $GENO
    awk -v "geno=$GENO" '$6 == geno' ../$ANNO > tmp.txt
    for TIME in $(cut -f5 tmp.txt | sort -u);
    do
	SAMPLE=${GENO}.${TIME}
	PER=$(grep "$SAMPLE" sampleCoverage.txt | cut -f3)
	echo $SAMPLE
	## bigWigToWig $SAMPLE /dev/stdout | awk -v "p=$PER" 'BEGIN{FS=OFS="\t"}$1 !~ "^\#"{print "chr" $1, $2, $3, int($4*p)}$1 ~ "^\#"{print}' | sed '/chr[JGM]/d' | wigToBigWig /dev/stdin mm10.chrSize $SAMPLE.tmp
	bigWigToWig $SAMPLE.fwd.bw /dev/stdout | awk -v "p=$PER" 'BEGIN{FS=OFS="\t"}$1 !~ "^\#"{print "chr" $1, $2, $3, int($4*p)}$1 ~ "^\#"{print}' 2> /dev/null | sed '/chr[JGM]/d' | wigToBigWig /dev/stdin mm10.chrSize $SAMPLE.fwd.tmp
	mv $SAMPLE.fwd.tmp $SAMPLE.fwd.bw
	###
	bigWigToWig $SAMPLE.rev.bw /dev/stdout | awk -v "p=$PER" 'BEGIN{FS=OFS="\t"}$1 !~ "^\#"{print "chr" $1, $2, $3, int($4*p)}$1 ~ "^\#"{print}' 2> /dev/null | sed '/chr[JGM]/d' | wigToBigWig /dev/stdin mm10.chrSize $SAMPLE.rev.tmp
	mv $SAMPLE.rev.tmp $SAMPLE.rev.bw
    done
done
