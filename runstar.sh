#~/bin/bash

for FILE in raw/*.fastq.gz 
do
	SAMPLE=$(basename $FILE .fastq.gz)
	mkdir aligned2/$SAMPLE
	STAR \
		--genomeDir ../ref3 \
		--readFilesIn $FILE \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
		--quantMode GeneCounts \
		--runThreadN 8 \
		--alignIntronMax 1 \
		--outFileNamePrefix aligned2/$SAMPLE/ 

	samtools sort aligned2/$SAMPLE/Aligned.out.bam -o aligned2/$SAMPLE/$SAMPLE.star.sorted.bam
	rm aligned2/$SAMPLE/Aligned.out.bam
	samtools index -b -@ 8 aligned2/$SAMPLE/$SAMPLE.star.sorted.bam > aligned2/$SAMPLE/$SAMPLE.star.sorted.bam.bai
done
