
#Convert fastqs PHRED format 
java -jar trimmomatic-0.36.jar PE -phred64 -trimlog $log $fq1 $fq2 $fq1P $fq1U $fq2P $fq2U TOPHRED33 

#BWA mapping (to the hs37d5 reference)
bwa mem -M -t 4 $ref $fq1P $fq2P | samtools sort -o $sorted.bam

#Add read groups
java -jar PICARD/1.95/AddOrReplaceReadGroups.jar INPUT=$sorted.bam OUTPUT=$RG.bam RGID=${sample}.${lane} RGLB=${sample} RGPL=Illumina RGPU=${lane} RGSM=${sample}

#Merge different lanes
samtools merge -f -@ 8 $merged.bam $bamFiles

#Remove duplicates
java -jar PICARD/1.95/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_RECORDS_IN_RAM=1500000  METRICS_FILE=$metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$merged.bam COMPRESSION_LEVEL=9 OUTPUT=$rmDup.bam

#GATK base quality recalibration
java -jar GATK/3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R $ref -I $rmDup.bam -knownSites $dbsnp -o $recalibration_data.table
java -jar GATK/3.6/GenomeAnalysisTK.jar -T PrintReads -R $ref -I $rmDup.bam -BQSR $recalibration_data.table -o $recal.bam

#GATK realignment
java -jar GATK/3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I $recal.bam --known $indels -o $realign1.intervals
java -jar GATK/3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I $recal.bam --known $indels2 -o $realign2.intervals
java -jar GATK/3.6/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $recal.bam -known $indels -targetIntervals $realign1.intervals -o $realgn1.bam
java -jar GATK/3.6/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $realgn1.bam -known $indels2 -targetIntervals $realign2.intervals -o $realgn.bam

#Remove secondary alignments
samtools view -b -F 256 $realgn.bam > $secalgn.bam
samtools index $secalgn.bam

#Intersect with exome target space
BEDTools/2.26.0/bin/intersectBed -abam $secalgn.bam -b exome_target.bed $exome.bam
samtools index $exome.bam

#Convert processed bam to cram
samtools view -T $ref -hC -o secalgn.cram secalgn.bam
