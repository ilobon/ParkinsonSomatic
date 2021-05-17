
#Make GVCFs
java -jar GATK/3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I $exome.bam -ploidy 2 -A StrandAlleleCountsBySample --emitRefConfidence GVCF -o $exome.g.vcf

#Genotype GVCFs
java -jar GATK/3.6/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $ref $allGVCFfiles -ploidy 2 -A StrandAlleleCountsBySample -o allSamples.p2.vcf

#GATK hard filter
java -jar GATK/3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $ref $allGVCFfiles -V allSamples.p2.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "HardFilter" -o allSamples.p2.hardFiltered.vcf
bgzip allSamples.p2.hardFiltered.vcf
tabix allSamples.p2.hardFiltered.vcf.gz

#Intersect with OMIM genes
intersectBed -a allSamples.p2.hardFiltered.vcf.gz -b parkinsonOMIM.bed >> allSamples.p2.hardFiltered.OMIMPDgenes.vcf

#Intersect with GWAS variants
intersectBed -a allSamples.p2.hardFiltered.vcf.gz -b parkinsonGWAS.bed >> allSamples.p2.hardFiltered.GWASPDgenes.vcf

#Annotate variants
java -jar snpEff/snpEff.jar eff -v -i vcf hg19 allSamples.p2.hardFiltered.vcf.gz > allSamples.p2.hardFiltered.snpeff.vcf
java -jar snpEff/SnpSift.jar dbnsfp -db snpEff/db/dbNSFP2.9.txt.gz allSamples.p2.hardFiltered.snpeff.vcf > allSamples.p2.hardFiltered.snpeff.snpsift.vcf
java -jar GATK/3.6/GenomeAnalysisTK.jar -T VariantFiltration -R $ref --variant allSamples.p2.hardFiltered.snpeff.snpsift.vcf --filterExpression "dbNSFP_CADD_phred > 15" --filterName CADD15 -o allSamples.p2.hardFiltered.snpeff.snpsift.CADD15.vcf

#Deleterious variants (CADD>15 & SIFT=="D" & 1000G_EUR_AF <0.1) (+not if only variant in DV2)
grep CADD15 DV.p2.hardFiltered.snpeff.snpsift.CADD15.vcf | awk '{split($8,a,"dbNSFP_SIFT_pred=");split(a[2],b,";");split($8,c,"dbNSFP_1000Gp1_EUR_AF=");split(c[2],d,";");if(b[1]~"D" && d[1]<0.1){print $0}}' | awk '{var=0;for(i=10;i<=NF;i=i+5){split($i,a,":");if(a[1]!="0/0"){var+=1}}split($20,a,":");if(var>0 && (a[1]=="0/0" || var>1))print $0}' > filteredVariantsNoDV2.vcf

#Get genenames for enrichment
awk '{split($8,e,"|");print e[4]}' filteredVariantsNoDV2.vcf | sort -V | uniq > deleteriousGenes.txt


##PCA

#Make plink files
vcftools --gzvcf allSamples.p2.hardFiltered.vcf.gz --plink --out allSamples.p2.hardFiltered

#Make .ind
paste <(awk '{print $1}' allSamples.p2.hardFiltered.ped) <(printf '0\n%.0s' {1..50}) <(awk '{print $1}' allSamples.p2.hardFiltered.ped | sed s'/.$//' ) > allSamples.p2.hardFiltered.ind

#Make .par
cp template.par .
sed s/IDname/allSamples.p2.hardFiltered/g template.par | head -n -2 > allSamples.p2.hardFiltered.par

#Run PCA
EIGENSOFT/7.2.1/bin/smartpca -p allSamples.p2.hardFiltered.par

#Plot PCA
./3.plot_germline_PCA.R


#Genes per ind for enrichment plot
for z in 10 15 $(seq 25 5 55);do while read line;do echo $z $line;done< <(awk -v z=$z '{split($z,a,":");if(a[1]!="0/0"){split($8,e,"|");print e[4]}}' filteredVariantsNoDV2.vcf) >> perInd.deleteriousGenes.txt;done

