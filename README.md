# ParkinsonSomatic

## 1. Data processing and mapping

FASTQ files are archived in 
/home/devel/irenel/scratch/parkinson/BGI

Scripts are in 
/home/devel/irenel/scratch/parkinson/scripts/

trimmomatic-0.36.jar PE -phred64 -trimlog $log $fq1 $fq2 $out1P $out1U $out2P $out2U TOPHRED33

(trimmomaticConvertQ.sh)


(addRG.sh)

(mergeBams.sh)

(rmDup.sh)

(baseRecalibration.sh)

(indelRealignment.sh)

(rmSecAlignments.sh)

(intersectExomeRegions.sh)

(coverageExome.sh)

## 2. Germline variants





