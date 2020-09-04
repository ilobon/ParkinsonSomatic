# ParkinsonSomatic

## Data processing and mapping

FASTQ files are archived in /home/devel/irenel/scratch/parkinson/BGI

#### 1.Quality conversion with trimmomatic
/home/devel/irenel/scratch/parkinson/scripts/trimmomaticConvertQ.sh

#### 2. Mapping and coverage distribution
./addRG.sh
./mergeBams.sh
./rmDup.sh
./baseRecalibration.sh
./indelRealignment.sh
./rmSecAlignments.sh
./intersectExomeRegions.sh
./coverageExome.sh

## Germline variants





