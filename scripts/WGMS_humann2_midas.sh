#! /bin/bash
#
#$ -S /bin/bash
#$ -o humann2_log
#$ -e humann2_log
#$ -r y
#$ -j y
#$ -R yes
#$ -pe smp 4
#$ -t 1-19
#$ -cwd

tasks1=(0	Q.71401.0103.2017.04.05.v2_R1.fastq.gz	Q.71401.0162.2017.11.14_R1.fastq.gz	Q.71401.0186.2018.04.10_R1.fastq.gz	Q.71401.0245.2018.09.08_R1.fastq.gz	Q.71601.0131.2018.04.01_R1.fastq.gz	Q.71601.0149.2018.07.26_R1.fastq.gz	Q.71701.0127.2017.10.23_R1.fastq.gz	Q.71702.0152.2018.01.30_R1.fastq.gz	Q.71801.0110.2017.10.09_R1.fastq.gz	Q.71802.0130.2018.01.10_R1.fastq.gz	Q1.71401.0002.12.27.15_R1.fastq.gz	Q1.71702.0015.8.1.16_R1.fastq.gz	Q2.71402.0073.11.29.16_R1.fastq.gz	Q2.71601.0010.8.10.16_R1.fastq.gz	Q2.71801.0030.9.29.16_R1.fastq.gz	S.71701.0043.10.27.16_R1.fastq.gz	S.71802.0004.4.17.16_R1.fastq.gz	S.71802.0013.6.22.16_R1.fastq.gz	UCSF.12_R1.fastq.gz)
tasks2=(0	Q.71401.0103.2017.04.05.v2_R2.fastq.gz	Q.71401.0162.2017.11.14_R2.fastq.gz	Q.71401.0186.2018.04.10_R2.fastq.gz	Q.71401.0245.2018.09.08_R2.fastq.gz	Q.71601.0131.2018.04.01_R2.fastq.gz	Q.71601.0149.2018.07.26_R2.fastq.gz	Q.71701.0127.2017.10.23_R2.fastq.gz	Q.71702.0152.2018.01.30_R2.fastq.gz	Q.71801.0110.2017.10.09_R2.fastq.gz	Q.71802.0130.2018.01.10_R2.fastq.gz	Q1.71401.0002.12.27.15_R2.fastq.gz	Q1.71702.0015.8.1.16_R2.fastq.gz	Q2.71402.0073.11.29.16_R2.fastq.gz	Q2.71601.0010.8.10.16_R2.fastq.gz	Q2.71801.0030.9.29.16_R2.fastq.gz	S.71701.0043.10.27.16_R2.fastq.gz	S.71802.0004.4.17.16_R2.fastq.gz	S.71802.0013.6.22.16_R2.fastq.gz	UCSF.12_R2.fastq.gz)
input1="${tasks1[$SGE_TASK_ID]}"
input2="${tasks2[$SGE_TASK_ID]}"
index1=$(sed s/.fastq.gz//g <<< ${input1})
index2=$(sed s/.fastq.gz//g <<< ${input2})
folder=$(sed s/_R1.*fastq.gz//g <<< ${input1})

date
hostname
echo R1=${input1} R2=${input2}

SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID
mkdir -p $SCRATCH_JOB

echo "Copy fastq files"
cp -r metagenomics_samples_raw_merge/${folder}* $SCRATCH_JOB/
cd $SCRATCH_JOB
mkdir $folder

echo "Combine forward and reverse reads"
cat ${input1} ${input2} > ${folder}.fastq.gz


echo "Run microbecensus"
run_microbe_census.py ${input1},${input2} ${folder}/${folder}\_microbcensus

echo "Run humann2"
source activate shogun
humann2 --input ${folder}.fastq.gz --output $folder --nucleotide-database databases/chocophlan --protein-database databases/uniref > $folder/${folder}.log

echo "Normalize the abundance"
if [ -f $folder/$folder\_genefamilies.tsv ] 
then
humann2_renorm_table --input $folder/$folder\_genefamilies.tsv --output $folder/$folder\_genefamilies-relab.tsv --units relab
humann2_renorm_table --input $folder/$folder\_pathabundance.tsv --output $folder/$folder\_pathabundance-relab.tsv --unite relab
fi

#echo "Run MIDAS"
#mkdir $folder\_midas
#run_midas.py species $folder\_midas -1 ${input1} -2 ${input2}
#run_midas.py genes $folder\_midas -1 ${input1} -2 ${input2}
#run_midas.py snps $folder\_midas -1 ${input1} -2 ${input2}

echo "Move output back to scratch"
mv -f $SCRATCH_JOB/${folder} humann2/
#mv -f $SCRATCH_JOB/$folder\_midas midas/

echo "Remove job specific scratch folder"
cd /scratch/
rm -rf $SCRATCH_JOB/
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
