#! /bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/scratch/xiaoyuan/humann2_log
#$ -e /wynton/scratch/xiaoyuan/humann2_log
#$ -r y
#$ -j y
#$ -l scratch=300G
#$ -l mem_free=16G
#$ -R yes
#$ -l h_rt=100:00:00
#$ -pe smp 4
#$ -t 1-9
#$ -cwd
tasks1=(0	71802.0066_R1.fastq.gz	71802.0070_R1.fastq.gz	71701.0009_R1.fastq.gz	71601.0061_R1.fastq.gz	71601.0141_R1.fastq.gz	71602.0124_R1.fastq.gz	71801.0142_R1.fastq.gz	71802.0084_R1.fastq.gz	71802.0125_R1.fastq.gz)
tasks2=(0	71802.0066_R2.fastq.gz	71802.0070_R2.fastq.gz	71701.0009_R2.fastq.gz	71601.0061_R2.fastq.gz	71601.0141_R2.fastq.gz	71602.0124_R2.fastq.gz	71801.0142_R2.fastq.gz	71802.0084_R2.fastq.gz	71802.0125_R2.fastq.gz)
input1="${tasks1[$SGE_TASK_ID]}"
input2="${tasks2[$SGE_TASK_ID]}"
index1=$(sed s/.fastq.gz//g <<< ${input1})
index2=$(sed s/.fastq.gz//g <<< ${input2})
folder=$(sed s/_R1.fastq.gz//g <<< ${input1})

date
hostname
echo R1=${input1} R2=${input2}

SCRATCH_JOB=/scratch/$USER/jobs/$JOB_ID
mkdir -p $SCRATCH_JOB

echo "Copy fastq files"
cp -r /wynton/scratch/xiaoyuan/metagenomics_subject/${input1} $SCRATCH_JOB/
cp -r /wynton/scratch/xiaoyuan/metagenomics_subject/${input2} $SCRATCH_JOB/
cd $SCRATCH_JOB

Source activate shogun
mkdir $folder\_midas
run_midas.py species $folder\_midas -1 ${input1} -2 ${input2} 
run_midas.py genes $folder\_midas -1 ${input1} -2 ${input2}  
run_midas.py snps $folder\_midas -1 ${input1} -2 ${input2}

mv -f $SCRATCH_JOB/$folder\_midas /wynton/scratch/xiaoyuan/midas_subject/

echo "Remove job specific scratch folder"
cd /scratch/
rm -rf $SCRATCH_JOB/
rmdir /scratch/$USER/jobs/$JOB_ID
rmdir /scratch/$USER/jobs
