#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/ED_runs/phase1a_runs.v2/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=24:00:00
#$ -N PostProcess
#cd /projectnb/dietzelab/paleon/ED_runs/phase1a_runs.v2/

out_dir=/projectnb/dietzelab/paleon/ED_runs/phase1a_runs.v2/post-process/
sites=(PHA PHO PUN PBL PDL PMB)

# make site dirs
if [ ! -d ${out_dir} ]
then
    mkdir ${out_dir}
fi

for site in ${sites[@]}
do
if [ ! -d ${out_dir}${site} ]
then
    mkdir ${out_dir}${site}
fi
done


R CMD BATCH paleon.formatting.R formatting.log