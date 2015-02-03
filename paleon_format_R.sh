#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/ED_runs/phase1a_runs.v2/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=10:00:00
#$ -N Format
#cd /projectnb/dietzelab/paleon/ED_runs/phase1a_runs.v2/
R CMD BATCH paleon.formatting.R formatting.log