#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/ED_runs/ED_Diagnostic_Scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=24:00:00
#$ -N ED_Diagnostics
#cd /projectnb/dietzelab/paleon/ED_runs/ED_Diagnostic_Scripts/
R CMD BATCH extract_output_general.R