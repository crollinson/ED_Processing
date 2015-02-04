#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/ED_runs/phase1a_spininitial.v2/PBL
#$ -j y 
#$ -S /bin/bash         
#$ -V 
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=24:00:00
#$ -N PBLspin
#cd /projectnb/dietzelab/paleon/ED_runs/phase1a_spininitial.v2/PBL
./ed_2.1-opt
