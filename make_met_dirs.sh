#!bin/bash
#This code makes the folders for new met data 
#Edits: Christy Rollinson, January 2015, crollinson@gmail.com

out_dir=/projectnb/dietzelab/paleon/ED_runs/met_drivers/phase1a_met/met_v4.1/

vars=(dlwrf prate pres sh tmp ugrd vbdsf)
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

for var in ${vars[@]}
do
if [ ! -d ${out_dir}${site}/${var} ]
then 
    mkdir ${out_dir}${site}/${var}/
fi
done
done
