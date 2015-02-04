#!bin/bash
#This code makes the folders for new met data 
#Edits: Christy Rollinson, January 2015, crollinson@gmail.com

out_dir=/projectnb/dietzelab/paleon/ED_runs/met_drivers/phase1a_met/met_v4.1/

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
		pushd ${out_dir}${site}/
		
		ln -s /usr2/postdoc/crolli/ED2/ED/build/ed_2.1-opt
	fi
done
