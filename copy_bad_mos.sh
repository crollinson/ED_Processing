#!bin/bash
#This code replaces years that are causing problems in the paleon runs with a different year
#Original, wonky met files are saved so that they can be examined in greater detail at a later point in time
#Edits: Christy Rollinson, January 2015, crollinson@gmail.com

met_dir=/projectnb/dietzelab/paleon/ED_runs/met_drivers/phase1a_met/met_v4.2/
orig_dir=/projectnb/dietzelab/paleon/ED_runs/met_drivers/phase1a_met/met_v4.2/originals/

vars=(dlwrf prate pres sh tmp ugrd vbdsf)
sites=(PHA PHO PUN PBL PDL PMB)
sites2=(PHO PUN PBL PDL PMB)

yrs_bad=(1852)
yr_good=(1853)
mos=(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC)


# make site dirs
if [ ! -d ${orig_dir} ]
then
    mkdir ${orig_dir}
fi



# # copying bad months to an archive file ("originals")
for site in ${sites2[@]}
do
	if [ ! -d ${orig_dir}${site} ]
	then
    	mkdir ${orig_dir}${site}
	fi

	for yr in ${yrs_bad[@]}
	do
		for var in ${vars[@]}
		do
			for mo in ${mos[@]}
			do
				cp ${met_dir}${site}/${var}/${site}_${var}_${yr}${mo}.h5 ${orig_dir}
			done
		done
	done
done 


# Replacing the bad year with a good year
for site in ${sites[@]}
do
	for yr in ${yrs_bad[@]}
	do
		for var in ${vars[@]}
		do
			for mo in ${mos[@]}
			do
				cp ${met_dir}${site}/${var}/${site}_${var}_${yr_good}${mo}.h5 ${met_dir}${site}/${var}/${site}_${var}_${yr}${mo}.h5
			done
		done
	done
done


