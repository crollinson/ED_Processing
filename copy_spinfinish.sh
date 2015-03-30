#!bin/bash
# This code copies the end of the spin file into the run directory & replaces the year as
#    the first day & year for the full runs 
 
spin_dir=/projectnb/dietzelab/paleon/ED_runs/phase1a_spinfinish.v2/
run_dir=/projectnb/dietzelab/paleon/ED_runs/phase1a_runs.v2/


spin_end=2000
run_start=1850
sites=(PHA PHO PUN PBL PDL PMB)

# Replacing the bad year with a good year
for site in ${sites[@]}
do
   cp ${spin_dir}${site}/histo/${site}spin-S-${spin_end}-01-01-000000-g01.h5 ${run_dir}${site}/histo/${site}spin-S-${run_start}-01-01-000000-g01.h5
done


