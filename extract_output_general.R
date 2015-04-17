# Script to extract monthly output from ED and put it into a netcdf 


source('~/Desktop/PalEON CR/ED_PalEON_processing/ED_post-process_ncdf.R', chdir = TRUE)

site = "PHA"
lat  = 42.54
lon  = -72.18

raw.dir <- "/projectnb/dietzelab/paleon/ED_runs/phase1a_spinfinish.v2/PHA/analy"
new.dir <- "/projectnb/dietzelab/paleon/ED_runs/phase1a_spinfinishl.v2/post-process"

model2netcdf.ED2.simple("PHA", start.run=1850, raw.dir=raw.dir, new.dir=new.dir, out.name="PHA", sitelat=lat, sitelon=lon)
