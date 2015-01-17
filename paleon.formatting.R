## --------------------------------------------------------------
## Script to loop through the 6 sites and extract data into the 100-yr files 
## --------------------------------------------------------------
source("/projectnb/dietzelab/paleon/ED_runs/phase1a_runs_HI2/model2netcdf.ED2.paleon.R")

raw.path <- "/projectnb/dietzelab/paleon/ED_runs/phase1a_runs_HI2/"
new.path <- "/projectnb/dietzelab/paleon/ED_runs/phase1a_runs_HI2/paleon_out/"
site.list <- c("PHA", "PHO", "PUN", "PBL", "PDL", "PMB")
site.lat <- c(42.54, 45.25, 46.22, 46.28, 47.17, 43.61)
site.lon <- c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83)

start.run <- as.Date("1850-01-01", "%Y-%m-%d") 
end.run <- as.Date("3010-01-01", "%Y-%m-%d") 


#for(s in unique(sites)){
site=""
  raw.dir <- paste0(raw.path, site, "/")
  new.dir <- paste0(new.path, site, "/")
  sitelat <- site.lat[which(site.list==site)]
  sitelon <- site.lon[which(site.list==site)]

  flist <- dir(file.path(raw.dir, "analy/"),"-E-") # Getting a list of what has been done
  
  # Getting a list of years that have been completed
  yr <- rep(NA,length(flist)) # create empty vector the same length as the file list
  for(i in 1:length(flist)){
    index <- gregexpr("-E-",flist[i])[[1]] # Searching for the monthly data marker (-E-); returns 3 bits of information: 1) capture.start (4); 2) capture.length (3); 3) capture.names (TRUE)
    index <- index[1] # indexing off of just where the monthly flag starts
    yr[i] <- as.numeric(substr(flist[i],index+3,index+6)) # putting in the Years, indexed off of where the year starts & ends
  }  
  
  start.loop <- as.Date(paste0(min(yr), "-01-01"), "%Y-%m-%d")
  end.loop <- as.Date(paste0(max(yr), "-01-01"), "%Y-%m-%d")
  bins <- c(as.numeric(strftime(start.loop, '%Y')), seq(from=as.numeric(paste0(substr(as.numeric(strftime(start.loop, "%Y"))+100, 1, 2), "00")), to=as.numeric(strftime(end.loop, '%Y')), by=100)) # Creating a vector with 100 year bins for the time period of interest

  print(paste0("----------  Processing Site: ", site, "  ----------")) 
  
  model2netcdf.ED2.paleon(site, raw.dir, new.dir, sitelat, sitelon, start.run, start.loop, end.run, end.loop, bins)
#}

