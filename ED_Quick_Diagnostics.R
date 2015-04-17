library(ncdf4)

setwd("/projectnb/dietzelab/paleon/ED_runs")
outdir <- "pha_pft_diagnostics"
pfts <- c(6,8,9,10,11)

vars <- c("AGB_PY", "DBH", "MMEAN_GPP_PY", "MMEAN_NPP_PY", "MMEAN_NEP_PY")
soil.depth <- c(2.33, 0.67, 0.4, 0.3, 0.2, 0.15, 0.15, 0.1, 0.8, 0.6, 0.6)

#for(p in 1:length(pfts)){
	pft.dir <- file.path(paste0("PHA_MsTMIP_met"), "analy")
	pft.files <- dir(pft.dir, "-E-")	
	
	pft <- pfts[p]
	for(i in 1:length(pft.files)){
		if(i == 1) year <- month <- dbh.max <- dbh.mean <- gpp <- npp <- nep <- nplant <- lai <- hite.mean <- hite.max <- elongf <- cbr <- soil_water <- bleaf <- broot <- bswood <- bhwood <- bstorage <- mort.age <- mort.cbr <- mort.frost <- mort.disturb <- vector()
		agb <- data.frame()
		nc <- nc_open(file.path(pft.dir, pft.files[i]))
		if(pfts[p] < 10){
			year <- c(year, substr(pft.files[i], 12, 15))
			month <- c(month, substr(pft.files[i], 17, 18))
		} else {
			year <- c(year, substr(pft.files[i], 13, 16))
			month <- c(month, substr(pft.files[i], 18, 19))
		}
		nplant.now <- ncvar_get(nc, "NPLANT")

		agb <- c(agb, sum(ncvar_get(nc, "AGB_PY")))
		dbh.max <- c(dbh.max, max(ncvar_get(nc, "DBH")))
		dbh.mean <- c(dbh.mean, mean(ncvar_get(nc, "HITE")))
		hite.max <- c(hite.max, max(ncvar_get(nc, "HITE")))
		hite.mean <- c(hite.mean, mean(ncvar_get(nc, "DBH")))
		lai <- c(lai, sum(ncvar_get(nc, "MMEAN_LAI_PY")))
		elongf <- c(elongf, mean(ncvar_get(nc, "ELONGF")*nplant.now/sum(nplant.now)))
		gpp <- c(gpp, ncvar_get(nc, "MMEAN_GPP_PY"))
		npp <- c(npp, ncvar_get(nc, "MMEAN_NPP_PY"))
		nep <- c(nep, ncvar_get(nc, "MMEAN_NEP_PY"))
		nplant <- c(nplant, sum(nplant.now))
		cbr <- c(cbr, sum(ncvar_get(nc, "CBR_BAR")*nplant.now/sum(nplant.now)))
		soil_water <- c(soil_water, sum(ncvar_get(nc, "MMEAN_SOIL_WATER_PY")*soil.depth/4.5))		
		bleaf <- c(bleaf, sum(ncvar_get(nc, "MMEAN_BLEAF_PY")))
		broot <- c(broot, sum(ncvar_get(nc, "MMEAN_BROOT_PY")))
		bstorage <- c(bstorage, sum(ncvar_get(nc, "MMEAN_BSTORAGE_PY")))
		bswood <- c(bswood, sum(ncvar_get(nc, "BSAPWOODA_PY"))+sum(ncvar_get(nc, "BSAPWOODB_PY")))
		bhwood <- c(bhwood, sum(ncvar_get(nc, "BDEAD_PY")))
		mort.age <- c(mort.age, mean(ncvar_get(nc, "MORT_RATE_CO")[1,]))
		mort.cbr <- c(mort.cbr, sum(ncvar_get(nc, "MORT_RATE_CO")[2,]*nplant.now/sum(nplant.now)))
		mort.frost <- c(mort.frost, sum(ncvar_get(nc, "MORT_RATE_CO")[4,]*nplant.now/sum(nplant.now)))
		mort.disturb <- c(mort.disturb, sum(ncvar_get(nc, "MORT_RATE_CO")[5,]*nplant.now/sum(nplant.now)))

		nc_close(nc)
	}
	year.month <- paste(year, month, sep="-")
	write.csv(cbind(pft, year, month, year.month, agb, dbh.max, dbh.mean, hite.max, hite.mean, nplant, lai, elongf, gpp, npp, nep, cbr, soil_water, bleaf, broot, bstorage, bswood, bhwood, mort.age, mort.cbr, mort.frost, mort.disturb), file.path(outdir, paste0("Quick_Diagnostics_PFT", pfts[p], ".csv")), row.names=F)
}
