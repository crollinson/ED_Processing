#Checking Annual and Monthly CO2 met driver patterns
library(ncdf4)

# -----------------------------------------
# Raw Env Drivers
# -----------------------------------------
env.driver <- "~/Dropbox/PalEON CR/phase1a_env_drivers/phase1a_env_drivers_v4/paleon_co2"


co2.ann <- nc_open(file.path(env.driver, "paleon_annual_co2.nc"))
co2.ann.dat <- ncvar_get(co2.ann, "co2")
summary(co2.ann$var)
nc_close(co2.ann)

plot(co2.ann.dat, type="l", lwd=3)

co2.mo <- nc_open(file.path(eng.driver, "paleon_monthly_co2.nc"))
co2.mo.dat <- ncvar_get(co2.mo, "co2")
summary(co2.mo$var)
nc_close(co2.mo)

plot(co2.mo.dat, type="l")
plot(co2.mo.dat[1:240], type="l")

length(co2.mo.dat)
length(co2.ann.dat)


# -----------------------------------------
ed.driver <- "~/Dropbox/PalEON CR/ED_Misc/spininit_INIT0/ed_CO2_driver"

ed.dir<- dir(ed.driver)

ed.co2 <- vector()
years <- 1850:1859
months <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

for(y in 1:length(years)){
	for(m in 1:length(months)){
		nc <- nc_open(file.path(ed.driver, paste0("co2_", years[y], months[m],".h5")))
		ed.co2 <- c(ed.co2, ncvar_get(nc, "co2"))
		nc_close(nc)
	}
}

summary(ed.co2)
plot(ed.co2, type='l', ylim=c(279:280))
plot(co2.mo.dat[1:120], type="l", ylim=c(279:280))


# -----------------------------------------
# ED Met Drivers
# -----------------------------------------
ed.driver <- "~/Dropbox/PalEON CR/ED_Misc/spininit_INIT0/ed_met_drivers"

ed.dlwrf <- ed.prate <- ed.pres <- ed.tmp <- ed.sh <- ed.ugrd <- ed.vbdsf <- vector()
years <- 1850:1859
months <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

for(y in 1:length(years)){
	for(m in 1:length(months)){
		dlwrf <- nc_open(file.path(ed.driver, "dlwrf", paste0("PHA_dlwrf_", years[y], months[m],".h5")))
		prate <- nc_open(file.path(ed.driver, "prate", paste0("PHA_prate_", years[y], months[m],".h5")))
		pres  <- nc_open(file.path(ed.driver, "pres",  paste0("PHA_pres_", years[y], months[m],".h5")))
		sh    <- nc_open(file.path(ed.driver, "sh",    paste0("PHA_sh_", years[y], months[m],".h5")))
		tmp   <- nc_open(file.path(ed.driver, "tmp",   paste0("PHA_tmp_", years[y], months[m],".h5")))
		ugrd  <- nc_open(file.path(ed.driver, "ugrd",  paste0("PHA_ugrd_", years[y], months[m],".h5")))
		vbdsf <- nc_open(file.path(ed.driver, "vbdsf", paste0("PHA_vbdsf_", years[y], months[m],".h5")))

		ed.dlwrf <- c(ed.dlwrf, ncvar_get(dlwrf, "dlwrf"))
		ed.prate <- c(ed.prate, ncvar_get(prate, "prate"))
		ed.pres  <- c(ed.pres,  ncvar_get(pres, "pres"))
		ed.tmp   <- c(ed.tmp,   ncvar_get(sh, "sh"))
		ed.sh    <- c(ed.sh,    ncvar_get(tmp, "tmp"))
		ed.ugrd  <- c(ed.ugrd,  ncvar_get(ugrd, "ugrd"))
		ed.vbdsf <- c(ed.vbdsf, ncvar_get(vbdsf, "vbdsf"))

		nc_close(dlwrf); nc_close(prate); nc_close(pres); nc_close(sh); nc_close(tmp); nc_close(ugrd); nc_close(vbdsf)
	}
}

summary(ed.tmp)

plot(ed.dlwrf, type='l')
plot(ed.dlwrf[1:(30*4)], type='l')
plot(ed.dlwrf[1:24], type='l')

plot(ed.prate, type='l')
plot(ed.prate[1:(30*4)], type='l')

plot(ed.pres, type='l')
plot(ed.pres[1:(30*4)], type='l')

plot(ed.sh, type='l')
plot(ed.sh[1:(30*4)], type='l')

plot(ed.tmp, type='l')
plot(ed.tmp[1:(30*4)], type='l')

plot(ed.ugrd, type='l')
plot(ed.ugrd[1:(30*4)], type='l')

plot(ed.vbdsf, type='l')
plot(ed.vbdsf[1:(30*4)], type='l')
plot(ed.vbdsf[1:24], type='l')


# -----------------------------------------
# MsTMIP Met Drivers
# -----------------------------------------
mstmip.driver <- "~/Dropbox/PalEON CR/ED_Misc/spininit_INIT0/MsTMIP_met_driver"

years <- 1991:2006
months <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

m.dlwrf <- m.prate <- m.pres <- m.sh <- m.tmp <- m.ugrd <- m.vbdsf<- m.vddsf <- m.nbdsf <- m.nddsf <- vector()

for(y in 1:length(years)){
	for(m in 1:length(months)){
		mstmip <- nc_open(file.path(mstmip.driver, paste0("US-Ha1_", years[y], months[m],".h5")))

		m.dlwrf <- c(m.dlwrf, ncvar_get(mstmip, "dlwrf"))
		m.prate <- c(m.prate, ncvar_get(mstmip, "prate"))
		m.pres  <- c(m.pres,  ncvar_get(mstmip, "pres"))
		m.tmp   <- c(m.tmp,   ncvar_get(mstmip, "sh"))
		m.sh    <- c(m.sh,    ncvar_get(mstmip, "tmp"))
		m.ugrd  <- c(m.ugrd,  ncvar_get(mstmip, "ugrd"))
		m.vbdsf <- c(m.vbdsf, ncvar_get(mstmip, "vbdsf"))
		m.vddsf <- c(m.vddsf, ncvar_get(mstmip, "vddsf"))
		m.vddsf <- c(m.vddsf, ncvar_get(mstmip, "vddsf"))
		m.nbdsf <- c(m.nbdsf, ncvar_get(mstmip, "nbdsf"))
		m.nddsf <- c(m.nddsf, ncvar_get(mstmip, "nddsf"))

		nc_close(mstmip)
	}
}

time.step <- seq(from=1, to=length(m.dlwrf), by=6)
length(time.step)

m2.dlwrf <- m2.prate <- m2.pres <- m2.sh <- m2.tmp <- m2.ugrd <- m2.vbdsf <- m2.vddsf <- m2.nbdsf <- m2.nddsf <- vector()
for(i in 1:length(time.step)){
	m2.dlwrf[i] <- mean(m.dlwrf[time.step[i]:(time.step[i]+5)])
	m2.prate[i] <- mean(m.prate[time.step[i]:(time.step[i]+5)])
	m2.pres[i] <- mean(m.pres[time.step[i]:(time.step[i]+5)])
	m2.sh[i] <- mean(m.sh[time.step[i]:(time.step[i]+5)])
	m2.tmp[i] <- mean(m.tmp[time.step[i]:(time.step[i]+5)])
	m2.ugrd[i] <- mean(m.ugrd[time.step[i]:(time.step[i]+5)])
	m2.vbdsf[i] <- mean(m.vbdsf[time.step[i]:(time.step[i]+5)])
	m2.vddsf[i] <- mean(m.vddsf[time.step[i]:(time.step[i]+5)])
	m2.nbdsf[i] <- mean(m.nbdsf[time.step[i]:(time.step[i]+5)])
	m2.nddsf[i] <- mean(m.nddsf[time.step[i]:(time.step[i]+5)])
}
length(ed.dlwrf)
length(m2.dlwrf)

plot(ed.dlwrf~m2.dlwrf[1:length(ed.dlwrf)], main="dlwrf")
	abline(0,1, col="red", lwd=2)

plot(ed.prate~m2.prate[1:length(ed.prate)], main="prate")
	abline(0,1, col="red", lwd=2)

plot(ed.pres~m2.pres[1:length(ed.pres)], main="pres")
	abline(0,1, col="red", lwd=2)

plot(ed.tmp~m2.tmp[1:length(ed.tmp)], main="tmp")
	abline(0,1, col="red", lwd=2)

plot(ed.sh~m2.sh[1:length(ed.sh)], main="actually tmp")
	abline(0,1, col="red", lwd=2)


plot(ed.ugrd~m2.ugrd[1:length(ed.ugrd)], main="ugrd")
	abline(0,1, col="red", lwd=2)

plot(ed.vbdsf~m2.vbdsf[1:length(ed.vbdsf)], main="vbdsf")
	abline(0,1, col="red", lwd=2)

radiation <- m2.vbdsf[1:length(ed.vbdsf)]+m2.vddsf[1:length(ed.vbdsf)]+m2.nbdsf[1:length(ed.vbdsf)]+m2.nddsf[1:length(ed.vbdsf)]
plot(ed.vbdsf~radiation, main="vbdsf", cex=0.1)
	abline(0,1, col="red", lwd=2)
