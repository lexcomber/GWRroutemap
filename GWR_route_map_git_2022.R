##### Code developed by Lex Comber and Paul Harris #####
##### August 2019                                  #####
##### gwr.multiscale_T.R developed by Binbin Lu	   #####
##### GWmodel now updated with this.               #####
##### Contact: a.comber@leeds.ac.uk                #####

##### The code is is organised by the 3 steps in   #####
##### the paper, with figures and tables written   #####
##### in sequence. See R_Session_Info.txt for      #####
##### details of package versions etc.             #####

##### Load packages
library(GISTools)
library(raster)
library(rgdal)
library(spdep)
library(GWmodel)
library(tidyverse)
library(repmis)
library(OpenStreetMap)
library(car)
library(broom)
require(gridExtra)
library(nlme)
# Multiscale GWR function with t-values
# update 2022 = not needed with current GWmodel 
# source("gwr.multiscale_T.R") 

##### Load in the variable data from GitHub 
source_data("https://github.com/lexcomber/GWRroutemap/blob/master/Liudaogou.RData?raw=True")
source_data("https://github.com/lexcomber/GWRroutemap/blob/master/boundary.RData?raw=True")
source_data("https://github.com/lexcomber/GWRroutemap/blob/master/MyMap.RData?raw=True")
coords <- data[,3:2]
proj. <- CRS("+proj=tmerc +lat_0=0 +lon_0=108 +k=1 +x_0=500000 +y_0=0 +ellps=krass +units=m +no_defs ")
data.sp <- SpatialPointsDataFrame(coords, data = data.frame(data), proj4string = proj.)
# transform the data
data.orig <- data
# Data tidy
data.sp@data$TNPC <- log(data.sp@data$TNPC+0.0001)
data.sp@data$TPPC <- (data.sp@data$TPPC)^0.5
data.sp@data$SOCgkg <- log(data.sp@data$SOCgkg)
data.sp@data$ClayPC <- (data.sp@data$ClayPC)^0.5
data.sp@data$NO3Ngkg <- log(abs(data.sp@data$NO3Ngkg))
data.sp@data$NH4Ngkg <- log(data.sp@data$NH4Ngkg)

##### Figure 1
crs.val <- CRS("+proj=longlat ")
# Map the data
fac = 0.001 # in case a bit is needed
tmp <- spTransform(data.sp, crs.val)
ul <- as.vector(cbind(bbox(tmp)[2,2]+fac, bbox(tmp)[1,1]-fac))
lr <- as.vector(cbind(bbox(tmp)[2,1]-fac, bbox(tmp)[1,2]+fac))
# MyMap <- openmap(ul,lr,zoom = NULL, type = "bing")
# save("MyMap", file = "MyMap.RData")
png(filename = "F1.png", w = 5, h = 6, units = "in", res = 100)
plot(MyMap, removeMargin=T) 
plot(spTransform(data.sp, osm()),add = T, col="white", pch = 19, cex =1)
scalebar(1000, type = "bar", xy = c(12286385, 4690886), divs = 4, col = "white", below = "m")
dev.off()
# watershed location
#ul+(ul - lr)/2
#ul-(ul - lr)/2

##### Step 1 Linear Regression
# Models for Analysts A to D
regmodA = as.formula(TNPC ~ SOCgkg + ClayPC + SiltPC  + NO3Ngkg + NH4Ngkg) # Analyst A = regmod2
regmodB = as.formula(TNPC ~ SandPC + NO3Ngkg) # Analyst B = regmod3
regmodC = as.formula(TNPC ~ SOCgkg + NH4Ngkg) # Analyst C = regmod4
regmodD = as.formula(TNPC ~ SOCgkg + SandPC + NO3Ngkg) # Analyst D = regmod5
# Create Linear Regresssion models
moda <- lm(regmodA, data = data.sp@data)
modb <- lm(regmodB, data = data.sp@data)
modc <- lm(regmodC, data = data.sp@data)
modd <- lm(regmodD, data = data.sp@data)

##### Table 2: OLS summaries
taba <- round(summary(moda)$coefficients[, c(1,4)], 4)
tabb <- round(summary(modb)$coefficients[, c(1,4)], 4)
tabc <- round(summary(modc)$coefficients[, c(1,4)], 4)
tabd <- round(summary(modd)$coefficients[, c(1,4)], 4)
tab2.out <- matrix(ncol = 8, nrow = 7)
# create a vector of the factors
vars = c("(Intercept)", "SOCgkg", "ClayPC", "SiltPC", "SandPC", "NO3Ngkg", "NH4Ngkg")
rownames(tab2.out) <- vars
tab2.out[match(rownames(taba), vars),1:2] <- taba
tab2.out[match(rownames(tabb), vars),3:4] <- tabb
tab2.out[match(rownames(tabc), vars),5:6] <- tabc
tab2.out[match(rownames(tabd), vars),7:8] <- tabd
colnames(tab2.out) = paste0(paste0("Analyst", c(1,1,2,2,3,3,4,4), c("Est", "p-val")))
write.csv(tab2.out, file = "tab2.csv")

##### Table 3. Moran's I and Fitness tests
# Moran's I test of regression residuals
# Autocorrelation tests on response and LM error data
# needs some kind of adjacency / neighbour network
# used Voronoi approach to define areas around points from which a neighbourhood can be computed
# http://carsonfarmer.com/2009/09/voronoi-polygons-with-r/
voronoipolygons = function(layer) {
    require(deldir)
    crds = layer@coords
    z = deldir(crds[,1], crds[,2])
    w = tile.list(z)
    polys = vector(mode='list', length=length(w))
    require(sp)
    for (i in seq(along=polys)) {
        pcrds = cbind(w[[i]]$x, w[[i]]$y)
        pcrds = rbind(pcrds, pcrds[1,])
        polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
    }
    SP = SpatialPolygons(polys)
    row.names=sapply(slot(SP, 'polygons'), function(x) slot(x, 'ID'))
    voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1], 
        y=crds[,2], row.names=row.names))
}
pt.zone <- voronoipolygons(data.sp)
# check with a simple plot
# plot(pt.zone)
# plot(data.sp, add = T, pch = 19, cex = 0.5)
# Now determine the neighbourhood list
gnb = poly2nb(pt.zone)
glw = nb2listw(gnb)
# apply Moran Test (UNBIASED)
mora <- lm.morantest(moda, listw=glw)
morb <- lm.morantest(modb, listw=glw)
morc <- lm.morantest(modc, listw=glw)
mord <- lm.morantest(modd, listw=glw)
# Create Table 3 output - Part 1 - Morans's I and p-value
tab3.out <- matrix(ncol = 2, nrow = 4)
rownames(tab3.out) <- paste0("Analyst ", c("A","B","C","D"))
tab3.out[1, ] <- c(unlist(mora[3])[1],unlist(mora[2]) )
tab3.out[2, ] <- c(unlist(morb[3])[1],unlist(morb[2]) )
tab3.out[3, ] <- c(unlist(morc[3])[1],unlist(morc[2]) )
tab3.out[4, ] <- c(unlist(mord[3])[1],unlist(mord[2]) )
# Fitness tests - to get AIC and R2
# use GWR with a boxcar and max distance to get standard forms of AICc and R2
# ie AICc calulated in the same way acorss models
# essentially GWR is used here to created global linear models
EUDM <- gw.dist(coordinates(data.sp))
bw = max(EUDM)
moda2 <- gwr.basic(regmodA, data = data.sp, bw = bw, kernel = "boxcar", dMat = EUDM)
modb2 <- gwr.basic(regmodB, data = data.sp, bw = bw, kernel = "boxcar", dMat = EUDM)
modc2 <- gwr.basic(regmodC, data = data.sp, bw = bw, kernel = "boxcar", dMat = EUDM)
modd2 <- gwr.basic(regmodD, data = data.sp, bw = bw, kernel = "boxcar", dMat = EUDM)
# extract fitness measures
moda.fit <- c(moda2$GW.diagnostic$AICc, moda2$GW.diagnostic$gw.R2)
modb.fit <- c(modb2$GW.diagnostic$AICc, modb2$GW.diagnostic$gw.R2)
modc.fit <- c(modc2$GW.diagnostic$AICc, modc2$GW.diagnostic$gw.R2)
modd.fit <- c(modd2$GW.diagnostic$AICc, modd2$GW.diagnostic$gw.R2)
# Create Table 3 output - Part 2 - AIC and R2
tab3.out.b <- matrix(nrow = 4, ncol = 2)
tab3.out.b[1,] <- moda.fit
tab3.out.b[2,] <- modb.fit
tab3.out.b[3,] <- modc.fit
tab3.out.b[4,] <- modd.fit
tab3.out <- cbind(tab3.out, tab3.out.b)
colnames(tab3.out) = c("Moran’s I",	"p-value", "AICc", "R2")
write.csv(round(tab3.out, 3), file = "tab3.csv")

##### Step 2 MS-GWR - Multiscale GWR

# Here first the bandwidths are found using scaled predictor data
# Then MS-GWR models are re-fit with specified / estimated bandwidths to unscaled predictor data

# The reasons are complex but are summarised in Oshan et al (2019) and the response in Lu et al (2019)
# Taylor Oshan, Levi John Wolf, A. Stewart Fotheringham, Wei Kang, Ziqi Li & Hanchen Yu (2019) A comment on geographically weighted regression with parameter-specific distance metrics, International Journal of Geographical Information Science, 33:7, 1289-1299
# Binbin Lu, Chris Brunsdon, Martin Charlton & Paul Harris (2019) A response to ‘A comment on geographically weighted regression with parameter-specific distance metrics’, International Journal of Geographical Information Science, 33:7, 1300-1312

# Note that Initial bandwidth values (bws0) are passed to the function 
# these are approximately 1/4 to 1/3 of the distance (1250m) or number of data points (200)

# Multiscale GWR takes time especially if the predictor data are not scaled to mean zero

# The code below runs the MS-GWR and saves the results locally - this takes time 
# When you have run it you can change the flag below to FALSE 
# This will then pick up the saved files 

do.flag = TRUE
if(do.flag) {
	gw.ma <- gwr.multiscale(regmodA, data = data.sp, 
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(1250,1250,1250,1250,1250,1250),
	                        #dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
	                        verbose = F, predictor.centered=rep(T, 5))
	# Adaptive version of Analyst A
	gw.ma2 <- gwr.multiscale(regmodA, data = data.sp, adaptive = T,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(200,200,200,200,200,200),
	                        #dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
	                        verbose = F, predictor.centered=rep(T, 5))
	gw.mb <- gwr.multiscale(regmodB, data = data.sp,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(1250,1250,1250),
	                        #dMats=list(EUDM,EUDM,EUDM),
	                        verbose = F, predictor.centered=rep(T, 2))
	gw.mc <- gwr.multiscale(regmodC, data = data.sp, 
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(1250,1250,1250),
	                        #dMats=list(EUDM,EUDM,EUDM),
	                        verbose = F, predictor.centered=rep(T, 2))
	gw.md <- gwr.multiscale(regmodD, data = data.sp, 
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(1250,1250,1250,1250),
	                        #dMats=list(EUDM,EUDM,EUDM,EUDM),
	                        verbose = F, predictor.centered=rep(T, 3))
	save(list = c("gw.ma", "gw.ma2", "gw.mb", "gw.mc", "gw.md"), 
		file = "MGWR_part1.RData")
}
if(!do.flag){
	load("MGWR_part1.RData")
}
# END do.flag

##### Table 4 - the table of MS-GWR Bandwidths
tab <- matrix(nrow = 5, ncol = 7)
colnames(tab) <- c("Intercept","SOCgkg","ClayPC","SiltPC","SandPC","NO3Ngkg","NH4Ngkg")
tab[1, -5] <- round(gw.ma[[2]]$bws,1)
tab[2, -5] <- round(gw.ma2[[2]]$bws,0)
tab[3, c(1,5:6)] <- round(gw.mb[[2]]$bws,1)
tab[4, c(1,2,7)] <- round(gw.mc[[2]]$bws,1)
tab[5, c(1,2,5,6)] <- round(gw.md[[2]]$bws,1)
write.csv(data.frame(tab), file = "tab4.csv")

# Now re-run multiscale GWR with bandwidths found from above but with predictors NOT centered
# assign bandwidths 
mbwa <- round(gw.ma[[2]]$bws,1)
mbwa2 <- round(gw.ma2[[2]]$bws,0)
mbwb <- round(gw.mb[[2]]$bws,1)
mbwc <- round(gw.mc[[2]]$bws,1)
mbwd <- round(gw.md[[2]]$bws,1)
# change flag as before 
do.flag = TRUE
if(do.flag) {
  gw.ma_2 <- gwr.multiscale(regmodA, data = data.sp,
                          max.iterations = 1000,
                          criterion="CVR", kernel = "bisquare", 
                          bws0=c(mbwa),
                          bw.seled=rep(T, 6),
                          #dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
                          verbose = F, predictor.centered=rep(F, 5))
  # the adaptive version Analyst A is never used again
  # so not run here
  gw.mb_2 <- gwr.multiscale(regmodB, data = data.sp,
                          max.iterations = 1000,
                          criterion="CVR", kernel = "bisquare", 
                          bws0=c(mbwb),
                          bw.seled=rep(T, 3),
                          #dMats=list(EUDM,EUDM,EUDM),
                          verbose = F, predictor.centered=rep(F, 2))
  gw.mc_2 <- gwr.multiscale(regmodC, data = data.sp, 
                          max.iterations = 1000,
                          criterion="CVR", kernel = "bisquare", 
                          bws0=c(mbwc),
                          bw.seled=rep(T, 3),
                          #dMats=list(EUDM,EUDM,EUDM),
                          verbose = F, predictor.centered=rep(F, 2))
  gw.md_2 <- gwr.multiscale(regmodD, data = data.sp, 
                          max.iterations = 1000,
                          criterion="CVR", kernel = "bisquare", 
                          bws0=c(mbwd),
                          bw.seled=rep(T, 4),
                          #dMats=list(EUDM,EUDM,EUDM,EUDM),
                          verbose = F, predictor.centered=rep(F, 3))
  save(list = c("gw.ma_2", "gw.mb_2", "gw.mc_2", "gw.md_2"), 
       file = "MGWR_part2.RData")
}
# END do.flag
if(!do.flag){
  load("MGWR_part2.RData")
}

# Bandwidths check
round(gw.ma_2[[2]]$bws,1)
round(gw.mb_2[[2]]$bws,1)
round(gw.mc_2[[2]]$bws,1)
round(gw.md_2[[2]]$bws,1)

##### Table 5: MS-GWR Moran's I (BIASED), AIC R2 etc
# Moran's i
gw.ma.i = moran.test(gw.ma_2$SDF$residual, listw=glw)
gw.mb.i = moran.test(gw.mb_2$SDF$residual, listw=glw)
gw.mc.i = moran.test(gw.mc_2$SDF$residual, listw=glw)
gw.md.i = moran.test(gw.md_2$SDF$residual, listw=glw)
tab5.out <- matrix(ncol = 2, nrow = 4)
rownames(tab5.out) <- rownames(tab3.out)
tab5.out[1, ] <- c(unlist(gw.ma.i[3])[1],unlist(gw.ma.i[2]) )
tab5.out[2, ] <- c(unlist(gw.mb.i[3])[1],unlist(gw.mb.i[2]) )
tab5.out[3, ] <- c(unlist(gw.mc.i[3])[1],unlist(gw.mc.i[2]) )
tab5.out[4, ] <- c(unlist(gw.md.i[3])[1],unlist(gw.md.i[2]) )
# R2 and AIC
moda.fit <- c(gw.ma_2$GW.diagnostic$AICc, gw.ma_2$GW.diagnostic$R2.val)
modb.fit <- c(gw.mb_2$GW.diagnostic$AICc, gw.mb_2$GW.diagnostic$R2.val)
modc.fit <- c(gw.mc_2$GW.diagnostic$AICc, gw.mc_2$GW.diagnostic$R2.val)
modd.fit <- c(gw.md_2$GW.diagnostic$AICc, gw.md_2$GW.diagnostic$R2.val)
# create a table to link to 
tab5.out.b <- matrix(nrow = 4, ncol = 2)
tab5.out.b[1,] <- moda.fit
tab5.out.b[2,] <- modb.fit
tab5.out.b[3,] <- modc.fit
tab5.out.b[4,] <- modd.fit
tab5.out <- cbind(tab5.out, tab5.out.b)
colnames(tab5.out) = c("Moran’s I",	"p-value", "AICc", "R2")
write.csv(round(tab5.out, 3), file = "tab5.csv")

##### Step 3 Final Models
# format boundary for mapping with ggplot
b <- gBuffer(boundary, width = 25)
b <- tidy(b)
# define some mapping functions
# Map function for GWR coefficients that are either entirely positive or negative
plot.func <- function(df = df, pal = "YlOrRd", dir = 1,
						var = df$SandPC, tval = df$SandPC_TV, tit = "Sand %") {
	index = tval < -1.96 | tval > 1.96
	df.index = df[index,]
	p = ggplot() + 
		geom_polygon(data = b, aes(x = long, y = lat), fill = "lightgrey") + 
		geom_point(data = df, aes(x = Longitude, y = Latitude, colour = var), size = 2.5) +
		 		 labs(x = NULL, y = NULL)+
		 		 coord_fixed() +	
		 		 scale_colour_distiller(type = "seq", palette = pal, name = tit, 
		 		 						direction = dir, na.value = "lightgrey",
		 		 						#breaks = seq(0, 1, 0.2), limits=c(0, 1),
		 		 						guide = guide_colourbar(title.position = "top",
		 		 					    		title.hjust = 0.5,
		 		 					        title.vjust = 0.8,
		 		 					        barwidth = 10)) +
		geom_point(data = df.index, aes(x = Longitude, y = Latitude), colour = "black", 
			shape = 1, size = 3.5)+
	 		 theme(axis.line = element_blank(),
		  		 	   axis.text = element_blank(),
		  		 	   axis.ticks = element_blank(),
		  		 	   legend.position = "bottom",
		  		 	   legend.direction = "horizontal",
		  		 	   panel.grid.major = element_blank(),
		  		 	   panel.grid.minor = element_blank(),
		  		 	   panel.border = element_blank(),
		  			   panel.background = element_blank())
	return(p)
}
# Map function for GWR coefficients that go from negative to positive 
plot.func.div <- function(df = df, var = df$SandPC,tval = df$SandPC_TV, tit = "Sand %") {
	index = tval < -1.96 | tval > 1.96
	df.index = df[index,]
	p = ggplot() + 
		geom_polygon(data = b, aes(x = long, y = lat), fill = "lightgrey") + 
		geom_point(data = df, aes(x = Longitude, y = Latitude, colour = var), size = 2.5) +
		 		 labs(x = NULL, y = NULL)+
		 		 coord_fixed() +	
	    scale_colour_gradient2(low = "#CB181D",
                         mid = "white",
                         high = "#2171B5",
                         midpoint = 0, name = tit,
                         guide = guide_colourbar(title.position = "top",
		 		 					    		title.hjust = 0.5,
		 		 					        title.vjust = 0.8,
		 		 					        barwidth = 10)) +				    
		geom_point(data = df.index, aes(x = Longitude, y = Latitude), colour = "black", 
			shape = 1, size = 3.5)+
		 		 theme(axis.line = element_blank(),
		  		 	   axis.text = element_blank(),
		  		 	   axis.ticks = element_blank(),
		  		 	   legend.position = "bottom",
		  		 	   legend.direction = "horizontal",
		  		 	   panel.grid.major = element_blank(),
		  		 	   panel.grid.minor = element_blank(),
		  		 	   panel.border = element_blank(),
		  			   panel.background = element_blank())
	return(p)
}

##### Fixed coefficient Linear Regression and SAM with REML - Analyst C
# reminder of the LR
summary(modc)
tab <- round(summary(modc)$coefficients, 4)
# REML as as SAM
data.sp$dummy <- rep(1, nrow(data.sp@data))
reml.mod <- lme(fixed = regmodC, data = data.sp@data, random = ~ 1 | dummy, method = "REML")
reml.up <- update(reml.mod, correlation = corExp(c(1000,0.03),
	form = ~ Longitude+Latitude, nugget = T), 
	method = "REML")
AIC(modc, reml.up)
pv = anova(reml.up)[,4]
cf = reml.up$coefficients$fixed
# Table 6:
tab6_out = cbind(tab[,c(1,4)], cbind(cf, pv))
write.csv(round(tab6_out, 3), "tab6.csv")

##### Mixed and Multiscale GWR - Analyst A
# Here using multiscale GWR and fixing the bandwidths to mimic a mixed GWR
# again change flag if already run
do.flag = TRUE
if(do.flag) {
gw.mix <- gwr.multiscale(regmodA, data = data.sp,
	                      max.iterations = 1000,
	                      criterion="CVR", kernel = "bisquare", 
	                      bws0=c(700, 100000,100000,700,700,100000),
                          bw.seled=rep(T, 6),
	                      #dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
                          verbose = F, predictor.centered=rep(F, 5))
save(gw.mix, file = "MGWR_part3.RData")
}
# END do.flag
if(!do.flag){
	load("MGWR_part3.RData")
}
# check the AICc
gw.mix$GW.diagnostic$AICc
# And map the results 
# here it is helpful to determine the proprties of the coefficients
# and then to select the plot function and shading schemes
summary(gw.mix$SDF)
# Coefficient Maps of MX-GWR as the selected approach
df <- as.data.frame(gw.mix$SDF)
p1 = plot.func(df = df, pal = "YlOrRd", dir = -1,
						var = df$Intercept, tval = df$Intercept_TV, tit = "Intercept")
p2 = plot.func(df = df, pal = "GnBu", dir = -1,
               var = df$SiltPC, tval = df$SiltPC_TV, tit = "Silt %") 
p3 <- plot.func.div(df = df, var = df$NO3Ngkg, tval = df$NO3Ngkg_TV, 
		tit = expression(paste("NO"[3],"N")))
##### Figure 2
png(filename = "F2.png", w = 12, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()
# Coefficient Maps for comparison with Multiscale GWR maps
gw.ma_2$GW.diagnostic$AICc
gw.mix$GW.diagnostic$AICc
# Maps of MS-GWR for comparison
df <- as.data.frame(gw.ma_2$SDF)
summary(df)
p1 = plot.func(df = df, pal = "YlOrRd", dir = -1,
               var = df$Intercept, tval = df$Intercept_TV, tit = "Intercept") 
p2 = plot.func(df = df, pal = "GnBu", dir = -1,
               var = df$SOCgkg, tval = df$SOCgkg_TV, tit = "SOC")
p3 = plot.func(df = df, pal = "YlOrRd", dir = -1,
               var = df$ClayPC, tval = df$ClayPC_TV, tit = "Clay %")
p4 = plot.func(df = df, pal = "GnBu", dir = -1,
               var = df$SiltPC, tval = df$SiltPC_TV, tit = "Silt %")
p5 = plot.func.div(df = df, var = df$NO3Ngkg, tval = df$NO3Ngkg_TV, tit = expression(paste("NO"[3],"N")))
p6 = plot.func(df = df, pal = "YlOrRd", dir = -1,
               var = df$NH4Ngkg, tval = df$NH4Ngkg_TV, tit = expression(paste("NH"[4],"N")))
##### Figure 3
png(filename = "F3.png", w = 12, h = 12, units = "in", res = 100)
grid.arrange(p1, p4, p5, p2, p3, p6, ncol=3, nrow=2)
dev.off()

##### Multiscale GWR - Analyst D
gw.md_2
gw.md_2$GW.diagnostic$AICc
# Coefficient Maps of MS-GWR as the selected approach
df <- as.data.frame(gw.md_2$SDF)
summary(df)
p1 = plot.func(df = df, pal = "YlOrRd", dir = -1,
						var = df$Intercept, tval = df$Intercept_TV, tit = "Intercept") 
p2 = plot.func(df = df, pal = "GnBu", dir = -1,
						var = df$SOCgkg, tval = df$SOCgkg_TV, tit = "SOC") 
p3 = plot.func(df = df, pal = "YlOrRd", dir = -1,
						var = df$SandPC, tval = df$SandPC_TV, tit = "Sand %") 
p4 = plot.func.div(df = df, var = df$NO3Ngkg, tval = df$NO3Ngkg_TV, tit = expression(paste("NO"[3],"N")))
##### Figure 4
png(filename = "F4.png", w = 16, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, p4, ncol=4)
dev.off()

##### Standard and Multiscale GWR - Analyst B
bw = bw.gwr(regmodB, data = data.sp, approach="AICc",kernel="bisquare",
                 adaptive=F, p=2, theta=0, longlat=F,dMat = EUDM)
gw = gwr.basic(regmodB, data = data.sp, bw = bw, kernel="bisquare", cv = F,
                 adaptive=F, p=2, theta=0, longlat=F,dMat = EUDM)
# Bandwidth function - evaluates different bandwidths
# GWR - AIC approach - provides a useful visual check on its behaviour.....
n.min <- 200 # user-specified minimum bandwidth
n.max <- max(EUDM)+100 # user-specified maximum bandwidth
interval.size <- 20
fixed <- seq(n.min,n.max,by=interval.size)
b.func.gwr <- matrix(nrow=length(fixed),ncol=1)
for(i in 1:length(fixed)) {
	g.gwr <- gwr.aic(fixed[i], Y = as.matrix(data.sp@data$TNPC), 
							  X = as.matrix(data.sp@data[, c("SandPC","NO3Ngkg")]), 
							  kernel="bisquare", adaptive=F, 
							  dp.locat = coordinates(data.sp), dMat=EUDM)
	b.func.gwr[i] <- g.gwr[1]
	if(i%%10 ==0) cat(i, "\t")
}
fixed[which.min(b.func.gwr)]
xy <- data.frame(x = fixed,y = b.func.gwr)
##### Figure 5
png(filename = "F5.png", w = 9, h = 3.5, units = "in", res = 100)
ggplot() + 
    geom_point(data = xy, aes(x=x, y=y), size = 0.7, alpha = 0.5) +
    geom_line(data = xy, aes(x=x, y=y)) +
    labs(
    	#subtitle = "GWR bandwidth function", 
    	x = "Bandwidth size (m)", 
    	y = "AIC")
dev.off()
# Coefficient Maps of the final GWR as the selected model
gw$GW.diagnostic$AICc
df <- as.data.frame(gw$SDF)
summary(df)
p1 = plot.func.div(df = df, var = df$Intercept, tval = df$Intercept_TV, tit = "Intercept") 
p2 = plot.func.div(df = df, var = df$SandPC, tval = df$SandPC_TV, tit = "Sand %")
p3 = plot.func(df = df, pal = "GnBu", dir = -1, var = df$NO3Ngkg, tval = df$NO3Ngkg_TV, tit = expression(paste("NO"[3],"N")))
##### Figure 6
png(filename = "F6.png", w = 12, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()
# Coefficient Maps of the final Multiscale GWR as the alternative selected model
df <- as.data.frame(gw.mb_2$SDF)
summary(df)
p1 = plot.func(df = df, pal = "YlOrRd", dir = -1,
               var = df$Intercept, tval = df$Intercept_TV, tit = "Intercept")  
p2 = plot.func(df = df, pal = "YlOrRd", dir = -1,
               var = df$SandPC, tval = df$SandPC_TV, tit = "Sand %") 
p3 = plot.func(df = df, pal = "GnBu", dir = -1,
               var = df$NO3Ngkg, tval = df$NO3Ngkg_TV, tit = expression(paste("NO"[3],"N")))
##### Figure 7
png(filename = "F7.png", w = 12, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

##### AIC table - Table 7 
tab.aic = data.frame("OLS" = c(AIC(moda), AIC(modb), AIC(modc), AIC(modd)),
          "MS-GWR" = c(gw.ma_2$GW.diagnostic$AICc,gw.mb_2$GW.diagnostic$AICc,
          	gw.mc_2$GW.diagnostic$AICc, gw.md_2$GW.diagnostic$AICc),
          "Chosen" = c(gw.mix$GW.diagnostic$AICc, gw.mb_2$GW.diagnostic$AICc, 
          	AIC(modc), gw.md_2$GW.diagnostic$AICc))
tab.aic
write.csv(round(tab.aic, 3), "tab7.csv")

##### Save image
save.image("GWR_route_map_ALL.RData")

##### END #####