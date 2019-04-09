
### Load data and libraries
### set up data.sp etc
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

### load in the variable data

# load data from GitHub 
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

#### Step 1 OLS

terms <- names(data)[c(4:11)]
## OLS N percentage
regmod <- paste(terms[1], "~")
for ( i in 3: length(terms) ) {
	if ( i != length(terms) ) regmod <- paste(regmod, terms[i],  "+")
	if ( i == length(terms) ) regmod <- paste(regmod, terms[i])
}
regmod.n <- as.formula(regmod)
summary(lm (regmod.n, data.sp@data))

## OLS P percentage
regmod <- paste(terms[2], "~")
for ( i in 3: length(terms) ) {
	if ( i != length(terms) ) regmod <- paste(regmod, terms[i],  "+")
	if ( i == length(terms) ) regmod <- paste(regmod, terms[i])
}
regmod.p <- as.formula(regmod)
summary(lm (regmod.p, data.sp@data))

mN <- lm(regmod.n, data = data.sp@data)
mP <- lm(regmod.p, data = data.sp@data)
tab1 <- round(summary(mN)$coefficients, 4)
tab2 <- round(summary(mP)$coefficients, 4)

write.csv(cbind(tab1, tab2), file = "tab1.csv")

#### Moran's I test of regression residuals
# Autocorrelation tests on response and LM error data
# moran.resid.i <- lm.morantest(ols.model.i, listw=glw) # regression residuals
# moran.pval.i <- moran.resid.i[[2]]
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
plot(pt.zone)
vor.mask <- poly.outer(pt.zone,boundary,extend=100)

# Moran Test
gnb = poly2nb(pt.zone)
glw = nb2listw(gnb)
moran.resid.n <- lm.morantest(mN, listw=glw) 
moran.resid.n
moran.resid.p <- lm.morantest(mP, listw=glw) 
moran.resid.p

# Local Moran
lI <- localmoran(mN$residuals,glw)
lI[is.na(lI)] <- 0
# Mappings 
par(mar = c(0,0,1,0))
sids.shade <- auto.shading(mN$residuals, cols=brewer.pal(5,"PRGn"))
choropleth(pt.zone,mN$residuals,shading=sids.shade, border = NA)
add.masking(vor.mask)
choro.legend(c("topright"),sh=sids.shade,fmt="%5.2f", cex = 0.7, box.lwd = 0) # Add
title("Nitrogen", cex.main=1) # Add title

#### VIFs
t(round(vif(lm(mN)),3))
r1 <- (round(vif(lm(mP)),3))
r2 <- rep(0, length(r1))
r3 <- vif(lm(TPPC ~ SOCgkg + ClayPC + SiltPC + NO3Ngkg + NH4Ngkg, data = data.sp@data))
r2 <- r3[match(names(r1), names(r3))] 
rbind(r1, r2)

write.csv(rbind(r1, r2), file = "tab2.csv")

#### Residuals - FIGURE 2
## define plot functions
# define areas
b <- gBuffer(boundary, width = 25)
b <- tidy(b)
# resid plot with legend
resid.plot = function(n.cols, cols, tit = "TSN Residuals") {
	Residuals = factor(n.cols, levels = c("< -2 sd","-2 < sd > +2","> +2 sd"))
	p = ggplot() + 
		geom_polygon(data = b, aes(x = long, y = lat), fill = "lightgrey") + 
		geom_point(data = data.sp@data, aes(x = Longitude, y = Latitude, 
			colour = Residuals), 
			size = 2.5) +
		labs(x = NULL, y = NULL)+
		coord_fixed() +	
		scale_colour_manual(values = cols, 
			guide = guide_legend(title.position = "top"),
			name = tit)+
		theme(axis.line = element_blank(),
		  		 	   axis.text = element_blank(),
		  		 	   axis.ticks = element_blank(),
		  		 	   legend.position = "bottom",
		  		 	   legend.direction = "horizontal",
		  		 	   panel.grid.major = element_blank(),
		  		 	   panel.grid.minor = element_blank(),
		  		 	   panel.border = element_blank(),
		  			   panel.background = element_blank(),
		  			   legend.title.align=0.5)
	return(p)
}
# local moran plot with legend
plot.moran = function(morani = morani.n[,1], tit = "Local Moran's I (TSN)") {
	p = ggplot() + 
		geom_polygon(data = b, aes(x = long, y = lat), fill = "lightgrey") + 
		geom_point(data = data.sp@data, aes(x = Longitude, y = Latitude, col = morani),
			size = 2.5)+
		labs(x = NULL, y = NULL)+
		coord_fixed() +		
		scale_colour_gradient2(low = brewer.pal(5, "PRGn")[5], high = brewer.pal(5, "PRGn")[1],
			name = tit, na.value = "lightgrey",
			guide = guide_colourbar(title.position = "top",
			title.hjust = 0.5, title.vjust = 0.8, barwidth = 10)) +
		theme(axis.line = element_blank(),
		  		 	   axis.text = element_blank(),
		  		 	   axis.ticks = element_blank(),
		  		 	   legend.position = "bottom",
		  		 	   legend.direction = "horizontal",
		  		 	   panel.grid.major = element_blank(),
		  		 	   panel.grid.minor = element_blank(),
		  		 	   panel.border = element_blank(),
		  			   panel.background = element_blank(),
		  			   legend.title.align=0.5)
	return(p)
}
## calulate standardized residuals 
summary(mN$residuals)
n.resids = rstandard(mN)summary(n.resids)
p.resids = rstandard(mP)
summary(p.resids)

resid.shades = shading(c(-1.96,1.96),c("< -2 sd","-2 < sd > +2","> +2 sd"))
n.cols = resid.shades$cols[1 + findInterval(n.resids, resid.shades$breaks)]
p.cols = resid.shades$cols[1 + findInterval(p.resids, resid.shades$breaks)]
cols = c("red", "darkgrey", "blue")
# plot N
p1 = resid.plot(n.cols, cols, "TSN Standardized Residuals") 
# plot P
p3 = resid.plot(p.cols, cols, "TSP Standardized Residuals") 

## calculate local moran's I of residuals
morani.n <- localmoran(n.resids,glw)
morani.p <- localmoran(p.resids,glw)
p2 = plot.moran(morani.n[,1], tit = "TSN Local Moran's I") 
p4 = plot.moran(morani.p[,1], tit = "TSP Local Moran's I") 

png(filename = "F2.png", w = 18, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, p4, ncol=4)
dev.off()

#### Step 2 MS-GWR
regmod.n
regmod.p
regmod.n2 <- as.formula(TNPC ~ SOCgkg + ClayPC + SiltPC + NO3Ngkg + NH4Ngkg)
regmod.n2 <- as.formula(TNPC ~ SOCgkg + ClayPC + SiltPC + NO3Ngkg + NH4Ngkg)
regmod.n3 <- as.formula(TNPC ~ SandPC + NO3Ngkg)
regmod.n4 <- as.formula(TNPC ~ SOCgkg + ClayPC + NH4Ngkg) 
regmod.n5 <- as.formula(TNPC ~ SOCgkg + SandPC + NO3Ngkg) 
regmod.n6 <- as.formula(TNPC ~ NO3Ngkg)

EUDM <- gw.dist(coordinates(data.sp))

# this takes time
# change the flag below to TRUE to run and save the results
# then change back to FALSE and it will pick up the saved file
do.flag = FALSE
if(do.flag) {
	gw.m1 <- gwr.multiscale(regmod.n, data = data.sp,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m1a <- gwr.multiscale(regmod.n, data = data.sp, adaptive = T,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m2 <- gwr.multiscale(regmod.n2, data = data.sp,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m2a <- gwr.multiscale(regmod.n2, data = data.sp, adaptive = T,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m3 <- gwr.multiscale(regmod.n3, data = data.sp,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m3a <- gwr.multiscale(regmod.n3, data = data.sp, adaptive = T,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m4 <- gwr.multiscale(regmod.n4, data = data.sp, 
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m4a <- gwr.multiscale(regmod.n4, data = data.sp, adaptive = T,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m5 <- gwr.multiscale(regmod.n5, data = data.sp, 
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
	gw.m5a <- gwr.multiscale(regmod.n5, data = data.sp, adaptive = T,
	                        max.iterations = 1000,
	                        criterion="CVR", kernel = "bisquare", 
	                        bws0=c(50,50,50,50),
	                        dMats=list(EUDM,EUDM,EUDM,EUDM),
	                        verbose = F)
save.image("rm_v1_part1.RData")
}
# END do.flag
if(!doFlag){
	load("rm_v1_part1.RData")
}
tab <- matrix(nrow = 10, ncol = 7)
colnames(tab) <- c("Intercept","SOCgkg","ClayPC","SiltPC","SandPC","NO3Ngkg","NH4Ngkg")
tab[1,] <- round(gw.m1[[2]]$bws,1)
tab[2, ] <- round(gw.m1a[[2]]$bws,0)
tab[3, -5] <- round(gw.m2[[2]]$bws,1)
tab[4, -5] <- round(gw.m2a[[2]]$bws,1)
tab[5, c(1,5:6)] <- round(gw.m3[[2]]$bws,1)
tab[6, c(1,5:6)] <- round(gw.m3a[[2]]$bws,1)
tab[7, c(1,2,3,7)] <- round(gw.m4[[2]]$bws,1)
tab[8, c(1,2,3,7)] <- round(gw.m4a[[2]]$bws,1)
tab[9, c(1,2,5,6)] <- round(gw.m5[[2]]$bws,1)
tab[10,c(1,2,5,6)] <- round(gw.m5a[[2]]$bws,1)
tab <- data.frame(tab)

write.csv(tab, file = "tab3.csv")

#### Step 3 Final Model
plot.func <- function(df = df, pal = "YlOrRd", dir = 1,
						var = df$SandPC_L, tit = "Sand %") {
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
plot.func.div <- function(df = df,var = df$SandPC_L, tit = "Sand %") {
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
### OLS 
regmod.n4
ols.mod <- lm(regmod.n4, data = data.sp@data)
summary(ols.mod)
tab <- round(summary(ols.mod)$coefficients, 4)

write.csv(tab, "tab4.csv")

## REML
data.sp$dummy <- rep(1, nrow(data.sp@data))
reml.mod <- lme(fixed = regmod.n4, data = data.sp@data, 
				random = ~ 1 | dummy, method = "REML")
reml.up <- update(reml.mod, correlation = corExp(form = ~Longitude+Latitude), method = "REML")
names(reml.up)$coefficients
coef(reml.up)
AIC(ols.mod, reml.up)

MuMIn::r.squaredGLMM(reml.up)

### Mixed GWR
# Model
# fix intercept series of 1s
gw.mix <- gwr.mixed(regmod.n, data = data.sp, 
					fixed.vars = c("SOCgkg", "ClayPC", "SiltPC", "NH4Ngkg"),
                    intercept.fixed=FALSE, bw = 700, diagnostic=T, kernel="bisquare",
                    adaptive=F, p=2, theta=0, longlat=F,dMat =EUDM)
# coeffs
tab <- apply(gw.mix$SDF@data, 2, summary)

write.csv(tab, "tab5.csv")
(gw.mix$aic)
(gw.mix$r.ss)

# maps
df <- as.data.frame(gw.mix$SDF)
p1 <- plot.func(df = df, pal = "GnBu",dir = -1, var = df$Intercept_L, tit = "Intercept")
p2 <- plot.func.div(df = df, var = df$SandPC_L, tit = "Sand %")
p3 <- plot.func.div(df = df, var = df$NO3Ngkg_L, tit = expression(paste("NO"[3],"N")))

##### Figure 3
png(filename = "F3.png", w = 15, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

### Multiscale GWR
# Model
gw.m5a
# coeffs
tab <- apply(gw.m5a$SDF@data, 2, summary)
tab <- rbind(tab[,1:4], gw.m5a$GW.arguments$bws)

write.csv(tab, "tab6.csv")
gw.m5a$GW.diagnostic$AICc

# maps
df <- as.data.frame(gw.m5a$SDF)
summary(df)
p1 <- plot.func(df = df, pal = "GnBu",dir = -1, var = df$Intercept, tit = "Intercept")
p2 <- plot.func(df = df, pal = "YlOrRd",var = df$SOCgkg, tit = "SOC")
p3 <- plot.func(df = df,var = df$SandPC, pal = "GnBu",dir = -1, tit = "Sand %")
p4 <- plot.func.div(df = df,var = df$NO3Ngkg, tit = expression(paste("NO"[3],"N")))

##### Figure 4

png(filename = "F4.png", w = 18, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, p4, ncol=4)
dev.off()

### Ordinary GWR
# Model
bw = bw.gwr(regmod.n3, data = data.sp, approach="CV",kernel="bisquare",
                 adaptive=F, p=2, theta=0, longlat=F,dMat = EUDM)

gw = gwr.basic(regmod.n3, data = data.sp, bw = bw, kernel="bisquare", cv = F,
                 adaptive=T, p=2, theta=0, longlat=F,dMat = EUDM)
# Bandwidth function
# GWR - CV approach only - provides a useful visual check on its behaviour.....
n.min <- 50 # user-specified minimum bandwidth
n.max <- max(EUDM) # user-specified maximum bandwidth
interval.size <- 20
fixed <- seq(n.min,n.max,by=interval.size)
b.func.gwr <- matrix(nrow=length(fixed),ncol=1)
for(i in 1:length(fixed)) {
	g.gwr <- gwr.cv(fixed[i], Y = as.matrix(data.sp@data$TNPC), 
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
    geom_vline(xintercept = bw, colour = "red") +
    scale_x_continuous(breaks = seq(100, max(EUDM), 200))+
    labs(
    	#subtitle = "GWR bandwidth function", 
    	x = "Kernel size (m)", 
    	y = "CV score")
dev.off()               
# coeffs
tab <- apply(gw$SDF@data, 2, summary)

write.csv(tab[,1:3], "tab7.csv")
gw$GW.diagnostic$AICc

# maps
df <- as.data.frame(gw$SDF)
summary(df)
p1 <- plot.func.div(df = df, var = df$Intercept, tit = "Intercept")
p2 <- plot.func(df = df, pal = "GnBu",dir = -1, var = df$SandPC, tit = "Sand %")
p3 <- plot.func(df = df, pal = "YlOrRd",var = df$NO3Ngkg, 
		tit = expression(paste("NO"[3],"N")))

##### Figure 6
png(filename = "F6.png", w = 15, h = 6, units = "in", res = 100)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

save.image("GWR_route_map.RData")


##### END #####